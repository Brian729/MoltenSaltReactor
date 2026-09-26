#!/usr/bin/env python3
"""Generate (and optionally run) the case matrix.

  python3 scripts/run_matrix.py generate      # writes cases/<name>/
  python3 scripts/run_matrix.py run [-j 8] [names...]   # runs Allrun in parallel (skips finished cases)
"""
import os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ITERS = 4000


def matrix():
    cases = []
    def add(name, **kw):
        cases.append((name, kw))
    powers = {0.005: [1, 3, 10, 20, 30, 40], 0.010: [1, 2, 3, 5, 10]}
    for d, qs in powers.items():
        for q in qs:
            for grav, tag in ((1, "conv"), (0, "cond")):
                add(f"d{int(d*1000):02d}mm_q{q:g}_{tag}", d_f=d, q_avg=q * 1e6, gravity=grav, **({} if grav else {"iters": 2000}))
    # coolant 600 C uniform variant
    for d, qs in ((0.005, [10, 20]), (0.010, [3, 5])):
        for q in qs:
            for grav, tag in ((1, "conv"), (0, "cond")):
                add(f"d{int(d*1000):02d}mm_q{q:g}_{tag}_Tc600", d_f=d, q_avg=q * 1e6, gravity=grav, coolant="uniform",
                    **({} if grav else {"iters": 2000}))
    # finite coolant-side heat transfer coefficient (Robin BC) variant
    for d, q in ((0.005, 10), (0.010, 3)):
        for h in (1000, 300):
            add(f"d{int(d*1000):02d}mm_q{q:g}_conv_h{h}", d_f=d, q_avg=q * 1e6, gravity=1, h=h)
    # mesh check (d = 5 mm, q = 20 MW/m3, convection)
    add("d05mm_q20_conv_meshCoarse", d_f=0.005, q_avg=20e6, gravity=1, nxf=20, nxg=4, nz=408)
    add("d05mm_q20_conv_meshFine", d_f=0.005, q_avg=20e6, gravity=1, nxf=80, nxg=16, nz=1630, iters=6000)
    return cases


def generate():
    for name, kw in matrix():
        args = [sys.executable, f"{ROOT}/scripts/make_case.py", f"{ROOT}/cases/{name}", "--solver", "steady",
                "--iters", str(kw.pop("iters", ITERS))]
        for k, v in kw.items():
            args += [f"--{k}", str(v)]
        subprocess.run(args, check=True)
    print(f"generated {len(matrix())} cases")


def run_one(name):
    c = f"{ROOT}/cases/{name}"
    if os.path.exists(f"{c}/done") or os.path.exists(f"{c}/log.blockMesh"):
        return name, "skipped"
    r = subprocess.run(["bash", f"{c}/Allrun"], capture_output=True, text=True)
    if r.returncode == 0:
        open(f"{c}/done", "w").write("ok\n")
    return name, f"rc={r.returncode} {r.stderr[-300:]}"


if __name__ == "__main__":
    if sys.argv[1] == "generate":
        generate()
    else:
        j = 8
        names = sys.argv[2:]
        if names[:1] == ["-j"]:
            j, names = int(names[1]), names[2:]
        names = names or [n for n, _ in matrix()]
        # longest first
        names.sort(key=lambda n: ("meshFine" not in n, n))
        with ThreadPoolExecutor(j) as ex:
            for name, st in ex.map(run_one, names):
                print(name, st, flush=True)
