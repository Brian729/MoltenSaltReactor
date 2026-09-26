"""Minimal ascii OpenFOAM field/mesh readers for the fuel-slot cases (numpy only)."""
import os, re
import numpy as np
from set_fields import read_mesh_centres  # noqa: F401  (re-exported)


def _list_after(txt, key_pos):
    m = re.compile(r"(uniform\s+(\([^()]*\)|[-+.eE\d]+))|(nonuniform\s+List<(\w+)>\s*(\d+)\s*\()").search(txt, key_pos)
    return m


def read_internal(path, ncell=None):
    txt = open(path).read()
    i = txt.index("internalField")
    m = _list_after(txt, i)
    if m.group(1):
        v = m.group(2)
        if v.startswith("("):
            return np.tile(np.array(v[1:-1].split(), float), (ncell, 1))
        return np.full(ncell, float(v))
    typ, n = m.group(4), int(m.group(5))
    body = txt[m.end():]
    if typ == "vector":
        return np.array([[float(x) for x in s.split()] for s in re.findall(r"\(([^()]*)\)", body)[:n]])
    return np.array(body[:body.index(")")].split(), float)[:n]


def read_patch(path, patch):
    txt = open(path).read()
    i = txt.index("boundaryField")
    j = re.compile(r"\n\s*" + re.escape(patch) + r"\s*\n\s*\{").search(txt, i)
    k = re.compile(r"\n\s*value\s").search(txt, j.end()).end()
    m = _list_after(txt, k)
    if m.group(1):
        return float(m.group(2))
    n = int(m.group(5))
    body = txt[m.end():]
    return np.array(body[:body.index(")")].split(), float)[:n]


def latest_time(case):
    ts = []
    for d in os.listdir(case):
        try:
            if float(d) > 0 and os.path.isdir(f"{case}/{d}/fuel"):
                ts.append((float(d), d))
        except ValueError:
            pass
    return max(ts)[1] if ts else None


def all_times(case):
    ts = []
    for d in os.listdir(case):
        try:
            if float(d) > 0 and os.path.isfile(f"{case}/{d}/fuel/T"):
                ts.append((float(d), d))
        except ValueError:
            pass
    return [t[1] for t in sorted(ts)]
