"""Graphite reflector thickness scan (radial AND axial) for the default HALEU core, case A geometry.
  python reflector_scan.py scan          -> R = 85 cm, t = 0,10,20,30,45,60 cm      (kind='scan')
  python reflector_scan.py crit 30 R1 R2 -> reflected core, t = 30 cm, radii R1..   (kind='crit')
Results are upserted into results/reflector_scan.csv.  Env NP/NB/NI override statistics; THICK = comma list."""
import os, sys, time
import pandas as pd, openmc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import build_model, analytic_metrics, make_materials, HALEU_U235_WT_PCT, UF4_MOLPCT_DEFAULT

WORK = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
RES = os.path.join(WORK, "results"); os.makedirs(RES, exist_ok=True)
GEO = dict(slot_depth=1.0, flat_width=0.5, web_thickness=1.5, wall_thickness=0.5, round_location="both_sides", stacking="plates")
ENR, UF4, T, THREADS = HALEU_U235_WT_PCT, UF4_MOLPCT_DEFAULT, 922.0, 8
RUN = dict(particles=int(os.environ.get("NP", 10000)), batches=int(os.environ.get("NB", 100)),
           inactive=int(os.environ.get("NI", 40)), temperature=T)

def run(R, t, kind):
    m = build_model(**GEO, mode="cylinder", core_radius=R, reflector_thickness=t, reflector_axial=t,
                    enrichment=ENR, uf4_mol_pct=UF4, **RUN)
    t0 = time.time()
    sp_path = m.run(cwd=os.path.join(WORK, "reflector", f"R{R:.1f}_t{t:.1f}"), output=False, threads=THREADS)
    with openmc.StatePoint(sp_path) as sp:
        gt = sp.global_tallies
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in gt["name"]]
        k, s, L = sp.keff.nominal_value, sp.keff.std_dev, float(gt["mean"][names.index("leakage")])
    am = analytic_metrics(core_radius=R, **GEO)
    f = make_materials(temperature=T, enrichment=ENR, uf4_mol_pct=UF4)["fuel"]
    Vf = am["core_fuel_volume_l"] * 1000
    row = dict(kind=kind, core_radius=R, core_height=2 * R, reflector_radial=t, reflector_axial=t, keff=k, keff_std=s,
               leakage_fraction=L, outer_diameter=2 * (R + t), outer_height=2 * (R + t),
               reflector_graphite_volume_m3=(3.141592653589793 * (R + t) ** 2 * 2 * (R + t) - 3.141592653589793 * R ** 2 * 2 * R) / 1e6,
               fuel_salt_volume_m3=Vf / 1e6, fuel_salt_mass_kg=Vf * f.density / 1000, u235_mass_kg=Vf * f.get_mass_density("U235") / 1000,
               runtime_s=time.time() - t0, enrichment=ENR, uf4_mol_pct=UF4, **{k_: RUN[k_] for k_ in ("particles", "batches", "inactive")})
    print({k_: (round(v, 5) if isinstance(v, float) else v) for k_, v in row.items()}, flush=True)
    return row

cmd = sys.argv[1]
if cmd == "scan":
    rows = [run(85.0, t, "scan") for t in [float(x) for x in os.environ.get("THICK", "0,10,20,30,45,60").split(",")]]
else:
    t = float(sys.argv[2]); rows = [run(float(R), t, "crit") for R in sys.argv[3:]]
df = pd.DataFrame(rows); path = os.path.join(RES, "reflector_scan.csv")
if os.path.exists(path):
    old = pd.read_csv(path); keys = ["kind", "core_radius", "reflector_radial"]
    old = old[~old.set_index(keys).index.isin(df.set_index(keys).index)]
    df = pd.concat([old, df], ignore_index=True)
df.to_csv(path, index=False)
