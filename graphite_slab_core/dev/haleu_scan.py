"""HALEU (19.75 wt% U-235) study: unit-cell k-inf for cases A/B, bare-cylinder R scan, critical radius.
  python haleu_scan.py kinf            -> results/haleu_kinf.csv
  python haleu_scan.py rscan A B       -> results/haleu_R_scan.csv  (appends/replaces rows for the given cases)
  python haleu_scan.py confirm A=Rc .. -> adds confirmation rows (kind='confirm') at the given radii
  python haleu_scan.py uf4 A           -> results/haleu_uf4_kinf.csv (unit-cell k-inf vs UF4 mol%, UF4 taken from LiF; not used)
Env: NP/NB/NI override particles/batches/inactive; UF4X = UF4 mol% (default 4.0)."""
import os, sys, time, math
import numpy as np, pandas as pd, openmc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import build_model, analytic_metrics, make_materials, DEFAULT_FUEL, HALEU_U235_WT_PCT, UF4_MOLPCT_DEFAULT

WORK = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
RES = os.path.join(WORK, "results"); os.makedirs(RES, exist_ok=True)
GEO = dict(slot_depth=1.0, flat_width=0.5, web_thickness=1.5, wall_thickness=0.5, round_location="both_sides", stacking="plates")
CASES = {"A": dict(GEO), "B": dict(GEO, coolant_depth=0.6)}
ENR, T, THREADS = HALEU_U235_WT_PCT, 922.0, 8
UF4 = float(os.environ.get("UF4X", UF4_MOLPCT_DEFAULT))   # UF4 mol% (default 4.0; 0.83 = original MSRE salt)
RUN = dict(particles=int(os.environ.get("NP", 10000)), batches=int(os.environ.get("NB", 100)),
           inactive=int(os.environ.get("NI", 40)), temperature=T)

def run(m, cwd):
    t0 = time.time()
    sp_path = m.run(cwd=cwd, output=False, threads=THREADS)
    with openmc.StatePoint(sp_path) as sp:
        gt = sp.global_tallies
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in gt["name"]]
        return sp.keff.nominal_value, sp.keff.std_dev, float(gt["mean"][names.index("leakage")]), time.time() - t0

def masses(case, R, fuel=None):
    """core volume, fuel-salt volume/mass and U-235 mass in the H=2R cylinder (fraction x volume)."""
    am = analytic_metrics(core_radius=R, **CASES[case])
    f = make_materials(fuel=fuel, temperature=T, enrichment=ENR, uf4_mol_pct=UF4)["fuel"]
    Vf_cc = am["core_fuel_volume_l"] * 1000
    return dict(core_volume_m3=am["core_volume_l"] / 1000, fuel_salt_volume_m3=Vf_cc / 1e6,
                fuel_salt_mass_kg=Vf_cc * f.density / 1000, u235_mass_kg=Vf_cc * f.get_mass_density("U235") / 1000,
                u_mass_kg=Vf_cc * sum(f.get_mass_density(n) for n in f.get_nuclides() if n.startswith("U")) / 1000)

def upsert(path, df, keys):
    if os.path.exists(path):
        old = pd.read_csv(path)
        old = old[~old.set_index(keys).index.isin(df.set_index(keys).index)]
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(path, index=False); return df

cmd = sys.argv[1]
if cmd == "kinf":
    rows = []
    for c, g in CASES.items():
        k, s, _, dt = run(build_model(**g, mode="unit_cell", enrichment=ENR, uf4_mol_pct=UF4, **RUN),
                          os.path.join(WORK, "haleu", f"{c}_kinf_uf4_{UF4:g}"))
        am = analytic_metrics(**g)
        rows.append(dict(case=c, coolant_depth=am["coolant_depth"], enrichment=ENR, uf4_mol_pct=UF4, fuel_density=make_materials(enrichment=ENR, uf4_mol_pct=UF4)["fuel"].density, kinf=k, kinf_std=s,
                         fuel_vf=am["fuel_vf"], coolant_vf=am["coolant_vf"], graphite_vf=am["graphite_vf"],
                         graphite_to_fuel=am["graphite_to_fuel"], runtime_s=dt, **{k_: RUN[k_] for k_ in ("particles", "batches", "inactive")}))
        print(rows[-1], flush=True)
    upsert(os.path.join(RES, "haleu_kinf.csv"), pd.DataFrame(rows), ["case", "uf4_mol_pct"])
elif cmd in ("rscan", "confirm"):
    if cmd == "rscan":
        radii = [float(r) for r in os.environ.get("RADII", "70,100,130,160,200").split(",")]
        jobs = [(c, R) for c in sys.argv[2:] for R in radii]
    else:
        jobs = [(a.split("=")[0], float(a.split("=")[1])) for a in sys.argv[2:]]
    rows = []
    for c, R in jobs:
        m = build_model(**CASES[c], mode="cylinder", core_radius=R, reflector_thickness=0.0, enrichment=ENR, uf4_mol_pct=UF4, **RUN)
        k, s, L, dt = run(m, os.path.join(WORK, "haleu", f"{c}_R{R:.1f}"))
        rows.append(dict(case=c, kind="scan" if cmd == "rscan" else "confirm", uf4_mol_pct=UF4, enrichment=ENR, core_radius=R, core_height=2 * R, keff=k, keff_std=s,
                         leakage_fraction=L, lattice_shape=str(m.params["lattice_shape"]), runtime_s=dt, **masses(c, R),
                         **{k_: RUN[k_] for k_ in ("particles", "batches", "inactive")}))
        print({k_: (round(v, 5) if isinstance(v, float) else v) for k_, v in rows[-1].items()}, flush=True)
    upsert(os.path.join(RES, "haleu_R_scan.csv"), pd.DataFrame(rows), ["case", "kind", "core_radius"])
elif cmd == "uf4":
    base = DEFAULT_FUEL["composition_mol"]
    rows = []
    for c in sys.argv[2:]:
        for x in [float(v) for v in os.environ.get("UF4", "0.83,1.5,2.5,4.0").split(",")]:
            comp = dict(base, UF4=x, LiF=base["LiF"] - (x - base["UF4"]))     # extra UF4 taken from LiF
            fuel = dict(composition_mol=comp, density_ref_mol=base)
            k, s, _, dt = run(build_model(**CASES[c], mode="unit_cell", enrichment=ENR, fuel=fuel, **RUN),
                              os.path.join(WORK, "haleu", f"{c}_uf4_{x}"))
            f = make_materials(fuel=fuel, temperature=T, enrichment=ENR)["fuel"]
            rows.append(dict(case=c, uf4_mol_pct=x, lif_mol_pct=comp["LiF"], fuel_density=f.density,
                             u235_g_per_cc_fuel=f.get_mass_density("U235"), kinf=k, kinf_std=s, runtime_s=dt))
            print(rows[-1], flush=True)
    upsert(os.path.join(RES, "haleu_uf4_kinf.csv"), pd.DataFrame(rows), ["case", "uf4_mol_pct"])
