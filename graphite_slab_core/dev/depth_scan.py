"""One-at-a-time slot-depth scans (flat 0.5, web 1.5, wall 0.5 cm fixed).
  python depth_scan.py common   -> fuel and coolant depth together (d = 0.6..1.2)  -> results/depth_scan.csv
  python depth_scan.py coolant  -> coolant depth only, fuel depth 1.0              -> results/coolant_depth_scan.csv
Runs unit-cell k-inf and bare-cylinder k-eff for each depth.  Plots are made in notebook section 9."""
import os, sys, time, json, math
import numpy as np, pandas as pd, openmc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import build_model, analytic_metrics, make_materials, MSRE_ISOTOPICS

WORK = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
SCAN = sys.argv[1] if len(sys.argv) > 1 else "common"
assert SCAN in ("common", "coolant")
DEPTHS = [0.6, 0.8, 1.0, 1.2]
FUEL_DEPTH = 1.0
FIXED = dict(flat_width=0.5, web_thickness=1.5, wall_thickness=0.5)
OPTS = dict(round_location="both_sides", stacking="plates")
R, ENR, T = 70.0, 33.477, 922.0
RUN = dict(particles=int(os.environ.get("NP", 10000)), batches=int(os.environ.get("NB", 100)),
           inactive=int(os.environ.get("NI", 40)), temperature=T)
THREADS = 8

def run(m, cwd):
    t0 = time.time()
    sp_path = m.run(cwd=cwd, output=False, threads=THREADS)
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
        gt = sp.global_tallies
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in gt["name"]]
        leak = float(gt["mean"][names.index("leakage")])
    return k.nominal_value, k.std_dev, leak, time.time() - t0

FUEL = dict(fuel=MSRE_ISOTOPICS, uf4_mol_pct=0.83)   # depth scans: original MSRE fuel (0.83 mol% UF4, MSRE uranium)
rho_fuel = make_materials(temperature=T, enrichment=ENR, **FUEL)["fuel"].density   # g/cc
rows = []
for d in DEPTHS:
    geo = dict(slot_depth=d, **FIXED) if SCAN == "common" else dict(slot_depth=FUEL_DEPTH, coolant_depth=d, **FIXED)
    am = analytic_metrics(core_radius=R, **geo, **OPTS)
    base = os.path.join(WORK, "depth_scan" if SCAN == "common" else "coolant_depth_scan", f"d{d:.1f}")
    kinf, kinf_s, _, t1 = run(build_model(**geo, **OPTS, mode="unit_cell", enrichment=ENR, **FUEL, **RUN), base + "_kinf")
    m = build_model(**geo, **OPTS, mode="cylinder", core_radius=R, reflector_thickness=0.0, enrichment=ENR, **FUEL, **RUN)
    keff, keff_s, leak, t2 = run(m, base + "_cyl")
    r = dict(slot_depth=geo["slot_depth"], coolant_depth=am["coolant_depth"], flat_width=FIXED["flat_width"],
             slot_width=am["slot_width"], coolant_slot_width=am["coolant_slot_width"], web_thickness=FIXED["web_thickness"],
             wall_thickness=FIXED["wall_thickness"], pitch_x=am["pitch_x"], pitch_y=am["pitch_y"], pitch_z=am["pitch_z"],
             kinf=kinf, kinf_std=kinf_s, keff=keff, keff_std=keff_s, leakage_fraction=leak,
             fuel_vf=am["fuel_vf"], coolant_vf=am["coolant_vf"], graphite_vf=am["graphite_vf"],
             graphite_to_fuel=am["graphite_to_fuel"],
             core_fuel_volume_l=am["core_fuel_volume_l"], core_coolant_volume_l=am["core_coolant_volume_l"],
             core_fuel_mass_kg=am["core_fuel_volume_l"] * rho_fuel,
             fuel_density_g_cc=rho_fuel, lattice_shape=str(m.params["lattice_shape"]),
             particles=RUN["particles"], batches=RUN["batches"], inactive=RUN["inactive"],
             runtime_kinf_s=t1, runtime_keff_s=t2)
    print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
    rows.append(r)
df = pd.DataFrame(rows)
df.to_csv(os.path.join(WORK, "results", "depth_scan.csv" if SCAN == "common" else "coolant_depth_scan.csv"), index=False)
print(df.round(5).to_string())
