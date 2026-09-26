"""Fuel slots with FIXED total width 1.0 cm (both_sides profile: flat = 1.0 - 2 d_f, so d_f <= 0.5); coolant slots at the
default (depth 1.0, flat 0.5 -> width 2.5 cm); wall 0.5 cm; HALEU 19.75 wt%, 4 mol% UF4.
  python fuel_depth_w1_scan.py depth  -> d_f = 0.20..0.50: unit-cell k-inf + bare R = 85 cm (H = 2R) k-eff  -> results/fuel_depth_w1_scan.csv
  python fuel_depth_w1_scan.py web    -> d_f = 0.5, web = 0.5..3.0 cm (extended with WEBS=4.0,5.0): unit-cell k-inf only                -> results/fuel_w1_web_scan.csv
Env: DEPTHS, WEBS, NP/NB/NI."""
import os, sys, time
import pandas as pd, openmc
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core import build_model, analytic_metrics, make_materials, HALEU_U235_WT_PCT, UF4_MOLPCT_DEFAULT

WORK = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
RES = os.path.join(WORK, "results"); os.makedirs(RES, exist_ok=True)
W_FUEL, R = 1.0, 85.0
ENR, UF4, T, THREADS = HALEU_U235_WT_PCT, UF4_MOLPCT_DEFAULT, 922.0, 8
RUN = dict(particles=int(os.environ.get("NP", 10000)), batches=int(os.environ.get("NB", 100)),
           inactive=int(os.environ.get("NI", 40)), temperature=T, enrichment=ENR, uf4_mol_pct=UF4)

def geo(d_f, web=1.5):
    return dict(slot_depth=d_f, flat_width=round(W_FUEL - 2 * d_f, 9), coolant_depth=1.0, coolant_flat_width=0.5,
                web_thickness=web, wall_thickness=0.5, round_location="both_sides", stacking="plates")

def run(m, cwd):
    t0 = time.time()
    sp_path = m.run(cwd=cwd, output=False, threads=THREADS)
    with openmc.StatePoint(sp_path) as sp:
        gt = sp.global_tallies
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in gt["name"]]
        return sp.keff.nominal_value, sp.keff.std_dev, float(gt["mean"][names.index("leakage")]), time.time() - t0

fuel = make_materials(temperature=T, enrichment=ENR, uf4_mol_pct=UF4)["fuel"]
def common(g, am):
    Vf = am["core_fuel_volume_l"] * 1000
    return dict(fuel_depth=g["slot_depth"], fuel_flat=g["flat_width"], fuel_slot_width=am["slot_width"], web_thickness=g["web_thickness"],
                coolant_slot_width=am["coolant_slot_width"], pitch_x=am["pitch_x"], pitch_y=am["pitch_y"], pitch_z=am["pitch_z"],
                fuel_vf=am["fuel_vf"], coolant_vf=am["coolant_vf"], graphite_vf=am["graphite_vf"], graphite_to_fuel=am["graphite_to_fuel"],
                fuel_salt_volume_m3=Vf / 1e6, fuel_salt_mass_kg=Vf * fuel.density / 1000, u235_mass_kg=Vf * fuel.get_mass_density("U235") / 1000,
                rel_conduction_power_limit=(1.0 / g["slot_depth"]) ** 2)

mode = sys.argv[1] if len(sys.argv) > 1 else "depth"
rows = []
if mode == "depth":
    for d in [float(x) for x in os.environ.get("DEPTHS", "0.20,0.25,0.30,0.35,0.40,0.45,0.50").split(",")]:
        g = geo(d); am = analytic_metrics(core_radius=R, **g)
        base = os.path.join(WORK, "fuel_w1", f"df{d:.2f}")
        ki, ks, _, t1 = run(build_model(**g, mode="unit_cell", **RUN), base + "_kinf")
        kb, kbs, L, t2 = run(build_model(**g, mode="cylinder", core_radius=R, reflector_thickness=0.0, **RUN), base + "_bare")
        rows.append(dict(**common(g, am), kinf=ki, kinf_std=ks, keff_bare=kb, keff_bare_std=kbs, leakage_bare=L,
                         core_radius=R, core_height=2 * R, runtime_s=t1 + t2))
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(RES, "fuel_depth_w1_scan.csv"), index=False)
else:
    for web in [float(x) for x in os.environ.get("WEBS", "0.5,1.0,1.5,2.0,2.5,3.0").split(",")]:
        g = geo(0.5, web); am = analytic_metrics(core_radius=R, **g)
        ki, ks, _, t1 = run(build_model(**g, mode="unit_cell", **RUN), os.path.join(WORK, "fuel_w1", f"web{web:.2f}_kinf"))
        rows.append(dict(**common(g, am), kinf=ki, kinf_std=ks, runtime_s=t1))
        print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rows[-1].items()}, flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(RES, "fuel_w1_web_scan.csv"), index=False)
