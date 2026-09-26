import nbformat as nbf
import re as _re
_src = open("core.py").read()
_parts = _re.split(r"# =+\n# (\d)\. (.*)\n# =+\n", _src)
chunks = {}
for i in range(1, len(_parts), 3):
    chunks[_parts[i]] = f"# {_parts[i]}. {_parts[i+1]}\n" + _parts[i+2].strip("\n") + "\n"
_pl = open("plots.py").read()
PLOT_HELPERS = _pl[_pl.index("MAT_COLORS ="):_pl.index("\ndef sketch_profiles")].strip() + "\n"
SKETCH = _pl[_pl.index("def sketch_profiles"):].strip() + "\n"
_lh = open("lhs.py").read(); lhs = _lh[_lh.index("PARAM_NAMES ="):]

cells = []
md = lambda s: cells.append(nbf.v4.new_markdown_cell(s.strip("\n")))
code = lambda s: cells.append(nbf.v4.new_code_cell(s.strip("\n")))

md(r'''
# Slotted graphite-plate MSR core — OpenMC CSG model + Latin-hypercube study

**Concept (Brian):** graphite plates with machined slots; plates are stacked with every other plate rotated 90°, so
rows of **vertical fuel slots** alternate with rows of **horizontal coolant slots**. All slots have the **same
machined profile and size**: depth `d`, both side walls **quarter-rounds of radius `d`**, flat bottom `flat_width`
between them → slot width **`w = 2d + flat_width`**. Neighbouring slots in a row are separated by graphite **webs**
(`web_thickness`); each slot row is closed by the solid backing of the next plate — **one shared wall**
(`wall_thickness`) between every fuel row and coolant row.

This notebook:
1. defines the geometry with OpenMC CSG (`build_model(slot_depth, flat_width, web_thickness, wall_thickness, **opts)`),
2. plots it (profile sketch, full-core and close-up xy / xz / yz slices, unit cell; coloured by material),
3. checks the CSG volumes against closed-form area formulas,
4. runs the first model (bare finite cylinder, R = 70 cm, H = 2R, MSRE fuel salt) and a quick unit-cell k-inf,
5. sets up a Latin-hypercube sweep over `slot_depth`, `flat_width`, `web_thickness`, `wall_thickness`
   (`scipy.stats.qmc.LatinHypercube`; enrichment and core radius optional), transport switched by `RUN_TRANSPORT`.

Everything is in this notebook (no helper modules needed). Tested with OpenMC 0.16.0, SciPy 1.18, pandas 3.0.
''')

md(r'''
## Geometry interpretation & assumptions  (edit here first if anything is wrong)

**Axes:** `x` = plate-stacking direction, `y` = horizontal in the plate plane, `z` = vertical.

**Plate stack (confirmed by Brian):** slots are machined into one face of a graphite plate (depth `d`); the rest of the
plate is a solid backing of thickness `wall_thickness`. Plates are stacked along x, alternately rotated 90°, so the
rows alternate fuel / coolant and each row's open mouth is closed by the backing of the next plate. The backing is
therefore **one shared wall** between each fuel row and the adjacent coolant row (`stacking="plates"`, default):

```
 x ──►  |<-- fuel plate -->|<-- coolant plate -->|<-- fuel plate ...
        | F-row  |  wall   | C-row   |  wall     | F-row
        |<- d -->|<t_wall->|<- d --->|<t_wall-->|
   F-row = row of vertical fuel slots      (run along z, width along y, pitch w + t_web in y)
   C-row = row of horizontal coolant slots (run along y, width along z, pitch w + t_web in z)
   unit cell: P_x = 2 (d + t_wall),  P_y = P_z = w + t_web
```

**Shared machined-slot profile** (one function `milled_slot_region()` builds both slot types, just oriented differently):

```
  u (across width) ──►
  0     d            w-d     w
  +=====+=============+=======+   <- mouth: open face (closed by the next plate's backing wall)
   '.   |             |    .'
     '. |    salt     | .'         both side walls: quarter-rounds of radius R = d,
       '+-------------+'           centres on the mouth plane at u = d and u = w-d
        |<- flat_w -->|           <- flat bottom at depth v = d
  v (depth into graphite, +x)          w = 2d + flat_width   (flat_width = 0 -> half-round)
```
* Default `round_location="both_sides"` (full-radius corners on both sides). CSG:
  `slot = [cyl(R=d, centre u=d) ∩ (u ≤ d)] ∪ [rectangle flat_width × d] ∪ [cyl(R=d, centre u=w−d) ∩ (u ≥ w−d)]`, all ∩ (depth ≥ 0);
  cylinder axes = slot run direction. Area `A = flat_width·d + π d²/2`. Graphite = complement of all slots.
* Fuel slot: profile in the **x–y** plane, extruded along **z**. Coolant slot: the same profile in the **x–z** plane, extruded along **y**.
* **All slots identical** (fuel = coolant width, depth, profile) and **one web thickness everywhere**: `web_thickness` is the land
  between neighbouring slots of a row (the only web in the plate stack; rows are separated by the shared wall).
* Optional `coolant_depth=` (default `None` = same as `slot_depth`) gives the coolant slots their own depth (width `2 d_c + flat`); used only in the §9 coolant-only depth scan.
* Every sample is feasible by construction (`d > 0`, `flat_width ≥ 0`, `web > 0`, `wall > 0`).
* Legacy options (not default): `round_location="side"` (one side rounded, w = d + flat), `"bottom"` (U-groove, explicit
  `slot_width ≤ 2d`), `"none"` (square, w = flat); `stacking="interleaved"` (old multi-row slabs with webs between rows).
* **Not modelled:** the run-out radius at the *ends* of a slot; slots are infinite in the unit cell and cut off at the cylinder surface.

**Core and materials:**
* **Default model = bare finite cylindrical core** (`mode="cylinder"`, `reflector_thickness = 0`): a `RectLattice` of the unit-cell
  universe (odd count per axis, a unit cell centred on the axis) truncated by a `ZCylinder` of radius `R = core_radius`
  (default **70 cm**, ~MSRE core radius) and z-planes at ±R (**H = 2R = D**), vacuum boundary. No plenums, vessel, downcomer
  or headers — **the fuel salt outside the core is not modelled**. `mode="unit_cell"` = periodic unit cell → k-infinity (quick mode).
* Materials at 922 K (MSRE operating temperature): **fuel = MSRE 235U-operation fuel salt** 7LiF-BeF₂-ZrF₄-UF₄ 65.0-29.17-5.0-0.83 mol%,
  33.477 wt% U-235 (variable `ENRICHMENT`), 99.995 % Li-7, ρ = 2.575 − 5.13·10⁻⁴·T[°C] g/cc (ORNL sources cited in §2);
  coolant = MSRE coolant salt 7LiF-BeF₂ 66-34 mol%; graphite 1.87 g/cc + `c_Graphite` S(α,β).
  No tube walls/liners — the salt wets the graphite directly.

### Note on wall thickness (rule of thumb — not a qualified design)
* Start at **5 mm** (default `wall_thickness = 0.5` cm); **~3 mm** is a practical floor for machining and handling nuclear graphite.
* **Pressure is not limiting.** Treating the wall over one slot as a plate strip of span `w` and thickness `t` loaded by the
  fuel/coolant pressure difference `p` (< 5 psi ≈ 0.034 MPa): σ ≈ p·w²/(2t²) ≈ **0.4 MPa** for the defaults (w = 2.5 cm, t = 5 mm),
  and only ~6 MPa at the most extreme LHS corner (w = 5.5 cm, t = 3 mm) — compared with ~30–50 MPa flexural strength of nuclear graphite.
* The real limits are **machinability**, **irradiation-induced dimensional change and the resulting internal stress over ~5 years**,
  **salt permeation / fuel–coolant cross-leak** through a thin wall, and component tolerances. The wall's **thermal
  resistance is minor** compared with the salt-side film resistances.
''')

code(r'''
# ---- user switches -----------------------------------------------------------------------------------
import os, math, json, shutil, warnings, re, time, textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import openmc
from scipy.stats import qmc

WORK_DIR = os.path.abspath(os.environ.get("MSR_WORK_DIR", "."))      # PNGs / CSV / run dirs go here
RUN_KINF       = True      # quick unit-cell k-infinity run of the default geometry (~1 min)
RUN_BASELINE   = True      # finite-cylinder k-eff run of the default geometry (the "first model")
RUN_TRANSPORT  = False     # transport for every LHS sample (set True for the sweep; skipped without nuclear data)
PARTICLES, BATCHES, INACTIVE = 10000, 100, 40   # finite core: 600k active histories -> k std ~100 pcm
# all transport steps are skipped automatically if no cross_sections.xml is found

# ---- core --------------------------------------------------------------------------------------------------
MODE = "cylinder"          # "cylinder" (finite core, H = 2R, vacuum boundary) or "unit_cell" (periodic, k-inf)
CORE_RADIUS = 70.0         # cm; H = 2R = 140 cm.  (MSRE graphite core radius was 70.2 cm - IRPhEP benchmark)
REFLECTOR_THICKNESS = 0.0  # cm of graphite around the cylinder (side, top, bottom); 0 = bare core
THREADS = None             # None = all OpenMP threads

# ---- default geometry (cm) --------------------------------------------------------------------------------
DEFAULTS = dict(slot_depth=1.0, flat_width=0.5, web_thickness=1.5, wall_thickness=0.5)   # -> slot width w = 2.5 cm
GEOM_OPTS = dict(round_location="both_sides", stacking="plates")                          # Brian's design (defaults)

# ---- fuel ------------------------------------------------------------------------------------------------
ENRICHMENT = 33.477        # U-235 WEIGHT % of total uranium. Default = MSRE 235U operation (ORNL-4658 Table 2.8)
TEMPERATURE_K = 922.0      # material temperature (MSRE operating ~650 C); salt densities follow ORNL correlations

import sys
if shutil.which("openmc") is None:                 # make the conda env's `openmc` executable visible to model.run()
    os.environ["PATH"] = os.path.join(sys.prefix, "bin") + os.pathsep + os.environ.get("PATH", "")
USE_NUCLEAR_DATA = os.environ.get("MSR_USE_XS", "1") == "1"   # set False to force the no-data path

# nuclear data: $OPENMC_CROSS_SECTIONS or first existing candidate below
for _xs in ([] if not USE_NUCLEAR_DATA else [os.environ.get("OPENMC_CROSS_SECTIONS"),
            os.path.expanduser("~/nucdata/endfb-viii.0-hdf5/cross_sections.xml"),
            "/workspace/nucdata/endfb-viii.0-hdf5/cross_sections.xml"]):
    if _xs and os.path.exists(_xs):
        openmc.config["cross_sections"] = _xs
        break
if not USE_NUCLEAR_DATA:
    openmc.config.pop("cross_sections", None)
HAVE_XS = openmc.config.get("cross_sections") is not None and os.path.exists(str(openmc.config.get("cross_sections")))
print("OpenMC", openmc.__version__, "| cross sections:", openmc.config.get("cross_sections") if HAVE_XS else "NOT FOUND (transport skipped, python plot fallback)")
print("WORK_DIR =", WORK_DIR)
''')

md("## 1. Shared machined-slot profile (fuel **and** coolant): double-rounded slot, w = 2d + flat")
code(chunks["1"])
code(SKETCH + '''
fig = sketch_profiles(d=DEFAULTS["slot_depth"], flat=DEFAULTS["flat_width"], filename=os.path.join(WORK_DIR, "profile_sketch.png"))
plt.show()''')

md(r"""
## 2. Materials — default fuel = MSRE fuel salt (U-235 operation), per ORNL sources

| Quantity | Value used (default) | Source |
|---|---|---|
| Fuel salt composition | 7LiF-BeF₂-ZrF₄-UF₄ **65.0-29.17-5.0-0.83 mol%** | R. E. Thoma, *Chemical Aspects of MSRE Operations*, ORNL-4658 (1971), p. 10–11 (fuel for 235U operation; also Table 1.1: 65-29.2-5-0.83) |
| Uranium isotopics (start of power operation, run 4-1) | U-234 0.342, **U-235 33.477**, U-236 0.141, U-238 66.041 wt% | ORNL-4658 Table 2.8 |
| Li-7 in fuel carrier salt | **99.995 at%** (batch assays 99.994–99.996) | ORNL-4658 Table 2.11; same value used in the IRPhEP MSRE benchmark (Shen/Fratoni et al., PHYSOR 2020) |
| Fuel density | **ρ = 2.575 − 5.13×10⁻⁴·T(°C) g/cm³** (±1 %) → 2.242 g/cm³ at 649 °C (139.9 lb/ft³ at 650 °C) | ORNL-4658 Table 8.2 (Cantor molar-volume method; Table 8.3) |
| Coolant salt | 7LiF-BeF₂ 66-34 mol%, 99.992 % Li-7, ρ = 2.214 − 4.2×10⁻⁴·T(°C) → 1.941 g/cm³ at 649 °C | ORNL-4658 Tables 2.1, 8.1; ORNL-4616 |
| Graphite | 1.87 g/cm³ (MSRE grade CGB), pure C + `c_Graphite` | IRPhEP MSRE benchmark evaluation (1.87 ± 0.02 g/cm³) |

**Enrichment variable:** `ENRICHMENT` (top cell) = U-235 **weight %** of total uranium; also `build_model(..., enrichment=...)`,
`make_materials(enrichment=...)` and an optional 5th LHS dimension (`LHS_BOUNDS["enrichment"]`, off by default).
U-234/U-236 scale proportionally with enrichment (`minor_u="scale"`; exact MSRE isotopics at 33.477 wt%);
use `fuel={"minor_u": "none"}` for a pure U-235/U-238 vector.
The UF₄ mole fraction stays at 0.83 mol% when enrichment changes (so changing enrichment changes the U-235 loading).

**Alternatives found in ORNL sources (not used by default; switch via the `fuel=` dict):**
* Nominal design composition 65-29.1-5-0.9 mol% (ORNL-TM-728, Robertson 1965, Table 2.1; ORNL-4616; Haubenreich & Engel, *Nucl. Appl. Tech.* 8 (1970)); U "about 32 %" (ORNL-4616) / "33 %" (Haubenreich & Engel) enriched.
* Zero-power first-criticality salt 64.88-29.27-5.06-0.79 mol%, 1.408 wt% U-235 in salt, ρ = 2.3275 ± 0.016 g/cm³ at 638 °C (IRPhEP MSRE benchmark).
* Density: in-reactor inventory estimate 139.01 lb/ft³ (2.227 g/cm³) at start of power operation; electrical-probe correlation ρ = 2.848 − 7.69×10⁻⁴·T(°C) (2.35 g/cm³ at 650 °C); design value 141 lb/ft³ (2.26 g/cm³) (ORNL-4658 §8.2, Table 8.3).
* Li-7: "at least 99.99 %" (ORNL-4616), 99.994 % (ORNL-TM-728 Table 2.1).
""")

code(chunks["2"] + '''
_m = make_materials(temperature=TEMPERATURE_K, enrichment=ENRICHMENT)
for k, m in _m.items():
    print(f"{k:9s} {m.name:55s} {m.density:.4f} g/cc @ {m.temperature:.0f} K  nuclides: {', '.join(n.name for n in m.nuclides)}")
print("uranium vector (wt% of U):", {k: round(100 * v, 3) for k, v in uranium_wt_fractions(
      ENRICHMENT, DEFAULT_FUEL["u234_wt_pct"], DEFAULT_FUEL["u236_wt_pct"], DEFAULT_FUEL["minor_u"]).items()})
print(f"U-235 mass fraction in fuel salt: {100 * _m['fuel'].get_mass_density('U235') / _m['fuel'].density:.3f} wt%")''')

md('''## 3. Parameter validation, layer stack and analytic volume fractions
`resolve_params()` validates the inputs (positive d/web/wall, flat_width ≥ 0), derives `slot_width = 2d + flat_width` and returns the x-layer stack;
`analytic_metrics()` gives exact volume fractions / moderator-to-fuel ratio for the periodic cell (no transport needed).''')
code(chunks["3"] + '''
pd.Series(analytic_metrics()).to_frame("default geometry")''')

md('''## 4. OpenMC geometry: `build_model(...) -> openmc.Model`''')
code(chunks["4"])

code(r'''
RUN_OPTS = dict(particles=PARTICLES, batches=BATCHES, inactive=INACTIVE, temperature=TEMPERATURE_K)

GEOM_OPTS_CORE = dict(GEOM_OPTS, mode=MODE, reflector_thickness=REFLECTOR_THICKNESS)
model = build_model(**DEFAULTS, **GEOM_OPTS_CORE, core_radius=CORE_RADIUS, enrichment=ENRICHMENT, **RUN_OPTS)   # finite core
model_uc = build_model(**DEFAULTS, **GEOM_OPTS, mode="unit_cell", enrichment=ENRICHMENT, **RUN_OPTS)            # k-inf cell
p = model_uc.params
print("layers along x:", p["layers"], f"| slot width w = 2d + flat = {p['slot_width']:.3f} cm")
print(f"unit cell pitch  Px={p['pitch_x']:.3f}  Py={p['pitch_y']:.3f}  Pz={p['pitch_z']:.3f} cm")
if MODE == "cylinder":
    print(f"core: R = {CORE_RADIUS} cm, H = {2*CORE_RADIUS} cm, reflector = {REFLECTOR_THICKNESS} cm, "
          f"lattice {model.params['lattice_shape']} (x, y, z) unit cells")
    print({k: round(v, 1) for k, v in analytic_metrics(**DEFAULTS, **GEOM_OPTS, core_radius=CORE_RADIUS).items() if k.startswith("core_")})
for c in model_uc.geometry.get_all_cells().values():
    print(f"  cell {c.id:3d} {c.name:28s} fill={c.fill.name if hasattr(c.fill,'name') else c.fill}")
os.makedirs(os.path.join(WORK_DIR, "baseline"), exist_ok=True)
model.export_to_model_xml(os.path.join(WORK_DIR, "baseline", "model.xml"))
os.makedirs(os.path.join(WORK_DIR, "baseline_kinf"), exist_ok=True)
model_uc.export_to_model_xml(os.path.join(WORK_DIR, "baseline_kinf", "model.xml"))
''')

md('''## 5. Geometry plots (coloured by material)
OpenMC's plotter needs a `cross_sections.xml`; without one the helper falls back to a pure-Python renderer
(`geometry.find()` per pixel), so plots always work. Full-core xy/xz slices of the cylinder, then close-ups of the slot pattern at the core centre.''')
code(PLOT_HELPERS)
code(r'''
R_ = CORE_RADIUS + REFLECTOR_THICKNESS
Px, Py, Pz, d = p["pitch_x"], p["pitch_y"], p["pitch_z"], p["slot_depth"]
# full core
plot_slice(model, "xy", origin=(0, 0, 0), width=(2.1 * R_, 2.1 * R_), pixels=(1600, 1600),
           title=f"full core, xy (horizontal) slice at z=0 (R = {CORE_RADIUS} cm)", filename=os.path.join(WORK_DIR, "geom_core_xy.png"))
plot_slice(model, "xz", origin=(0, 0, 0), width=(2.1 * R_, 2.1 * R_), pixels=(1600, 1600),
           title=f"full core, xz (vertical) slice at y=0 (H = {2*CORE_RADIUS} cm)", filename=os.path.join(WORK_DIR, "geom_core_xz.png"))
# close-ups of the slot pattern at the core centre (a unit cell is centred on the axis)
W = (4 * Px + 1, 4 * Py + 1, 4 * Pz + 1)
x_fuel = layer_start(p, "F") + 0.35 * d      # inside the fuel-slot row of the central unit cell
x_cool = layer_start(p, "C") + 0.35 * d      # inside the coolant-slot row
plot_slice(model, "xy", origin=(0, 0, 0), width=(W[0], W[1]), pixels=(1200, int(1200 * W[1] / W[0])),
           title="close-up xy at z=0 (through coolant-slot centres): fuel-slot profiles, coolant slots lengthwise",
           filename=os.path.join(WORK_DIR, "geom_xy.png"))
plot_slice(model, "xz", origin=(0, 0, 0), width=(W[0], W[2]), pixels=(1200, int(1200 * W[2] / W[0])),
           title="close-up xz at y=0 (through fuel-slot centres): coolant-slot profiles, fuel slots lengthwise",
           filename=os.path.join(WORK_DIR, "geom_xz.png"))
plot_slice(model, "yz", origin=(x_fuel, 0, 0), width=(W[1], W[2]), pixels=(800, 800),
           title="close-up yz inside a fuel-slot row (vertical fuel slots)", filename=os.path.join(WORK_DIR, "geom_yz_fuel_row.png"))
plot_slice(model, "yz", origin=(x_cool, 0, 0), width=(W[1], W[2]), pixels=(800, 800),
           title="close-up yz inside a coolant-slot row (horizontal coolant slots)", filename=os.path.join(WORK_DIR, "geom_yz_coolant_row.png"))
plt.show()
''')
code(r'''
# zoom on ONE periodic unit cell (the k-inf model): profiles of both slot types
plot_slice(model_uc, "xy", origin=(0, 0, 0), width=(p["pitch_x"], p["pitch_y"]), pixels=(1000, int(1000 * p["pitch_y"] / p["pitch_x"])),
           title="unit cell, xy (fuel-slot profiles; coolant slots cut lengthwise)", filename=os.path.join(WORK_DIR, "geom_unitcell_xy.png"))
plot_slice(model_uc, "xz", origin=(0, 0, 0), width=(p["pitch_x"], p["pitch_z"]), pixels=(1000, int(1000 * p["pitch_z"] / p["pitch_x"])),
           title="unit cell, xz (coolant-slot profiles; fuel slots cut lengthwise)", filename=os.path.join(WORK_DIR, "geom_unitcell_xz.png"))
plt.show()
''')

md('''## 6. Check: CSG volumes (Monte-Carlo point sampling of the OpenMC geometry) vs. closed-form formulas
Pure Python, no nuclear data needed. Agreement within ~2σ confirms the CSG matches the intended profile.''')
code(r'''
rows = []
cases = [("both_sides", None, 0.5), ("both_sides", None, 0.0), ("both_sides", None, 1.5),
         ("side", 2.5, None), ("bottom", 2.0, None), ("none", 2.5, None)]
for loc, w, flat in cases:
    kw = dict(DEFAULTS, round_location=loc)
    if w is not None: kw["slot_width"] = w
    if flat is not None: kw["flat_width"] = flat
    w = resolve_params(**kw)["slot_width"]
    m = build_model(**kw, mode="unit_cell")
    mc = csg_volume_check(m, n=20000, seed=1)
    an = analytic_metrics(**kw)
    for k in ("fuel_vf", "coolant_vf", "graphite_vf"):
        rows.append(dict(round_location=loc, slot_width=w, quantity=k, analytic=an[k], csg_mc=mc[k][0],
                         sigma=mc[k][1], n_sigma=(mc[k][0] - an[k]) / mc[k][1]))
chk = pd.DataFrame(rows)
display(chk.round(4))
assert (chk["n_sigma"].abs() < 4).all(), "CSG and analytic volume fractions disagree!"
''')

md("""## 7. First model: finite-cylinder k-eff at the default parameters (+ quick unit-cell k-inf)""")
code(r"""
def run_and_report(m, subdir, label):
    t0 = time.time()
    sp_path = m.run(cwd=os.path.join(WORK_DIR, subdir), output=False, threads=THREADS)
    dt = time.time() - t0
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
        n_act = sp.n_batches - sp.n_inactive
        gt = sp.global_tallies
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in gt["name"]]
        leak = float(gt["mean"][names.index("leakage")]) if "leakage" in names else None
    res = dict(case=label, keff=k.nominal_value, keff_std=k.std_dev, keff_std_pcm=k.std_dev * 1e5,
               leakage_fraction=leak, particles=m.settings.particles, batches=m.settings.batches,
               inactive=m.settings.inactive, runtime_s=dt, library=str(openmc.config.get("cross_sections")),
               enrichment_u235_wt_pct=ENRICHMENT, temperature_K=TEMPERATURE_K, **DEFAULTS,
               **{k_: v for k_, v in m.params.items() if k_ in ("mode", "core_radius", "core_height", "reflector_thickness", "lattice_shape")})
    print(f"{label}: k = {res['keff']:.5f} +/- {res['keff_std']:.5f} ({res['keff_std_pcm']:.0f} pcm)"
          + (f", leakage = {leak:.4f}" if leak is not None else "")
          + f"  [{res['particles']} particles x {n_act} active / {res['batches']} batches, {dt:.0f} s]")
    return res

baseline = []
if HAVE_XS and RUN_KINF:
    baseline.append(run_and_report(model_uc, "baseline_kinf", "unit cell (periodic) k-inf"))
if HAVE_XS and RUN_BASELINE:
    baseline.append(run_and_report(model, "baseline", f"finite cylinder R={CORE_RADIUS} cm, H={2*CORE_RADIUS} cm k-eff"))
if baseline:
    bdf = pd.DataFrame(baseline)
    bdf.to_csv(os.path.join(WORK_DIR, "baseline_results.csv"), index=False)
    with open(os.path.join(WORK_DIR, "baseline_results.json"), "w") as fh:
        json.dump(baseline, fh, indent=1, default=str)
    display(bdf[["case", "keff", "keff_std_pcm", "leakage_fraction", "particles", "batches", "inactive", "runtime_s"]].round(5))
    if len(baseline) == 2:
        kinf, keff, L = baseline[0]["keff"], baseline[1]["keff"], baseline[1]["leakage_fraction"]
        print(f"k_eff / k_inf = {keff / kinf:.3f}   (1 - leakage = {1 - L:.3f}); "
              f"the bare R = {CORE_RADIUS} cm core is {'sub' if keff < 1 else 'super'}critical by {abs(keff - 1) * 1e5:.0f} pcm")
else:
    print("baseline transport skipped (HAVE_XS =", HAVE_XS, ")")
""")

md('''## 8. Latin-hypercube sampling of the four geometry parameters
* Sampled: `slot_depth` d, `flat_width`, `web_thickness`, `wall_thickness`; the slot width is **derived**, `w = 2d + flat_width`.
* `LHS_BOUNDS` — (low, high) in cm for each free parameter; `FIXED` — hold any parameter constant (it is removed from the hypercube).
* With the double-rounded profile every sample is feasible by construction (the feasibility filter is kept for legacy profiles).
* Each feasible sample gets its own directory `lhs_runs/sample_XXX/` with `model.xml` + `params.json` (+ statepoint if run).''')
code(lhs)
code(r'''
LHS_BOUNDS = {
    "slot_depth":     (0.8, 2.0),
    "flat_width":     (0.0, 1.5),     # slot width w = 2*slot_depth + flat_width  (1.6 - 5.5 cm)
    "web_thickness":  (0.8, 2.0),
    "wall_thickness": (0.3, 1.5),     # 3 mm practical floor (see wall-thickness note)
    # "enrichment":   (5.0, 33.477),   # optional: U-235 wt% of U (off -> ENRICHMENT for every sample)
    # "core_radius":  (50.0, 100.0),   # optional: cylinder radius R in cm, H = 2R (off -> CORE_RADIUS)
}
FIXED = {}                  # e.g. {"wall_thickness": 0.5}
N_SAMPLES, LHS_SEED = 12, 2026
LHS_DIR = os.path.join(WORK_DIR, "lhs_runs")

samples = lhs_samples({k: v for k, v in LHS_BOUNDS.items() if k not in FIXED}, N_SAMPLES, seed=LHS_SEED,
                      fixed=FIXED, geometry_opts=GEOM_OPTS, defaults={"enrichment": ENRICHMENT, "core_radius": CORE_RADIUS})
print(f"{samples['feasible'].sum()} / {len(samples)} samples feasible")
display(samples.round(3))

free = [k for k in ALL_PARAMS if k in LHS_BOUNDS and k not in FIXED]
fig, axs = plt.subplots(len(free), len(free), figsize=(2.3 * len(free), 2.3 * len(free)), squeeze=False)
for i, a in enumerate(free):
    for j, b in enumerate(free):
        ax = axs[i, j]
        if i == j:
            ax.hist(samples[a], bins=N_SAMPLES, range=LHS_BOUNDS[a], color="0.6")
        else:
            for f, c in [(True, "tab:green"), (False, "tab:red")]:
                s = samples[samples["feasible"] == f]
                ax.scatter(s[b], s[a], s=14, c=c, label="feasible" if f else "infeasible")
        if i == len(free) - 1: ax.set_xlabel(b, fontsize=8)
        if j == 0: ax.set_ylabel(a, fontsize=8)
axs[0, -1].legend(fontsize=7)
fig.suptitle("LHS design (one sample per row/column stratum in each 1-D projection)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(WORK_DIR, "lhs_design.png"), dpi=130); plt.show()
''')
code(r'''
do_run = RUN_TRANSPORT and HAVE_XS
if RUN_TRANSPORT and not HAVE_XS:
    print("RUN_TRANSPORT requested but no cross_sections.xml -> building models + analytic metrics only")
results = run_lhs(samples, LHS_DIR, geometry_opts=GEOM_OPTS_CORE, run_opts=RUN_OPTS, run_transport=do_run, threads=THREADS)
csv_path = os.path.join(WORK_DIR, "lhs_results.csv")
results.to_csv(csv_path)
print("saved", csv_path)
display(results[list(PARAM_NAMES) + ["slot_width"] + list(OPTIONAL_PARAMS) + ["fuel_vf", "coolant_vf", "graphite_vf", "graphite_to_fuel", "keff", "keff_std", "runtime_s", "status"]].round(4))
''')
code(r'''
have_k = results["keff"].notna().any()
XCOLS = list(PARAM_NAMES) + [c for c in OPTIONAL_PARAMS if results[c].nunique() > 1]
ycols = ["keff"] if have_k else []
ycols += ["fuel_vf", "graphite_to_fuel"]
fig, axs = plt.subplots(len(ycols), len(XCOLS) + 1, figsize=(3.2 * (len(XCOLS) + 1), 2.8 * len(ycols)), squeeze=False)
for r, y in enumerate(ycols):
    for c, x in enumerate(XCOLS + ["graphite_to_fuel"] if y != "graphite_to_fuel" else XCOLS + ["fuel_vf"]):
        ax = axs[r, c]
        if y == "keff":
            ax.errorbar(results[x], results["keff"], yerr=results["keff_std"], fmt="o", ms=4, capsize=2)
        else:
            ax.plot(results[x], results[y], "o", ms=4)
        ax.set_xlabel(x, fontsize=8); ax.set_ylabel(y, fontsize=8); ax.grid(alpha=0.3)
fig.suptitle("LHS results" + ("" if have_k else "  (no transport: analytic metrics only)"), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(WORK_DIR, "lhs_results.png"), dpi=130); plt.show()
''')
code(r'''
# keff vs each parameter (own figure) + summary table (markdown + PNG)
if have_k:
    xs = XCOLS + ["graphite_to_fuel"]
    fig, axs = plt.subplots(1, len(xs), figsize=(3.4 * len(xs), 3.4), sharey=True)
    for ax, x in zip(axs, xs):
        ax.errorbar(results[x], results["keff"], yerr=results["keff_std"], fmt="o", ms=5, capsize=3)
        for i, r in results.iterrows():
            ax.annotate(str(i), (r[x], r["keff"]), fontsize=6, xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel(x + (" [cm]" if x in PARAM_NAMES else " [U-235 wt%]" if x == "enrichment" else " [cm]" if x == "core_radius" else " (V_gr / V_fuel)")); ax.grid(alpha=0.3)
    axs[0].set_ylabel("k-eff" if MODE == "cylinder" else "k-inf")
    fig.suptitle(f"{'finite-cylinder k-eff' if MODE == 'cylinder' else 'k-inf'} vs parameters (LHS, {len(results)} feasible samples, {PARTICLES} particles x {BATCHES-INACTIVE} active batches)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(WORK_DIR, "keff_vs_params.png"), dpi=140); plt.show()

cols = list(PARAM_NAMES) + ["slot_width"] + list(OPTIONAL_PARAMS) + ["fuel_vf", "coolant_vf", "graphite_to_fuel", "keff", "keff_std", "runtime_s", "status"]
summ = results[cols].copy()
summ["keff_std_pcm"] = summ["keff_std"] * 1e5
summ = summ[list(PARAM_NAMES) + ["slot_width"] + list(OPTIONAL_PARAMS) + ["fuel_vf", "coolant_vf", "graphite_to_fuel", "keff", "keff_std_pcm", "runtime_s", "status"]]
fmt = {c: "{:.3f}" for c in PARAM_NAMES + ("slot_width",)} | {"enrichment": "{:.3f}", "core_radius": "{:.1f}", "runtime_s": "{:.0f}", "fuel_vf": "{:.4f}", "coolant_vf": "{:.4f}", "graphite_to_fuel": "{:.2f}",
                                          "keff": "{:.5f}", "keff_std_pcm": "{:.0f}"}
tbl = summ.copy()
for c, f in fmt.items():
    tbl[c] = [f.format(v) if pd.notna(v) else "-" for v in summ[c]]
md_path = os.path.join(WORK_DIR, "lhs_summary.md")
with open(md_path, "w") as fh:
    fh.write(f"# LHS summary ({len(results)} feasible of {len(samples)} samples, seed {LHS_SEED})\n\n")
    fh.write(f"Lengths in cm; slot_width = 2*slot_depth + flat_width. round_location={GEOM_OPTS['round_location']}, stacking={GEOM_OPTS['stacking']}, "
             f"mode={MODE}, R={CORE_RADIUS} cm (H=2R) unless varied, reflector={REFLECTOR_THICKNESS} cm, "
             f"enrichment {ENRICHMENT} wt% U-235 unless varied, T = {TEMPERATURE_K} K, "
             f"{PARTICLES} particles x {BATCHES} batches ({INACTIVE} inactive), library: {openmc.config.get('cross_sections') if HAVE_XS else 'none'}\n\n")
    hdr = ["sample"] + list(tbl.columns)                      # plain markdown table (no 'tabulate' dependency)
    fh.write("| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n")
    for i, r in tbl.iterrows():
        fh.write("| " + " | ".join([str(i)] + [str(v) for v in r.values]) + " |\n")
    fh.write("\n")
print(open(md_path).read())
fig, ax = plt.subplots(figsize=(14, 0.45 * (len(tbl) + 2)))
ax.axis("off")
t = ax.table(cellText=tbl.values, colLabels=list(tbl.columns), rowLabels=[str(i) for i in tbl.index], loc="center", cellLoc="center")
t.auto_set_font_size(False); t.set_fontsize(8); t.scale(1, 1.3)
ax.set_title("LHS results summary (lengths in cm)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(WORK_DIR, "lhs_summary_table.png"), dpi=150); plt.show()
''')

md(r'''## 9. Depth scans — does a shallower slot improve k-eff?
Two one-at-a-time scans, holding `flat_width = 0.5`, `web_thickness = 1.5`, `wall_thickness = 0.5` cm, with the same model as the
baseline (`both_sides` profile, plate stacking, bare cylinder R = 70 cm, H = 2R, default enrichment, 10 000 particles × 100 batches, 40 inactive):

* **Common depth scan:** fuel and coolant slot depth together, `d = 0.6, 0.8, 1.0, 1.2` cm (`w = 2d + flat`; the pitch shrinks with `d`).
* **Coolant-only depth scan:** fuel depth fixed at 1.0 cm, `coolant_depth = 0.6 … 1.2` cm (coolant width `w_c = 2 d_c + flat`).
  With `coolant_depth ≠ slot_depth` the unit cell stays consistent: the coolant layer is `d_c` thick
  (`P_x = d + d_c + 2·wall`), fuel slots repeat along y with `P_y = w + web`, coolant slots along z with `P_z = w_c + web`
  (one web thickness everywhere). `coolant_depth=None` (default) = identical slots.

The scans were run with `dev/depth_scan.py common|coolant`; this cell **loads** `results/depth_scan.csv` and
`results/coolant_depth_scan.csv` (set `RUN_DEPTH_SCAN = True` to recompute them here, ~10 min each).
Core fuel volume = fuel volume fraction × cylinder volume (2155 L); fuel mass uses the MSRE fuel-salt density at 922 K.''')
code(r'''
RUN_DEPTH_SCAN = False
SCAN_DEPTHS = [0.6, 0.8, 1.0, 1.2]
DS_FIXED = dict(flat_width=0.5, web_thickness=1.5, wall_thickness=0.5)
SCANS = {  # name: (csv, png, x column, geometry for depth x)
    "common":  ("depth_scan.csv", "depth_scan.png", "slot_depth", lambda x: dict(slot_depth=x, **DS_FIXED)),
    "coolant": ("coolant_depth_scan.csv", "coolant_depth_scan.png", "coolant_depth",
                lambda x: dict(slot_depth=1.0, coolant_depth=x, **DS_FIXED)),
}

def _run_k(m, cwd):
    sp_path = m.run(cwd=cwd, output=False, threads=THREADS)
    with openmc.StatePoint(sp_path) as sp:
        gt = sp.global_tallies
        names = [n.decode() if isinstance(n, bytes) else str(n) for n in gt["name"]]
        return sp.keff.nominal_value, sp.keff.std_dev, float(gt["mean"][names.index("leakage")])

def run_depth_scan(name):
    csv, _, xcol, geo_of = SCANS[name]
    rho_f = make_materials(temperature=TEMPERATURE_K, enrichment=ENRICHMENT)["fuel"].density
    rows = []
    for x_ in SCAN_DEPTHS:
        geo = geo_of(x_)
        am = analytic_metrics(core_radius=CORE_RADIUS, **geo, **GEOM_OPTS)
        base = os.path.join(WORK_DIR, "depth_scan" if name == "common" else "coolant_depth_scan", f"d{x_:.1f}")
        ki = _run_k(build_model(**geo, **GEOM_OPTS, mode="unit_cell", enrichment=ENRICHMENT, **RUN_OPTS), base + "_kinf")
        ke = _run_k(build_model(**geo, **GEOM_OPTS, mode="cylinder", core_radius=CORE_RADIUS, enrichment=ENRICHMENT, **RUN_OPTS), base + "_cyl")
        rows.append(dict(slot_depth=geo["slot_depth"], coolant_depth=am["coolant_depth"], **DS_FIXED, slot_width=am["slot_width"],
                         coolant_slot_width=am["coolant_slot_width"], kinf=ki[0], kinf_std=ki[1], keff=ke[0], keff_std=ke[1],
                         leakage_fraction=ke[2], fuel_vf=am["fuel_vf"], coolant_vf=am["coolant_vf"], graphite_vf=am["graphite_vf"],
                         graphite_to_fuel=am["graphite_to_fuel"], core_fuel_volume_l=am["core_fuel_volume_l"],
                         core_fuel_mass_kg=am["core_fuel_volume_l"] * rho_f))
    os.makedirs(os.path.join(WORK_DIR, "results"), exist_ok=True)
    pd.DataFrame(rows).to_csv(os.path.join(WORK_DIR, "results", csv), index=False)

def plot_depth_scan(ds, xcol, title, png):
    """phone-friendly: three stacked panels, large fonts"""
    with plt.rc_context({"font.size": 13}):
        fig, (a1, a2, a3) = plt.subplots(3, 1, figsize=(6, 10.5), sharex=True, gridspec_kw=dict(height_ratios=[1, 1, 0.8]))
        a1.errorbar(ds[xcol], ds.kinf, yerr=ds.kinf_std, fmt="o-", color="tab:red", capsize=4, lw=2, ms=7)
        a1.set_ylabel("unit-cell k-inf"); a1.grid(alpha=0.3); a1.set_title(title, fontsize=12)
        a2.errorbar(ds[xcol], ds.keff, yerr=ds.keff_std, fmt="s-", color="tab:blue", capsize=4, lw=2, ms=7)
        a2.set_ylabel(f"bare-cylinder k-eff\n(R = {CORE_RADIUS:.0f} cm, H = 2R)"); a2.grid(alpha=0.3)
        for x_, k_, L_ in zip(ds[xcol], ds.keff, ds.leakage_fraction):
            a2.annotate(f"leak {L_:.3f}", (x_, k_), xytext=(0, 9), textcoords="offset points", ha="center", fontsize=10)
        a2.margins(y=0.25)
        a3.plot(ds[xcol], ds.graphite_to_fuel, "D-", color="0.3", lw=2, ms=7)
        a3.set_ylabel("C / fuel ratio"); a3.grid(alpha=0.3)
        b_ = a3.twinx()
        b_.plot(ds[xcol], ds.fuel_vf, "^--", color="tab:orange", lw=1.5, ms=6, label="fuel")
        b_.plot(ds[xcol], ds.coolant_vf, "v:", color="tab:cyan", lw=1.5, ms=6, label="coolant")
        b_.set_ylabel("salt vol. fraction"); b_.legend(fontsize=9, loc="best")
        a3.set_xlabel(("slot depth d" if xcol == "slot_depth" else "coolant slot depth d_c") + " [cm]"); a3.set_xticks(ds[xcol])
        fig.tight_layout(); os.makedirs(os.path.dirname(png), exist_ok=True); fig.savefig(png, dpi=150); plt.show()

TITLES = {"common": "Slot depth scan (fuel = coolant slots)\nflat 0.5, web 1.5, wall 0.5 cm; w = 2d + flat",
          "coolant": "Coolant-only depth scan (fuel depth 1.0 cm)\nflat 0.5, web 1.5, wall 0.5 cm; w_c = 2d_c + flat"}
SHOW = ["slot_depth", "coolant_depth", "slot_width", "coolant_slot_width", "kinf", "kinf_std", "keff", "keff_std", "leakage_fraction",
        "fuel_vf", "coolant_vf", "graphite_vf", "graphite_to_fuel", "core_fuel_volume_l", "core_fuel_mass_kg"]
depth_scans = {}
for name, (csv, png, xcol, _) in SCANS.items():
    path = os.path.join(WORK_DIR, "results", csv)
    if RUN_DEPTH_SCAN and HAVE_XS:
        run_depth_scan(name)
    if not os.path.exists(path):
        print(f"[{name}] no results at {path} - run dev/depth_scan.py {name} or set RUN_DEPTH_SCAN = True"); continue
    ds = depth_scans[name] = pd.read_csv(path)
    print(f"\n{name} depth scan:"); display(ds[[c for c in SHOW if c in ds]].round(5))
    plot_depth_scan(ds, xcol, TITLES[name], os.path.join(WORK_DIR, "figures", png))
''')

md('''## Notes / next steps
* Raise `PARTICLES`/`BATCHES` (e.g. 20 000 × 150) and `N_SAMPLES` for real studies; each unit-cell run is independent → trivially parallel.
* Useful extra outputs to add per sample: fuel/coolant temperature coefficients (re-run at ±ΔT), conversion ratio (tally U-238 capture / U-235 absorption), spectrum.
* `mode="unit_cell"` gives the quick periodic k-inf model; `reflector_thickness` adds graphite around the cylinder.
* All slots are identical by design (Brian); the profile builder takes width/depth per call if that ever needs to change.
''')

nb = nbf.v4.new_notebook(); nb.cells = cells
nb.metadata["kernelspec"] = {"name": "python3", "display_name": "Python 3 (OpenMC)", "language": "python"}
nbf.write(nb, "/workspace/msr_slab/msr_graphite_slab_lhs.ipynb")
print("cells:", len(cells))
