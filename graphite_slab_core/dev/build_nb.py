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
  (default **70 cm**, ~MSRE core radius) and z-planes at ±R (**H = 2R = D**), vacuum boundary. Optional graphite reflector:
  `reflector_thickness` (radial) and `reflector_axial` (top and bottom; default = radial) — see §11. No plenums, vessel, downcomer
  or headers — **the fuel salt outside the core is not modelled**. `mode="unit_cell"` = periodic unit cell → k-infinity (quick mode).
* Materials at 922 K (MSRE operating temperature): **fuel = MSRE carrier salt with 4.0 mol% UF₄** (extra UF₄ taken from LiF),
  7LiF-BeF₂-ZrF₄-UF₄ **61.83-29.17-5.0-4.0 mol%** (variable `UF4_MOLPCT`), **HALEU 19.75 wt% U-235** (variable `ENRICHMENT`),
  99.995 % Li-7; ρ = MSRE correlation 2.575 − 5.13·10⁻⁴·T[°C] g/cc rescaled with additive molar volumes → 2.593 g/cc at 922 K (§2);
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
REFLECTOR_THICKNESS = 0.0  # cm of graphite around the cylinder (side; top/bottom too unless reflector_axial= is given); 0 = bare
THREADS = None             # None = all OpenMP threads

# ---- default geometry (cm) --------------------------------------------------------------------------------
DEFAULTS = dict(slot_depth=1.0, flat_width=0.5, web_thickness=1.5, wall_thickness=0.5)   # -> slot width w = 2.5 cm
GEOM_OPTS = dict(round_location="both_sides", stacking="plates")                          # Brian's design (defaults)

# ---- fuel ------------------------------------------------------------------------------------------------
ENRICHMENT = 19.75         # U-235 WEIGHT % of total uranium. Default = HALEU (MSRE 235U operation was 33.477, ORNL-4658 Table 2.8)
UF4_MOLPCT = 4.0           # mol% UF4 in the fuel salt; the difference from MSRE's 0.83 mol% is taken from LiF (0.83 = MSRE salt)
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
## 2. Materials — default fuel = MSRE carrier salt with 4.0 mol% UF₄ and HALEU (19.75 wt% U-235)

**Default fuel (Brian, Sep 2026):** 7LiF-BeF₂-ZrF₄-UF₄ **61.83-29.17-5.0-4.0 mol%** — the MSRE fuel salt with UF₄ raised from 0.83 to
**4.0 mol%, the difference taken from LiF** (`UF4_MOLPCT`, `build_model(..., uf4_mol_pct=...)`; `fuel_composition(x)` gives
LiF = 65.0 − (x − 0.83)). No thorium. Uranium = **HALEU 19.75 wt% U-235** (`ENRICHMENT`).

* **HALEU isotopics** (`minor_u="correlation"`, default): U-234 = 0.0089 × U-235 wt% = **0.176 wt%** (typical for enrichment from
  natural feed; ASTM C996 caps U-234 at 1.1 × 10⁴ µg/g U-235, i.e. 0.011 × e), **U-236 = 0** (fresh, non-recycled feed; set
  `fuel={"u236_trace_wt_pct": ...}` for downblended/recycled material), **U-238 = balance 80.074 wt%**.
  MSRE uranium is still available: `fuel=MSRE_ISOTOPICS` (U-234 0.342, U-235 33.477, U-236 0.141 wt%).
* **Density with 4 mol% UF₄:** the MSRE correlation (below) belongs to the 0.83 mol% salt, so it is rescaled by the ratio of
  (molar mass / molar volume) of the new and the MSRE composition, using **additive molar volumes** from S. Cantor,
  *Density and viscosity of several molten fluoride mixtures*, ORNL-TM-4308 (1973): LiF 13.24/13.77, BeF₂ 24.0/24.2,
  ZrF₄ 46/48, UF₄ 45.1/46.1 cm³/mol at 550/700 °C (linear in T). Cantor found additive volumes within ~2 % of measured molar
  volumes for LiF-BeF₂-(Zr,Th,U)F₄ melts; the ratio form keeps the measured MSRE density exact at 0.83 mol%.
  Result at 922 K: **ρ = 2.593 g/cm³** (vs 2.242 for the 0.83 mol% salt) → U-235 density 0.0960 g/cm³ of salt
  (MSRE: 0.0355). Uncertainty of the estimate ~±2 %. Liquidus/solubility of 4 mol% UF₄ in this carrier is **not** checked here.

| Quantity | Value used | Source |
|---|---|---|
| Reference fuel salt (MSRE) | 7LiF-BeF₂-ZrF₄-UF₄ **65.0-29.17-5.0-0.83 mol%** (`uf4_mol_pct=0.83`) | R. E. Thoma, *Chemical Aspects of MSRE Operations*, ORNL-4658 (1971), p. 10–11 (fuel for 235U operation; also Table 1.1: 65-29.2-5-0.83) |
| MSRE uranium isotopics (start of power operation, run 4-1) — option `MSRE_ISOTOPICS` | U-234 0.342, U-235 33.477, U-236 0.141, U-238 66.041 wt% | ORNL-4658 Table 2.8 |
| Li-7 in fuel carrier salt | **99.995 at%** (batch assays 99.994–99.996) | ORNL-4658 Table 2.11; same value used in the IRPhEP MSRE benchmark (Shen/Fratoni et al., PHYSOR 2020) |
| MSRE fuel density (reference for the rescaling) | **ρ = 2.575 − 5.13×10⁻⁴·T(°C) g/cm³** (±1 %) → 2.242 g/cm³ at 649 °C (139.9 lb/ft³ at 650 °C) | ORNL-4658 Table 8.2 (Cantor molar-volume method; Table 8.3) |
| Coolant salt | 7LiF-BeF₂ 66-34 mol%, 99.992 % Li-7, ρ = 2.214 − 4.2×10⁻⁴·T(°C) → 1.941 g/cm³ at 649 °C | ORNL-4658 Tables 2.1, 8.1; ORNL-4616 |
| Graphite | 1.87 g/cm³ (MSRE grade CGB), pure C + `c_Graphite` | IRPhEP MSRE benchmark evaluation (1.87 ± 0.02 g/cm³) |

**Enrichment variable:** `ENRICHMENT` (top cell) = U-235 **weight %** of total uranium; also `build_model(..., enrichment=...)`,
`make_materials(enrichment=...)` and an optional 5th LHS dimension (`LHS_BOUNDS["enrichment"]`, off by default).
With the default `minor_u="correlation"` U-234 follows 0.0089 × e and U-236 = 0; `minor_u="scale"` scales the MSRE U-234/U-236
(exact MSRE isotopics at 33.477 wt%); `fuel={"minor_u": "none"}` gives a pure U-235/U-238 vector.
The UF₄ mole fraction (`UF4_MOLPCT`) is held when enrichment changes (so changing enrichment changes the U-235 loading).

**Alternatives found in ORNL sources (not used by default; switch via the `fuel=` dict):**
* Nominal design composition 65-29.1-5-0.9 mol% (ORNL-TM-728, Robertson 1965, Table 2.1; ORNL-4616; Haubenreich & Engel, *Nucl. Appl. Tech.* 8 (1970)); U "about 32 %" (ORNL-4616) / "33 %" (Haubenreich & Engel) enriched.
* Zero-power first-criticality salt 64.88-29.27-5.06-0.79 mol%, 1.408 wt% U-235 in salt, ρ = 2.3275 ± 0.016 g/cm³ at 638 °C (IRPhEP MSRE benchmark).
* Density: in-reactor inventory estimate 139.01 lb/ft³ (2.227 g/cm³) at start of power operation; electrical-probe correlation ρ = 2.848 − 7.69×10⁻⁴·T(°C) (2.35 g/cm³ at 650 °C); design value 141 lb/ft³ (2.26 g/cm³) (ORNL-4658 §8.2, Table 8.3).
* Li-7: "at least 99.99 %" (ORNL-4616), 99.994 % (ORNL-TM-728 Table 2.1).
""")

code(chunks["2"] + '''
_m = make_materials(temperature=TEMPERATURE_K, enrichment=ENRICHMENT, uf4_mol_pct=UF4_MOLPCT)
print("fuel composition (mol%):", fuel_composition(UF4_MOLPCT))
for k, m in _m.items():
    print(f"{k:9s} {m.name:55s} {m.density:.4f} g/cc @ {m.temperature:.0f} K  nuclides: {', '.join(n.name for n in m.nuclides)}")
print("uranium vector (wt% of U):", {k: round(100 * v, 3) for k, v in uranium_wt_fractions(
      ENRICHMENT, DEFAULT_FUEL["u234_wt_pct"], DEFAULT_FUEL["u236_wt_pct"], DEFAULT_FUEL["minor_u"],
      DEFAULT_FUEL["u234_per_u235"], DEFAULT_FUEL["u236_trace_wt_pct"]).items()})
print(f"U-235 mass fraction in fuel salt: {100 * _m['fuel'].get_mass_density('U235') / _m['fuel'].density:.3f} wt%, "
      f"U-235 density {_m['fuel'].get_mass_density('U235'):.4f} g/cc")
print(f"reference MSRE salt (0.83 mol% UF4) density: {make_materials(uf4_mol_pct=0.83)['fuel'].density:.4f} g/cc")''')

md('''## 3. Parameter validation, layer stack and analytic volume fractions
`resolve_params()` validates the inputs (positive d/web/wall, flat_width ≥ 0), derives `slot_width = 2d + flat_width` and returns the x-layer stack;
`analytic_metrics()` gives exact volume fractions / moderator-to-fuel ratio for the periodic cell (no transport needed).''')
code(chunks["3"] + '''
pd.Series(analytic_metrics()).to_frame("default geometry")''')

md('''## 4. OpenMC geometry: `build_model(...) -> openmc.Model`''')
code(chunks["4"])

code(r'''
RUN_OPTS = dict(particles=PARTICLES, batches=BATCHES, inactive=INACTIVE, temperature=TEMPERATURE_K, uf4_mol_pct=UF4_MOLPCT)

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
               enrichment_u235_wt_pct=ENRICHMENT, uf4_mol_pct=UF4_MOLPCT, temperature_K=TEMPERATURE_K, **DEFAULTS,
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
             f"enrichment {ENRICHMENT} wt% U-235 unless varied, UF4 {UF4_MOLPCT} mol%, T = {TEMPERATURE_K} K, "
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
Two one-at-a-time scans, holding `flat_width = 0.5`, `web_thickness = 1.5`, `wall_thickness = 0.5` cm, with the same geometry model as the
baseline (`both_sides` profile, plate stacking, bare cylinder R = 70 cm, H = 2R, 10 000 particles × 100 batches, 40 inactive), but with the
**original MSRE fuel** (0.83 mol% UF₄, MSRE uranium 33.477 wt% U-235) — the scans were done before the switch to HALEU / 4 mol% UF₄:

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
DS_FUEL = dict(fuel=MSRE_ISOTOPICS, enrichment=MSRE_U235_WT_PCT)          # depth scans used the MSRE fuel ...
DS_RUN = dict(RUN_OPTS, uf4_mol_pct=0.83)                                 # ... with 0.83 mol% UF4
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
    rho_f = make_materials(temperature=TEMPERATURE_K, uf4_mol_pct=0.83, **DS_FUEL)["fuel"].density
    rows = []
    for x_ in SCAN_DEPTHS:
        geo = geo_of(x_)
        am = analytic_metrics(core_radius=CORE_RADIUS, **geo, **GEOM_OPTS)
        base = os.path.join(WORK_DIR, "depth_scan" if name == "common" else "coolant_depth_scan", f"d{x_:.1f}")
        ki = _run_k(build_model(**geo, **GEOM_OPTS, mode="unit_cell", **DS_FUEL, **DS_RUN), base + "_kinf")
        ke = _run_k(build_model(**geo, **GEOM_OPTS, mode="cylinder", core_radius=CORE_RADIUS, **DS_FUEL, **DS_RUN), base + "_cyl")
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

md(r'''## 10. HALEU critical radius (bare cylinder, H = 2R, no reflector)
Fuel: **7LiF-BeF₂-ZrF₄-UF₄ 61.83-29.17-5.0-4.0 mol%, HALEU 19.75 wt% U-235** (U-234 0.176 wt%, no U-236; ρ = 2.593 g/cc at 922 K, see §2).
Two geometries (d = 1.0, flat 0.5, web 1.5, wall 0.5 cm, plate stacking):
**A** = identical fuel/coolant slots; **B** = coolant slots only 0.6 cm deep (`coolant_depth=0.6`, w_c = 1.7 cm).

1. Unit-cell k-inf for A and B (the 0.83 mol% UF₄ HALEU values are kept as a reference point).
2. Bare-cylinder k-eff for R = 70 … 200 cm (H = 2R).
3. Critical radius from the one-group bare-core form `1/k_eff = a + b/(R + δ)²` (2-parameter fit with δ = 0, or 3-parameter if δ is physical
   and improves χ²); Rc solves `1/k = 1`; fit uncertainty by parametric Monte Carlo on the k-eff errors.
4. One confirmation run at the fitted Rc; the final Rc is corrected with the local slope dk/dR: `Rc = R_conf + (1 − k_conf)/(dk/dR)`,
   σ(Rc) = σ(k_conf)/(dk/dR).
5. Masses at Rc: fuel-salt volume = fuel volume fraction × πR²·2R (in-core salt only; no loops/plena), mass with ρ = 2.593 g/cc, U-235 from the
   material composition.

Run with `dev/haleu_scan.py kinf | rscan A B | confirm A=Rc B=Rc` and `dev/haleu_fit.py`; this cell loads `results/haleu_kinf.csv`,
`results/haleu_R_scan.csv`, `results/haleu_critical.csv` (set `RUN_HALEU_SCAN = True` to recompute here, ~30 min).''')
code(r'''
RUN_HALEU_SCAN = False
HALEU_CASES = {"A": dict(DEFAULTS, **GEOM_OPTS), "B": dict(DEFAULTS, **GEOM_OPTS, coolant_depth=0.6)}
HALEU_RADII = [70.0, 100.0, 130.0, 160.0, 200.0]
RES_DIR = os.path.join(WORK_DIR, "results")

def haleu_masses(case, R):
    am = analytic_metrics(core_radius=R, **HALEU_CASES[case])
    f = make_materials(temperature=TEMPERATURE_K, enrichment=ENRICHMENT, uf4_mol_pct=UF4_MOLPCT)["fuel"]
    Vf = am["core_fuel_volume_l"] * 1000
    return dict(core_volume_m3=am["core_volume_l"] / 1000, fuel_salt_volume_m3=Vf / 1e6, fuel_salt_mass_kg=Vf * f.density / 1000,
                u235_mass_kg=Vf * f.get_mass_density("U235") / 1000,
                u_mass_kg=Vf * sum(f.get_mass_density(n) for n in f.get_nuclides() if n.startswith("U")) / 1000)

if RUN_HALEU_SCAN and HAVE_XS:
    kr, rr = [], []
    for c, g in HALEU_CASES.items():
        k_, s_, _ = _run_k(build_model(**g, mode="unit_cell", enrichment=ENRICHMENT, **RUN_OPTS), os.path.join(WORK_DIR, "haleu", f"{c}_kinf"))
        kr.append(dict(case=c, uf4_mol_pct=UF4_MOLPCT, enrichment=ENRICHMENT, kinf=k_, kinf_std=s_))
        for R_ in HALEU_RADII:
            m_ = build_model(**g, mode="cylinder", core_radius=R_, enrichment=ENRICHMENT, **RUN_OPTS)
            k_, s_, L_ = _run_k(m_, os.path.join(WORK_DIR, "haleu", f"{c}_R{R_:.1f}"))
            rr.append(dict(case=c, kind="scan", uf4_mol_pct=UF4_MOLPCT, core_radius=R_, core_height=2 * R_, keff=k_, keff_std=s_,
                           leakage_fraction=L_, **haleu_masses(c, R_)))
    pd.DataFrame(kr).to_csv(os.path.join(RES_DIR, "haleu_kinf.csv"), index=False)
    pd.DataFrame(rr).to_csv(os.path.join(RES_DIR, "haleu_R_scan.csv"), index=False)
    # (confirmation runs at the fitted Rc: see dev/haleu_scan.py confirm)
''')
code(r'''
"""Critical-radius fit for the HALEU R scan + phone-friendly plot.
One-group bare-core form with H = 2R:  1/k_eff = a + b/(R + delta)^2   (a ~ 1/k_inf,eff, b ~ M^2 B^2 R^2 / k_inf).
Rc from 1/k = 1; uncertainty by parametric Monte Carlo on the k-eff statistical errors.
Writes results/haleu_critical.csv and figures/haleu_R_scan.png."""
import os, sys, json
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def inv_k(R, a, b, delta=0.0):
    return a + b / (R + delta) ** 2

def fit_case(df, n_mc=2000, seed=1):
    R, k, s = df.core_radius.values, df.keff.values, df.keff_std.values
    y, sy = 1 / k, s / k**2
    out = {}
    for name, f, p0 in (("2p", lambda R, a, b: inv_k(R, a, b), (0.8, 1e3)), ("3p", inv_k, (0.8, 1e3, 5.0))):
        if name == "3p" and len(R) < 4: continue
        try:
            p, cov = curve_fit(f, R, y, p0=p0, sigma=sy, absolute_sigma=True, maxfev=20000)
        except Exception:
            continue
        chi2 = float(np.sum(((f(R, *p) - y) / sy) ** 2)); dof = len(R) - len(p)
        rc = lambda pp: (np.sqrt(pp[1] / (1 - pp[0])) - (pp[2] if len(pp) > 2 else 0.0)) if pp[0] < 1 else np.nan
        rng = np.random.default_rng(seed); rcs = []
        for _ in range(n_mc):
            try:
                pp, _ = curve_fit(f, R, y + rng.normal(0, sy), p0=p, sigma=sy, maxfev=20000); rcs.append(rc(pp))
            except Exception:
                pass
        rcs = np.array(rcs); rcs = rcs[np.isfinite(rcs)]
        out[name] = dict(params=list(map(float, p)), chi2=chi2, dof=dof, Rc=float(rc(p)), Rc_std=float(rcs.std()) if len(rcs) else np.nan,
                         kinf_fit=1 / p[0])
    # prefer 3-parameter fit only if delta is physical (0..30 cm) and it lowers chi2/dof
    best = "2p"
    if "3p" in out and 0 <= out["3p"]["params"][2] <= 30 and out["3p"]["dof"] > 0 and \
            out["3p"]["chi2"] / out["3p"]["dof"] < out["2p"]["chi2"] / max(out["2p"]["dof"], 1):
        best = "3p"
    return out, best

if os.path.exists(os.path.join(WORK_DIR, "results", "haleu_R_scan.csv")):
    W = WORK_DIR
    scan = pd.read_csv(os.path.join(W, "results", "haleu_R_scan.csv"))
    kinf_all = pd.read_csv(os.path.join(W, "results", "haleu_kinf.csv"))
    UF4_SCAN = float(scan.uf4_mol_pct.iloc[0]) if "uf4_mol_pct" in scan else 4.0
    kinf = kinf_all[np.isclose(kinf_all.uf4_mol_pct, UF4_SCAN)].set_index("case")      # k-inf of the scanned salt
    rows, fits = [], {}
    for c, g in scan[scan.kind == "scan"].groupby("case"):
        g = g.sort_values("core_radius"); out, best = fit_case(g); fits[c] = (out, best)
        conf = scan[(scan.case == c) & (scan.kind == "confirm")]
        r = dict(case=c, fit=best, Rc_fit=out[best]["Rc"], Rc_fit_std=out[best]["Rc_std"], chi2=out[best]["chi2"], dof=out[best]["dof"],
                 fit_params=json.dumps(out[best]["params"]), Rc_2p=out["2p"]["Rc"], Rc_3p=out.get("3p", {}).get("Rc", np.nan),
                 kinf_unit_cell=kinf.loc[c, "kinf"], kinf_unit_cell_std=kinf.loc[c, "kinf_std"])
        if len(conf):
            cr = conf.iloc[-1]; p = out[best]["params"]
            # local slope dk/dR of the fit at the confirmation radius -> corrected Rc
            Rq = cr.core_radius; h = 0.5
            dkdR = (1 / inv_k(Rq + h, *p) - 1 / inv_k(Rq - h, *p)) / (2 * h)
            r.update(R_confirm=Rq, keff_confirm=cr.keff, keff_confirm_std=cr.keff_std, leakage_confirm=cr.leakage_fraction,
                     dkdR_per_cm=dkdR, Rc=Rq + (1 - cr.keff) / dkdR,
                     Rc_std=float(cr.keff_std / dkdR), H_c=2 * (Rq + (1 - cr.keff) / dkdR),
                     # volumes/masses scale with R^3 (H = 2R): evaluate at the corrected Rc
                     **{k: cr[k] * ((Rq + (1 - cr.keff) / dkdR) / Rq) ** 3
                        for k in ("core_volume_m3", "fuel_salt_volume_m3", "fuel_salt_mass_kg", "u235_mass_kg", "u_mass_kg")})
            # systematic: spread of Rc between fit forms and a local fit (R <= 100 cm) - lattice-edge granularity / fit form
            loc_out, _ = fit_case(g[g.core_radius <= 100], n_mc=200)
            cands = [out["2p"]["Rc"], out.get("3p", {}).get("Rc", np.nan)] + [v["Rc"] for v in loc_out.values()]
            cands = [x for x in cands if np.isfinite(x)]
            r.update(Rc_fit_spread=max(cands) - min(cands))
        rows.append(r)
    crit = pd.DataFrame(rows); crit.to_csv(os.path.join(W, "results", "haleu_critical.csv"), index=False)
    print("unit-cell k-inf (HALEU):"); display(kinf_all[["case", "uf4_mol_pct", "coolant_depth", "kinf", "kinf_std", "fuel_vf", "coolant_vf", "graphite_vf", "graphite_to_fuel"]].round(5))
    print("R scan:"); display(scan[["case", "kind", "core_radius", "keff", "keff_std", "leakage_fraction", "core_volume_m3", "fuel_salt_volume_m3", "fuel_salt_mass_kg", "u235_mass_kg"]].sort_values(["case", "core_radius"]).round(4))
    print("critical radius:"); display(crit.T)

    # ---- phone-friendly plot ----
    colors = {"A": "tab:blue", "B": "tab:red"}
    labels = {"A": "A: coolant depth 1.0 cm", "B": "B: coolant depth 0.6 cm"}
    with plt.rc_context({"font.size": 13}):
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 9.5), sharex=True, gridspec_kw=dict(height_ratios=[1.6, 1]))
        Rgrid = np.linspace(scan.core_radius.min() * 0.95, scan.core_radius.max() * 1.03, 200)
        for c, g in scan.groupby("case"):
            col = colors.get(c, "k"); s_ = g[g.kind == "scan"].sort_values("core_radius"); cf = g[g.kind == "confirm"]
            a1.errorbar(s_.core_radius, s_.keff, yerr=s_.keff_std, fmt="o", color=col, capsize=4, ms=7, label=labels.get(c, c))
            if c in fits:
                out, best = fits[c]; a1.plot(Rgrid, 1 / inv_k(Rgrid, *out[best]["params"]), "-", color=col, lw=1.5, alpha=0.8)
            if len(cf):
                a1.errorbar(cf.core_radius, cf.keff, yerr=cf.keff_std, fmt="*", color=col, ms=14, mec="k", capsize=4)
            if c in kinf.index:
                a1.axhline(kinf.loc[c, "kinf"], color=col, ls="--", lw=1.2)
                a1.text(Rgrid[0], kinf.loc[c, "kinf"], f" k-inf {c} = {kinf.loc[c, 'kinf']:.3f}", color=col, va="bottom", fontsize=10)
            a2.plot(s_.core_radius, s_.leakage_fraction, "o-", color=col, ms=6, lw=1.5)
        a1.axhline(1.0, color="k", lw=1)
        for _, r in crit.iterrows():
            Rc = r.get("Rc", r["Rc_fit"]); sd = r.get("Rc_std", r["Rc_fit_std"])
            a1.axvline(Rc, color=colors.get(r.case, "k"), ls=":", lw=1.2)
            a1.annotate(f"Rc {r.case} = {Rc:.1f} ± {sd:.1f} cm", (Rc, 1.0), xytext=(4, -18 if r.case == "A" else 8),
                        textcoords="offset points", color=colors.get(r.case, "k"), fontsize=10)
        a1.set_ylabel("bare-cylinder k-eff (H = 2R)"); a1.grid(alpha=0.3); a1.legend(fontsize=10, loc="lower right")
        a1.set_title(f"HALEU 19.75 wt% U-235, 7LiF-BeF2-ZrF4-UF4\n({100 - 34.17 - UF4_SCAN:.2f}-29.17-5.0-{UF4_SCAN:.1f} mol%), bare, H = 2R; ★ = check at Rc", fontsize=12)
        a2.set_ylabel("leakage fraction"); a2.set_xlabel("core radius R [cm]"); a2.grid(alpha=0.3)
        fig.tight_layout(); os.makedirs(os.path.join(W, "figures"), exist_ok=True)
        fig.savefig(os.path.join(W, "figures", "haleu_R_scan.png"), dpi=150); plt.show()

else:
    print("no HALEU results found in", os.path.join(WORK_DIR, "results"))

''')

md(r'''## 11. Reflector thickness (R = 85 cm, H = 2R)
Graphite reflector (same graphite as the core: 1.87 g/cc, `c_Graphite`, 922 K) of thickness **t on the side AND on top/bottom**
(`build_model(..., reflector_thickness=t, reflector_axial=None)`; `reflector_axial` sets a different top/bottom thickness,
default = radial). No vessel, downcomer or plena yet — vacuum outside the reflector. Default fuel/geometry: HALEU 19.75 wt%,
4.0 mol% UF₄, case A (d = 1.0, flat 0.5, web 1.5, wall 0.5 cm, identical slots). Statistics as before (10 000 × 100, 40 inactive).

* **Savings vs bare:** Δk = (k − k_bare)·10⁵ pcm and Δρ = (1/k_bare − 1/k)·10⁵ pcm; **marginal worth** = Δk/Δt between steps.
* **Equivalent bare radius** R_eq(t): the bare-core radius with the same k-eff (monotone PCHIP interpolation of the §10 bare scan, case A).
  R_eq − 85 cm is the "equivalent radius gain" at R = 85 cm; it is **not** the critical-radius reduction (reflector savings shrink with
  the core size), so it overestimates the savings at criticality.
* **Critical radius with t = 30 cm (direct search):** reflected runs at R = 45, 55, 65 cm plus one at the interpolated radius; Rc by
  linear interpolation of k-eff between the two bracketing radii (σ from the k-eff errors).

Run with `dev/reflector_scan.py scan` / `crit 30 R1 R2 …` and `dev/reflector_fit.py`; this cell loads `results/reflector_scan.csv`
(set `RUN_REFLECTOR_SCAN = True` to recompute the thickness scan here, ~6 min).''')
code(r'''
RUN_REFLECTOR_SCAN = False
REFL_R, REFL_T = 85.0, [0, 10, 20, 30, 45, 60]
if RUN_REFLECTOR_SCAN and HAVE_XS:
    rr = []
    for t_ in REFL_T:
        m_ = build_model(**DEFAULTS, **GEOM_OPTS, mode="cylinder", core_radius=REFL_R, reflector_thickness=t_, reflector_axial=t_,
                         enrichment=ENRICHMENT, **RUN_OPTS)
        k_, s_, L_ = _run_k(m_, os.path.join(WORK_DIR, "reflector", f"R{REFL_R:.1f}_t{float(t_):.1f}"))
        rr.append(dict(kind="scan", core_radius=REFL_R, core_height=2 * REFL_R, reflector_radial=t_, reflector_axial=t_,
                       keff=k_, keff_std=s_, leakage_fraction=L_))
    pd.DataFrame(rr).to_csv(os.path.join(WORK_DIR, "results", "reflector_scan.csv"), index=False)

# geometry check: xz slice of the reflected core (t = 30 cm)
_mr = build_model(**DEFAULTS, **GEOM_OPTS, mode="cylinder", core_radius=REFL_R, reflector_thickness=30.0, enrichment=ENRICHMENT, **RUN_OPTS)
_Wr = 2.1 * (REFL_R + 30.0)
plot_slice(_mr, "xz", origin=(0, 0, 0), width=(_Wr, _Wr), pixels=(1000, 1000),
           title=f"reflected core, xz slice (R = {REFL_R:.0f} cm, H = 2R, 30 cm graphite radially and axially)",
           filename=os.path.join(WORK_DIR, "figures", "geom_reflected_xz.png"))
plt.show()
''')
code(r'''
"""Reflector-scan analysis + phone-friendly plot.
Savings vs bare: dk = (k - k_bare) x 1e5 pcm and reactivity drho = (1/k_bare - 1/k) x 1e5 pcm; marginal worth = dk/dt between steps.
Equivalent bare radius R_eq(t): radius of the BARE core (monotone PCHIP interpolation of the HALEU bare R scan, case A) with the
same k-eff -> "equivalent radius gain" R_eq - 85 cm.  NOTE: this is NOT the critical-radius reduction (the savings shrink when
the core shrinks); the reflected critical radius comes only from the direct search.
Direct check: 'crit' runs at t = 30 cm, Rc by linear interpolation between the bracketing radii.
Writes results/reflector_summary.csv and figures/reflector_scan.png."""
import os
import numpy as np, pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
import matplotlib.pyplot as plt

W = WORK_DIR
if os.path.exists(os.path.join(W, "results", "reflector_scan.csv")):
    df = pd.read_csv(os.path.join(W, "results", "reflector_scan.csv"))
    sc = df[df.kind == "scan"].sort_values("reflector_radial").reset_index(drop=True)
    k0 = sc.loc[sc.reflector_radial == 0, "keff"].iloc[0]; s0 = sc.loc[sc.reflector_radial == 0, "keff_std"].iloc[0]
    sc["dk_pcm"] = (sc.keff - k0) * 1e5
    sc["dk_pcm_std"] = np.hypot(sc.keff_std, s0) * 1e5
    sc["drho_pcm"] = (1 / k0 - 1 / sc.keff) * 1e5
    sc["marginal_pcm_per_cm"] = np.r_[np.nan, np.diff(sc.keff) / np.diff(sc.reflector_radial) * 1e5]
    # equivalent bare radius from the bare HALEU scan (case A)
    hs = pd.read_csv(os.path.join(W, "results", "haleu_R_scan.csv"))
    b = hs[(hs.case == "A") & (hs.kind == "scan")].sort_values("core_radius")
    kb = PchipInterpolator(b.core_radius.values, b.keff.values)
    crit = pd.read_csv(os.path.join(W, "results", "haleu_critical.csv")).set_index("case")
    Rc_bare, Rc_bare_std = crit.loc["A", "Rc"], crit.loc["A", "Rc_std"]
    Rmin, Rmax = b.core_radius.min(), b.core_radius.max()
    def r_eq(k):
        return brentq(lambda R: kb(R) - k, Rmin, Rmax) if kb(Rmin) < k < kb(Rmax) else np.nan
    sc["R_equiv_bare"] = [r_eq(k) for k in sc.keff]
    sc["equiv_radius_gain_cm"] = sc.R_equiv_bare - sc.core_radius
    # direct critical search (kind='crit')
    cr = df[df.kind == "crit"].sort_values("core_radius")
    direct = []
    for t, g in cr.groupby("reflector_radial"):
        g = pd.concat([g, sc[sc.reflector_radial == t]]).sort_values("core_radius")   # include the R=85 point
        R_, k_, s_ = g.core_radius.values, g.keff.values, g.keff_std.values
        i = np.searchsorted(k_, 1.0)
        if 0 < i < len(k_):
            slope = (k_[i] - k_[i - 1]) / (R_[i] - R_[i - 1])
            Rc = R_[i - 1] + (1 - k_[i - 1]) / slope
            w = (Rc - R_[i - 1]) / (R_[i] - R_[i - 1])
            sk = np.hypot((1 - w) * s_[i - 1], w * s_[i])
            fvf = g.fuel_salt_volume_m3.iloc[0] / (np.pi * R_[0] ** 2 * 2 * R_[0] / 1e6)          # fuel volume fraction
            Vf = fvf * np.pi * Rc ** 2 * 2 * Rc / 1e6
            rho_f = g.fuel_salt_mass_kg.iloc[0] / g.fuel_salt_volume_m3.iloc[0]
            u5 = g.u235_mass_kg.iloc[0] / g.fuel_salt_volume_m3.iloc[0]
            direct.append(dict(reflector=t, Rc_direct=Rc, Rc_direct_std=sk / slope, dkdR_pcm_per_cm=slope * 1e5,
                               bracket=f"{R_[i - 1]:g}-{R_[i]:g}", savings_vs_bare_cm=Rc_bare - Rc, Rc_bare=Rc_bare, H_c=2 * Rc,
                               outer_diameter=2 * (Rc + t), outer_height=2 * (Rc + t),
                               core_volume_m3=np.pi * Rc ** 2 * 2 * Rc / 1e6, fuel_salt_volume_m3=Vf, fuel_salt_mass_kg=Vf * rho_f,
                               u235_mass_kg=Vf * u5,
                               reflector_graphite_m3=np.pi * ((Rc + t) ** 2 * 2 * (Rc + t) - Rc ** 2 * 2 * Rc) / 1e6))
    direct = pd.DataFrame(direct)
    sc.to_csv(os.path.join(W, "results", "reflector_summary.csv"), index=False)
    if len(direct):
        direct.to_csv(os.path.join(W, "results", "reflector_critical.csv"), index=False)

    with plt.rc_context({"font.size": 13}):
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 9.5), sharex=True, gridspec_kw=dict(height_ratios=[1.5, 1]))
        a1.errorbar(sc.reflector_radial, sc.keff, yerr=sc.keff_std, fmt="o-", color="tab:blue", capsize=4, ms=7, lw=2)
        for _, r in sc.iterrows():
            if r.reflector_radial > 0:
                a1.annotate(f"+{r.dk_pcm:.0f} pcm", (r.reflector_radial, r.keff), xytext=(0, -20), textcoords="offset points",
                            ha="center", fontsize=10)
        a1.axhline(1.0, color="k", lw=1)
        a1.set_ylabel("k-eff (R = 85 cm, H = 2R)"); a1.grid(alpha=0.3)
        a1.set_title("Graphite reflector scan (radial = axial = t)\nHALEU 19.75 %, 4 mol% UF4, case A geometry", fontsize=12)
        a1.margins(y=0.15)
        a2.plot(sc.reflector_radial, sc.leakage_fraction, "s-", color="tab:red", ms=6, lw=2, label="leakage fraction")
        a2.set_ylabel("leakage fraction", color="tab:red"); a2.tick_params(axis="y", colors="tab:red"); a2.grid(alpha=0.3)
        b2 = a2.twinx()
        b2.plot(sc.reflector_radial, sc.marginal_pcm_per_cm, "D--", color="0.3", ms=6, lw=1.5)
        b2.set_ylabel("marginal worth [pcm/cm]")
        a2.set_xlabel("reflector thickness t [cm]"); a2.set_xticks(sc.reflector_radial)
        if len(direct):
            d0 = direct.iloc[0]
            a1.text(0.02, 0.03, f"critical R with t = {d0.reflector:.0f} cm: {d0.Rc_direct:.1f} ± {d0.Rc_direct_std:.1f} cm\n"
                                f"(bare: {Rc_bare:.1f} cm; savings {d0.savings_vs_bare_cm:.1f} cm)",
                    transform=a1.transAxes, fontsize=10, va="bottom")
        fig.tight_layout(); os.makedirs(os.path.join(W, "figures"), exist_ok=True)
        fig.savefig(os.path.join(W, "figures", "reflector_scan.png"), dpi=150); plt.show()
    display(sc[["reflector_radial", "keff", "keff_std", "leakage_fraction", "dk_pcm", "drho_pcm", "marginal_pcm_per_cm",
              "R_equiv_bare", "equiv_radius_gain_cm"]].round(4))
    display(direct.T) if len(direct) else print("no direct critical search yet")

''')

md(r'''## 12. Thin fuel slots with fixed width 1.0 cm — fuel/graphite optimum
Thin fuel slots raise the stagnant-salt conduction power limit (∝ 1/d_f²). Here the **fuel slot total width is fixed at 1.0 cm**
(`both_sides` profile, so `d_f ≤ 0.5` and `flat = 1.0 − 2 d_f`; d_f = 0.5 is a half-round) and the fuel depth is scanned,
d_f = 0.20 … 0.50 cm. Coolant slots stay at the default (depth 1.0, flat 0.5 → 2.5 cm wide) via the new optional
`coolant_depth=` / `coolant_flat_width=` parameters; web 1.5, wall 0.5 cm; HALEU 19.75 wt%, 4 mol% UF₄.
The unit cell stays consistent: P_x = d_f + 1.0 + 2·0.5, P_y = 1.0 + web (fuel slots), P_z = 2.5 + web (coolant slots).

For each depth: unit-cell k-inf, graphite-to-fuel volume ratio, bare R = 85 cm (H = 2R) k-eff; fuel salt and U-235 in the core
(fuel fraction × core volume); relative conduction limit (1.0/d_f)² vs the 1.0 cm default. If k-inf peaks at an edge of the depth
range, a unit-cell web-thickness scan at d_f = 0.5 (web 0.5 … 5 cm) locates the k-inf optimum in C/fuel ratio
(note: a thicker web also lowers the coolant-salt fraction, which is a parasitic absorber).
(`settings.source_rejection_fraction` is lowered to 0.001 because the fuel can be < 5 % of the volume.)

Run with `dev/fuel_depth_w1_scan.py depth|web` and `dev/fuel_w1_fit.py`; this cell loads `results/fuel_depth_w1_scan.csv`
(+ `results/fuel_w1_web_scan.csv`).''')
code(r'''
"""Analysis + phone-friendly plot for the fixed-width (1.0 cm) fuel-slot depth scan and the optional web scan.
Optimum = maximum unit-cell k-inf; located with a parabola through the best point and its neighbours (if interior).
Writes results/fuel_w1_summary.csv and figures/fuel_depth_w1_scan.png."""
import os
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

W = WORK_DIR
p_d, p_w = os.path.join(W, "results", "fuel_depth_w1_scan.csv"), os.path.join(W, "results", "fuel_w1_web_scan.csv")

def optimum(x, k):
    """argmax of k(x): parabola through the best point and its neighbours; flags edge maxima."""
    x, k = np.asarray(x, float), np.asarray(k, float); i = int(np.argmax(k))
    if i in (0, len(x) - 1):
        return dict(x_opt=x[i], k_opt=k[i], at_edge=True)
    c = np.polyfit(x[i - 1:i + 2], k[i - 1:i + 2], 2)
    xo = -c[1] / (2 * c[0]) if c[0] < 0 else x[i]
    return dict(x_opt=xo, k_opt=np.polyval(c, xo), at_edge=False)

if os.path.exists(p_d):
    fd = pd.read_csv(p_d).sort_values("fuel_depth").reset_index(drop=True)
    wb = pd.read_csv(p_w).sort_values("web_thickness").reset_index(drop=True) if os.path.exists(p_w) else None
    fd["critical_bare_R85"] = fd.keff_bare - 2 * fd.keff_bare_std >= 1.0
    summ = [dict(scan="fuel depth (w_f = 1.0, web 1.5)", variable="fuel_depth", **optimum(fd.fuel_depth, fd.kinf))]
    summ[0]["graphite_to_fuel_at_opt"] = float(np.interp(summ[0]["x_opt"], fd.fuel_depth, fd.graphite_to_fuel))
    if wb is not None and len(wb) > 2:
        o = optimum(wb.web_thickness, wb.kinf)
        o["graphite_to_fuel_at_opt"] = float(np.interp(o["x_opt"], wb.web_thickness, wb.graphite_to_fuel))
        summ.append(dict(scan="web (d_f = 0.5, w_f = 1.0)", variable="web_thickness", **o))
    summ = pd.DataFrame(summ); summ.to_csv(os.path.join(W, "results", "fuel_w1_summary.csv"), index=False)

    n = 4 if wb is not None else 3
    with plt.rc_context({"font.size": 13}):
        fig, axs = plt.subplots(n, 1, figsize=(6, 3.4 * n + 0.8))
        a1, a2, a3 = axs[:3]
        a1.errorbar(fd.fuel_depth, fd.kinf, yerr=fd.kinf_std, fmt="o-", color="tab:purple", capsize=4, ms=7, lw=2)
        a1.set_ylabel("unit-cell k-inf"); a1.grid(alpha=0.3)
        b1 = a1.twinx(); b1.plot(fd.fuel_depth, fd.graphite_to_fuel, "D--", color="0.4", ms=5, lw=1.2)
        b1.set_ylabel("C / fuel (vol)", color="0.4")
        a1.set_title("Fuel slots 1.0 cm wide (flat = 1 − 2 d_f), coolant 1.0 × 2.5 cm\nHALEU 19.75 %, 4 mol% UF4, web 1.5, wall 0.5 cm", fontsize=12)
        a2.errorbar(fd.fuel_depth, fd.keff_bare, yerr=fd.keff_bare_std, fmt="s-", color="tab:blue", capsize=4, ms=7, lw=2)
        a2.axhline(1.0, color="k", lw=1); a2.set_ylabel("bare k-eff\nR = 85 cm, H = 2R"); a2.grid(alpha=0.3)
        for x_, k_, L_ in zip(fd.fuel_depth, fd.keff_bare, fd.leakage_bare):
            a2.annotate(f"{L_:.2f}", (x_, k_), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=9)
        a2.text(0.02, 0.9, "labels = leakage fraction", transform=a2.transAxes, fontsize=9)
        a3.plot(fd.fuel_depth, fd.u235_mass_kg, "o-", color="0.3", ms=6, lw=2); a3.set_ylabel("U-235 in core [kg]"); a3.grid(alpha=0.3)
        b3 = a3.twinx(); b3.plot(fd.fuel_depth, fd.rel_conduction_power_limit, "^--", color="tab:orange", ms=6, lw=1.5)
        b3.set_ylabel("rel. conduction limit\n(1.0/d_f)²", color="tab:orange"); b3.tick_params(axis="y", colors="tab:orange")
        for a in (a1, a2, a3):
            a.set_xticks(fd.fuel_depth); a.set_xlabel("fuel slot depth d_f [cm]", fontsize=11)
        if wb is not None:
            a4 = axs[3]
            a4.errorbar(fd.graphite_to_fuel, fd.kinf, yerr=fd.kinf_std, fmt="o-", color="tab:purple", capsize=3, ms=6, label="depth scan (web 1.5)")
            a4.errorbar(wb.graphite_to_fuel, wb.kinf, yerr=wb.kinf_std, fmt="s-", color="tab:green", capsize=3, ms=6, label="web scan (d_f 0.5)")
            for x_, k_, t_ in zip(wb.graphite_to_fuel, wb.kinf, wb.web_thickness):
                a4.annotate(f"web {t_:g}", (x_, k_), xytext=(4, -12), textcoords="offset points", fontsize=8, color="tab:green")
            a4.set_xlabel("graphite / fuel volume ratio"); a4.set_ylabel("unit-cell k-inf"); a4.grid(alpha=0.3); a4.legend(fontsize=9)
        fig.tight_layout(); os.makedirs(os.path.join(W, "figures"), exist_ok=True)
        fig.savefig(os.path.join(W, "figures", "fuel_depth_w1_scan.png"), dpi=150); plt.show()
    display(fd[["fuel_depth", "fuel_flat", "kinf", "kinf_std", "graphite_to_fuel", "fuel_vf", "coolant_vf", "keff_bare", "keff_bare_std",
              "leakage_bare", "critical_bare_R85", "fuel_salt_volume_m3", "fuel_salt_mass_kg", "u235_mass_kg", "rel_conduction_power_limit"]].round(4))
    if wb is not None:
        display(wb[["web_thickness", "pitch_y", "pitch_z", "kinf", "kinf_std", "graphite_to_fuel", "fuel_vf", "coolant_vf"]].round(4))
    display(summ.round(4))
else:
    print("no fixed-width fuel-depth results in", os.path.join(WORK_DIR, "results"))

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
