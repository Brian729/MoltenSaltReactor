# Graphite plate MSR core: OpenMC model with Latin hypercube sampling

`msr_graphite_slab_lhs.ipynb` builds a parameterized OpenMC model of a slotted-graphite molten salt reactor core:

- Slots are machined into one face of a graphite plate. Plates are stacked with every other plate rotated 90 degrees, so rows of vertical fuel slots alternate with rows of horizontal coolant slots. The solid backing of each plate is the single shared wall between a fuel row and the next coolant row (`stacking="plates"`: F | wall | C | wall).
- All slots (fuel and coolant) are identical. Default profile `round_location="both_sides"`: depth d, both side walls are quarter-rounds of radius d, flat bottom in between, so slot width w = 2d + flat_width. One web thickness everywhere (the land between slots in a row). Legacy profiles (`side`, `bottom`, `none`) and `stacking="interleaved"` are still available.
- Defaults: d = 1.0, flat = 0.5 (w = 2.5), web = 1.5, wall = 0.5 cm. The notebook includes a rule-of-thumb note on wall thickness (5 mm start, about 3 mm floor; pressure is not the limit).
- The plates fill a bare cylinder of radius R and height 2R (no reflector). A periodic unit-cell mode gives k-inf.
- Fuel: MSRE U-235 fuel salt, 7LiF-BeF2-ZrF4-UF4 65.0-29.17-5.0-0.83 mol% (ORNL-4658). Enrichment is set by the `ENRICHMENT` variable (U-235 wt%).
- Coolant: 7LiF-BeF2 66-34. Graphite: 1.87 g/cc.
- Latin hypercube sampling (scipy `qmc.LatinHypercube`) over slot depth, flat width, web thickness and wall thickness; the width is derived. Enrichment and R are optional extra dimensions. Set `RUN_TRANSPORT = True` to run the sweep.

Baseline (default geometry, R = 70 cm, bare, ENDF/B-VIII.0, 10k particles x 100 batches with 40 inactive):
- k-eff = 0.84793 +/- 0.00117, with leakage fraction 0.442.
- Unit-cell k-inf = 1.58746 +/- 0.00111.

Slot depth scans (notebook section 9, `dev/depth_scan.py common|coolant`; flat 0.5, web 1.5, wall 0.5 cm; same statistics as the baseline):
- **Common depth** (fuel and coolant slots together, d = 0.6 to 1.2 cm):
  - Shallower slots raise unit-cell k-inf (1.572 to 1.618).
  - They lower bare-core k-eff (0.851 to 0.823), because leakage rises (0.436 to 0.471).
  - Results: `results/depth_scan.csv`, plot: `figures/depth_scan.png`.
- **Coolant-only depth** (fuel depth 1.0 cm, coolant depth d_c = 0.6 to 1.2 cm, optional `coolant_depth=` parameter):
  - Shallower coolant slots raise both k-inf and k-eff (k-eff 0.890 at d_c = 0.6 vs 0.848 at 1.0).
  - Results: `results/coolant_depth_scan.csv`, plot: `figures/coolant_depth_scan.png`.

Requires OpenMC (tested with 0.16.0 from conda-forge) and `OPENMC_CROSS_SECTIONS` pointing to an ENDF/B-VIII.0 HDF5 library that includes `c_Graphite`.

Folders:
- `figures/`: geometry plots, including `profile_sketch.png`.
- `results/`: the baseline and LHS tables.
- `dev/`: the scripts that generate the notebook.
