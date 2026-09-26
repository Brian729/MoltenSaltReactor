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

Requires OpenMC (tested with 0.16.0 from conda-forge) and `OPENMC_CROSS_SECTIONS` pointing to an ENDF/B-VIII.0 HDF5 library that includes `c_Graphite`.

Folders:
- `figures/`: geometry plots, including `profile_sketch.png`.
- `results/`: the baseline and LHS tables.
- `dev/`: the scripts that generate the notebook.
