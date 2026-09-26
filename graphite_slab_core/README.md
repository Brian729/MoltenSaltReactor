# Graphite slab MSR core: OpenMC model with Latin hypercube sampling

`msr_graphite_slab_lhs.ipynb` builds a parameterized OpenMC model of a graphite-slab molten salt reactor core:

- Graphite slabs with vertical fuel slots alternating (through the slab thickness) with horizontal coolant slots, separated by graphite webs; slabs are isolated by a graphite wall.
- Both slot types share one machined profile: width x depth, flat at the mouth, rounded on one side with radius = slot depth (`round_location` = `side` | `bottom` | `none`).
- Slabs fill a cylinder of radius R and height 2R (bare by default; optional graphite reflector). A periodic unit-cell mode gives k-inf.
- Fuel: MSRE U-235 fuel salt, 7LiF-BeF2-ZrF4-UF4 65.0-29.17-5.0-0.83 mol% (ORNL-4658). Enrichment is the `ENRICHMENT` variable (U-235 wt%).
- Coolant: 7LiF-BeF2 66-34. Graphite 1.87 g/cc.
- Latin hypercube sampling (scipy `qmc.LatinHypercube`) over slot width, slot depth, web thickness, wall thickness; enrichment and R are optional extra dimensions. Set `RUN_TRANSPORT = True` to run the sweep.

Baseline (R = 70 cm, bare, ENDF/B-VIII.0, 10k particles x 100 batches): k-eff = 0.81888 +/- 0.00136; unit-cell k-inf = 1.63091 +/- 0.00107.

Requires OpenMC (tested with 0.16.0 from conda-forge) and `OPENMC_CROSS_SECTIONS` pointing to an ENDF/B-VIII.0 HDF5 library that includes `c_Graphite`.

`figures/` holds geometry plots, `results/` the baseline and LHS tables, `dev/` the scripts that generate the notebook.
