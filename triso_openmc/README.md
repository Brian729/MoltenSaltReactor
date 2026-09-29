# TRISO fuel in OpenMC: starter model and references

This folder starts the move from liquid-fuel slab cores to **TRISO solid fuel in a FLiBe-cooled graphite core**
(a KP-FHR / gFHR-type pebble bed).

| File | What it is |
|---|---|
| `triso_pebble_flibe.ipynb` | Executed notebook: explicit, randomly packed TRISO in a gFHR pebble, BCC pebble cell, FLiBe coolant, reflective BCs → k-inf |
| `make_nb.py` | Script that regenerates the notebook (run it, then `jupyter nbconvert --execute`) |
| `figures/triso_pebble_geometry.png` | Geometry plot (phone-friendly, portrait): the whole cell plus a zoom on the particles |

## Starter model results (OpenMC 0.16.0, ENDF/B-VIII.0, 900 K, fresh fuel)

| Case | TRISO per pebble | Fuel-zone packing | ≈ g U per pebble | **k-inf** | Wall time (8 threads) |
|---|---|---|---|---|---|
| Baseline, gFHR particle count | 11,660 | 28.4 % | 4.37 | **1.38258 ± 0.00149** | 271 s |
| 40 % packing variant | 16,403 | 40.0 % | 6.15 | **1.30663 ± 0.00114** | 352 s |

Each case ran 10,000 particles × 70 batches (20 inactive). Packing the spheres took 4 s at 28 % and 43 s at 40 %.
k-inf *falls* as the TRISO loading rises, so this fresh FLiBe pebble lattice is under-moderated at these loadings.
That points toward lower packing or smaller kernels if you want the most reactivity per gram of U.

The model:
* **TRISO**: UC0.5O1.5 kernel, 425 µm diameter, 10.5 g/cc, **19.75 % U-235**. Layers: buffer 100 µm (1.05 g/cc),
  IPyC 40 µm (1.90), SiC 35 µm (3.18), OPyC 40 µm (1.90). Outer radius 427.5 µm.
* **Pebble**: 4 cm diameter. Low-density graphite core r = 1.38 cm (1.41 g/cc), fuel zone to 1.80 cm (matrix 1.74 g/cc),
  shell to 2.0 cm (1.74 g/cc).
* **Lattice**: one pebble at the centre of a cube plus 8 corner eighths (BCC). The cube side (4.816 cm) gives 60 %
  pebble packing, with FLiBe (Li-7 99.995 %, 1.973 g/cc at 900 K) filling the gaps. All six faces are reflective.
* **Source of dimensions and densities**: the gFHR benchmark as published on the INL Virtual Test Bed
  ([reactor description, Tables 1–3](https://mooseframework.inl.gov/virtual_test_bed/pbfhr/g_fhr/reactor_description.html)).
  The TRISO radii match the OECD/NEA MHTGR-350 benchmark, NEA/NSC/R(2017)4, Table I.4.
  The one deliberate deviation is the enrichment: 19.75 % here vs 19.55 % in the gFHR spec.

Caveats: this is an infinite, isothermal, fresh-fuel lattice, so it is **not** comparable to the gFHR equilibrium-core k-eff.
The pebble arrangement is regular (BCC), not random, and each case uses one random TRISO realisation.
ENDF/B-VIII.0 has no FLiBe S(α,β) data, so FLiBe is treated as a free gas. The gFHR sources disagree on particles per
pebble: VTB gives 11,660, while ANL's SAM paper gives 9,022. The VTB core table's "3.4 g IHM" per pebble matches
9,022 particles (9,022 × 3.75e-4 g U ≈ 3.4 g), so 9,022 may be the intended value. Re-run with `n_triso=9022` to check.

## Reference models, ranked by usefulness as a starting point for a FLiBe-cooled TRISO core

All links were checked on 2026-09-29. The "OpenMC" column gives the version each model was built or used with.

1. **gFHR pebble-bed benchmark: OpenMC input files** (Rizwan Ali et al., Zenodo, 2026-09-02, CC-BY-4.0).
   <https://zenodo.org/records/22251994>
   Full-core FLiBe-cooled pebble bed (gFHR) in five double-heterogeneity representations: fully stochastic,
   HEX-HEX, HCP-HCP, HEX-stochastic and HCP-stochastic. Built with **OpenMC 0.15.3**. Files:
   `Geometry_Representations.zip`, `Spectrum_power_flux_tallies.zip`.
   *This is the closest public match to KP-FHR.* The record was confirmed through the Zenodo API. The HTML page blocks
   bots, and I have not downloaded or run the inputs.
2. **Cardinal `scripts/openmc_pebble_ped_model.py`** (ANL/INL Cardinal, see the repo LICENSE).
   <https://github.com/neams-th-coe/cardinal/blob/devel/scripts/openmc_pebble_ped_model.py>
   Builds an OpenMC FLiBe pebble-bed model (1.5 cm Mk1 PB-FHR-style pebbles, 19.9 % U, Li-7 99.995 %,
   optional random TRISO) from a file of pebble centres. Uses the current `openmc.model` API.
   *Best existing OpenMC code for scaling up from a pebble to a bed.*
3. **gFHR specification on the INL Virtual Test Bed** (the VTB repo is CC-BY-4.0).
   <https://mooseframework.inl.gov/virtual_test_bed/pbfhr/g_fhr/reactor_description.html> ·
   models: <https://github.com/idaholab/virtual_test_bed/tree/devel/pbfhr/gFHR>
   Public tables for the core, pebble and TRISO that this starter model uses. The VTB gFHR and Mk1 models themselves
   are Griffin/Pronghorn (MOOSE), **not OpenMC**, so treat them as a spec and multiphysics reference.
4. **OpenMC "Modeling TRISO Particles" notebook** (openmc-dev/openmc-notebooks; no license file in the repo).
   <https://github.com/openmc-dev/openmc-notebooks/blob/main/triso.ipynb>
   A 1 cm reflective box at 30 % packing that shows `pack_spheres`, `TRISO` and `create_triso_lattice`.
   It uses the current API and runs on 0.16.
   The rendered docs page exists only for old versions, e.g. <https://docs.openmc.org/en/v0.10.0/examples/triso.html>,
   which uses the old `pack_trisos` API. API docs:
   [pack_spheres](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.model.pack_spheres.html),
   [TRISO](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.model.TRISO.html),
   [create_triso_lattice](https://docs.openmc.org/en/stable/pythonapi/generated/openmc.model.create_triso_lattice.html).
5. **UIUC ARFC OpenMC model of the OECD-NEA FHR (AHTR plate-fuel) benchmark** (BSD-3, OpenMC 0.12-dev).
   <https://github.com/pep8speaks/fhr-benchmark>
   A fork of the original `gwenchee/fhr-benchmark`, which is no longer online. Covers Phase I-A/B, FLiBe-cooled
   **plate-type** TRISO fuel element cases 1a–7a, plus depletion.
   Spec: <https://www.oecd-nea.org/jcms/pl_58295/benchmark-specifications-for-the-fluoride-salt-high-temperature-reactor-fhr-reactor-physics-calculations>
   (the full specification must be requested from the NEA). Useful if a plate or prismatic FHR is on the table.
   Expect minor API updates for 0.16.
6. **Cardinal `gas_compact` tutorial** (a TRISO compact unit cell) and **VTB `htgr/assembly`** (a prismatic assembly). Both
   use OpenMC with `pack_spheres` and `create_triso_lattice` and are helium-cooled.
   <https://github.com/neams-th-coe/cardinal/blob/devel/tutorials/gas_compact/unit_cell.py> ·
   <https://cardinal.cels.anl.gov/tutorials/gas_compact.html> ·
   <https://cardinal.cels.anl.gov/tutorials/gas_assembly.html> ·
   <https://github.com/idaholab/virtual_test_bed/tree/devel/htgr/assembly>
   These are clean, current-API **prismatic** templates (hex periodic BCs, axial layers, TRISO lattice) and are the best
   template for a prismatic FHR if you swap He for FLiBe in the channels.
   For gas-cooled microreactors (gcmr), the VTB `microreactors/gcmr` models are Griffin-based, not OpenMC; Cardinal's
   `gas_assembly` is the OpenMC analogue.
7. **MIT CRPG `openmc-reactor-examples` VHTR notebook** (no license stated).
   <https://github.com/mit-crpg/openmc-reactor-examples>
   A prismatic VHTR template with TRISO compacts at 30 % packing and a 500 µm kernel. A simple teaching model.
8. **Benchmark specs with no public OpenMC model found**: the MHTGR-350 OECD benchmark
   (<https://www.oecd-nea.org/upload/docs/application/pdf/2020-01/dir1/nsc-r2017-4.pdf>; the search index has it, but a
   direct download from the box timed out) and HTR-10 / HTTR / PBMR-400. VTB `htgr/htr10`, `httr` and `pbmr400` are
   Griffin models. For HTR-10 in OpenMC there are user-built models in an OpenMC forum thread
   (<https://openmc.discourse.group/t/htr-10-neutron-flux/6150>) and a BATAN journal paper
   (<https://jurnal.batan.go.id/index.php/jpen/article/download/6104/5331>, found by search; the fetch timed out).
   `mit-crpg/benchmarks` has no TRISO benchmarks.

## Recommendation

For a KP-FHR-like core, start from **this notebook**, which gives a verified pebble with current API, data and HALEU.
Then build out using **#1, the gFHR OpenMC inputs**: these are FLiBe pebble beds with a published verification, and they
let you compare double-heterogeneity treatments. Use **#2, the Cardinal pebble-bed script**, as a second working example
of turning pebble centres into an OpenMC bed. If you go prismatic or plate instead, start from **#6 (Cardinal/VTB compact
and assembly)** and **#5 (the UIUC OECD FHR plate benchmark)**.
