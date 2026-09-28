# 1D heat-path sweep vs fuel gap thickness (pad-supported flat-plate stack)

A quick steady 1D model of the heat path, with no CFD. It is for choosing the fuel gap thickness `g` in a flat-plate core.
In this core, graphite plates are held apart by a hex array of round graphite pads (1 cm diameter on a 3.5 cm pitch).

```
... | fuel gap g (stagnant fuel salt, q''') | graphite t_w | coolant gap g_c (flowing FLiBe) | graphite t_w | fuel gap g | ...
```

Each fuel gap is cooled through both faces, and each coolant gap is heated from both walls. By symmetry, the repeat unit is
`g + g_c + 2 t_w`.

Run: `python3 heat_path_1d.py` (needs numpy + matplotlib). It writes `results/*.csv` and `figures/*.png` in about 2 s.

## Equations

These are local values, taken away from the pads, with coolant bulk temperature `T_b` at the location.

| term | expression |
|---|---|
| heat flux per fuel-gap face | `q'' = q''' g / 2` |
| fuel conduction (uniform source, both faces cooled, **no convection**) | `ΔT_fuel = q''' g² / (8 k_f)` |
| graphite plate | `ΔT_gr = q'' t_w / k_g` |
| coolant film | `ΔT_film = q'' / h_c` |
| allowed peak source | `q'''_allowed = B / [ g²/(8 k_f) + (g/2)(t_w/k_g + 1/h_c) ]` |
| fuel volume fraction | `φ_f = g (1 − f_pad) / (g + g_c + 2 t_w)` |
| core power | `P = q'''_allowed / F_peak × V_core × φ_f` |

Coolant side (parallel plates, `D_h = 2 g_c`):

- Laminar (Re < 2300): `Nu = 8.235`. This is the fully developed value for both walls at uniform heat flux (Shah & London, *Laminar Flow Forced Convection in Ducts*, 1978). Friction is `f = 96/Re`.
- Turbulent (Re ≥ 3000): Gnielinski `Nu = (f/8)(Re−1000)Pr / (1 + 12.7 √(f/8)(Pr^{2/3}−1))`, with the Petukhov smooth-tube friction factor `f = (0.790 ln Re − 1.64)^−2`.
- 2300–3000: linear blend between the two.
- Pressure drop: `dP/L = f/D_h · ρv²/2`. A **×2** column gives a rough allowance for pad obstruction and form losses. The model has no pad pressure-drop term.

## Inputs

| quantity | value | source / note |
|---|---|---|
| t_w | 5 mm | design |
| k_g | 30 W/m·K (unirradiated), 15 W/m·K (irradiated) | typical nuclear graphite at about 600 °C |
| k_f (fuel salt) | 1.05 W/m·K | `openfoam_fuel_slot/README.md` |
| pads | 1 cm diameter, 3.5 cm hex pitch → `f_pad = (π/4 d²)/(√3/2 p²) = 0.0740` | design |
| V_core | 3.4 m³ (R = 81.5 cm, H = 2R) | reference slab-core size |
| F_peak (radial × axial) | 3.2 (bare, 2.3 × 1.4) and 1.5 (flattened) | assumed |
| budget B (coolant bulk → fuel centreline) | 50 K (fuel ≤ 700 °C at the 650 °C coolant hot end) and 100 K (coolant about 600 °C mid-core) | constraint |
| FLiBe at 600 °C (873 K) | ρ = 2413 − 0.488 T = 1987 kg/m³; μ = 1.16e-4 exp(3755/T) = 8.55 mPa·s; k = 1.1 W/m·K; c_p = 2386 J/kg·K; Pr = 18.5 | D.F. Williams et al., ORNL/TM-2006/12 (2006); M.R. Romatoski & L.W. Hu, "Fluoride salt coolant properties for nuclear reactor applications: A review", *Ann. Nucl. Energy* 109 (2017) 635–647 |

Sweep: g = 1, 1.5, 2, 3, 4, 5, 6, 8, 10 mm × g_c = 3, 5, 10 mm × v = 0.5, 1, 2 m/s × k_g × budget × peaking.

**Pad treatment.** Pads take 7.4% of the fuel volume (and of the plate area). The local ΔT budget is evaluated in the open
region between pads, so pads do not change `q'''_allowed`; they only reduce `φ_f`. Pads would also conduct some heat
directly from fuel-gap wall to wall, but as solid graphite in the fuel gap they carry no fuel. Their heat-path effect is
**ignored**. Ignoring it is conservative for the open-region peak, which is the location that sets the limit.

## Results

### Coolant hydraulics (`results/coolant_hydraulics.csv`)

| g_c (mm) | v (m/s) | Re | regime | h (W/m²K) | dP/L (kPa/m) | ×2 pads | dP across 1.63 m (kPa) |
|---|---|---|---|---|---|---|---|
| 3 | 0.5 | 697 | laminar | 1510 | 5.7 | 11.4 | 9.3 |
| 3 | 1 | 1394 | laminar | 1510 | 11.4 | 22.8 | 18.6 |
| 3 | 2 | 2787 | blend | 4450 | 29.4 | 58.8 | 47.9 |
| 5 | 0.5 | 1161 | laminar | 906 | 2.05 | 4.1 | 3.3 |
| 5 | 1 | 2323 | blend (≈laminar) | 989 | 4.16 | 8.3 | 6.8 |
| 5 | 2 | 4646 | turbulent | 5776 | 15.7 | 31.4 | 25.6 |
| 10 | 0.5 | 2323 | blend (≈laminar) | 494 | 0.52 | 1.0 | 0.85 |
| 10 | 1 | 4646 | turbulent | 2888 | 1.96 | 3.9 | 3.2 |
| 10 | 2 | 9292 | turbulent | 5824 | 6.39 | 12.8 | 10.4 |

In laminar flow, h does not depend on velocity and scales as 1/g_c. That is why g_c = 3 mm at 0.5 and 1 m/s gives
identical results. A narrow coolant gap gives the best laminar h, and the jump in h comes when flow reaches the turbulent
regime.

### Allowed q''' and core power vs g (50 K budget, k_g = 30, F_peak = 3.2)

Entries are `q'''_allowed (MW/m³) / P (MW)`. Full data is in `results/heat_path_sweep.csv`.

| g_c, v | g = 1 | 1.5 | 2 | 3 | 4 | 5 | 6 | 8 | 10 mm |
|---|---|---|---|---|---|---|---|---|---|
| 3 mm, ≤1 m/s | 94/6.6 | 56/5.7 | 38/5.0 | 22/4.0 | 14/3.3 | 9.9/2.7 | 7.4/2.3 | 4.6/1.7 | 3.1/1.3 |
| 3 mm, 2 m/s | 159/11.2 | 89/9.1 | 58/7.6 | 30/5.6 | 19/4.3 | 13/3.5 | 9.2/2.9 | 5.4/2.0 | 3.6/1.5 |
| 5 mm, 0.5 m/s | 66/4.1 | 41/3.7 | 29/3.3 | 17/2.8 | 11/2.3 | 8.1/2.0 | 6.2/1.7 | 3.9/1.4 | 2.7/1.1 |
| 5 mm, 1 m/s | 71/4.3 | 43/3.9 | 30/3.5 | 18/2.9 | 12/2.4 | 8.4/2.1 | 6.4/1.8 | 4.1/1.4 | 2.8/1.1 |
| 5 mm, 2 m/s | 173/10.6 | 96/8.6 | 61/7.1 | 32/5.2 | 19/4.0 | 13/3.2 | 9.4/2.7 | 5.6/1.9 | 3.7/1.5 |
| 10 mm, 0.5 m/s | 41/1.9 | 26/1.8 | 19/1.7 | 12/1.5 | 8.0/1.3 | 5.9/1.2 | 4.6/1.1 | 3.1/0.9 | 2.2/0.7 |
| 10 mm, 1 m/s | 133/6.2 | 77/5.3 | 51/4.5 | 27/3.5 | 17/2.8 | 12/2.3 | 8.6/2.0 | 5.2/1.5 | 3.5/1.1 |
| 10 mm, 2 m/s | 174/8.1 | 96/6.6 | 61/5.5 | 32/4.1 | 19/3.2 | 13/2.6 | 9.4/2.1 | 5.6/1.6 | 3.7/1.2 |

Scaling rules:
- The **100 K budget doubles every number.**
- **F_peak = 1.5** multiplies P by 3.2/1.5 = 2.13.
- **k_g = 15** lowers P by about 5–15%. The loss is largest where the film is thin, i.e. the turbulent cases at small g.

### There is no interior optimum: thinner is always better thermally

With g_c fixed, P(g) **decreases monotonically**. Dividing numerator and denominator by g gives

```
P ∝ (1 − f_pad) B / { [ g/(8 k_f) + (t_w/k_g + 1/h)/2 ] (g + g_c + 2 t_w) }
```

Both factors in the denominator grow with g. The thermal limit as g → 0 is finite:
`P0 = 2 B (1−f_pad) V_core / [F_peak (t_w/k_g + 1/h)(g_c + 2 t_w)]`.
In practice, g is set from below by **neutronics** (enough fuel volume for criticality), by manufacturing and pad
tolerances, and by fuel inventory and fill. So `results/design_points.csv` reports these, rather than a meaningless "optimum = thinnest":

- `P_limit_g0_MW`: the thermal limit at g → 0.
- `g_at_90pct_mm`, `g_at_50pct_mm`: the gap at which P falls to 90% and 50% of that limit.
- `g_ref_fraction_mm`: the gap that gives the reference critical-core fuel volume fraction 0.586/3.4 = **0.172**, with P, q''' and the dominant term there.

The table below uses the 50 K budget, F_peak = 3.2 and k_g = 30. For 100 K, double P. For F_peak = 1.5, multiply P by 2.13.

| g_c | v | g at φ_f = 0.172 | q''' there (MW/m³) | **P there (MW)** | dominant (share) | P0 as g → 0 (MW) | g at 50% P0 |
|---|---|---|---|---|---|---|---|
| 3 mm | 0.5–1 m/s | 2.97 mm | 21.9 | **4.0** | fuel 46% / film 43% | 9.1 | 2.4 mm |
| 3 mm | 2 m/s | 2.97 mm | 30.6 | **5.6** | fuel 64% | 19.3 | 1.3 mm |
| 5 mm | 0.5 m/s | 3.43 mm | 14.0 | **2.6** | film 53% | 5.2 | 3.4 mm |
| 5 mm | 1 m/s | 3.43 mm | 14.6 | **2.7** | film 51% | 5.6 | 3.2 mm |
| 5 mm | 2 m/s | 3.43 mm | 25.2 | **4.6** | fuel 71% | 19.3 | 1.2 mm |
| 10 mm | 0.5 m/s | 4.57 mm | 6.7 | **1.2** | film 62% | 2.2 | 5.3 mm |
| 10 mm | 1 m/s | 4.57 mm | 13.6 | **2.5** | fuel 68% | 9.6 | 1.8 mm |
| 10 mm | 2 m/s | 4.57 mm | 15.3 | **2.8** | fuel 76% | 14.5 | 1.3 mm |

At the reference fuel fraction, the best combination is **g_c = 3 mm with g ≈ 3 mm**. It gives about 4 MW for a bare core at a 50 K budget,
about 8.5 MW flattened, and about 17 MW for a flattened core at 100 K, even with laminar coolant. The coolant pressure drop at 0.5 m/s is
about 9 kPa across the core (about 19 kPa with pads ×2).

### Which resistance dominates

See `figures/dT_budget_breakdown.png`.

- **Graphite** is always a minor term: 5–15% of the budget at k_g = 30, rising to 10–26% at k_g = 15.
- **Laminar film** (h ≈ 500–1500 W/m²K): the film dominates for g below about 3–5 mm. The crossover g, where the film share equals the fuel share, is
  `g* = 4 k_f / h`. That is 2.8 mm for g_c = 3, 4.6 mm for g_c = 5, and 8.5 mm for g_c = 10. Above g*, fuel conduction dominates.
- **Turbulent film** (h ≈ 3000–6000 W/m²K): fuel conduction dominates from about 1.5 mm up.
- Because ΔT_fuel ∝ g² while the other terms go as g, **fuel conduction always wins at large g**. At g = 10 mm it takes 65–90% of the budget.

## Figures (phone-friendly, single column)

- `figures/q_allowed_and_power_vs_g_50K.png`: allowed peak q''' (log scale), core power (left axis F = 3.2, right axis F = 1.5), and fuel volume fraction vs g. One curve per g_c/v, with the 0.172 reference line.
- `figures/q_allowed_and_power_vs_g_100K.png`: same, for the 100 K budget.
- `figures/dT_budget_breakdown.png`: stacked ΔT budget (film / graphite / fuel) vs g for g_c = 5 mm, v = 1 m/s, at k_g = 30 and 15, labelled with q'''_allowed.

## Caveats

1. **No fuel convection** (conservative). Per `openfoam_fuel_slot/`, natural convection is negligible at 5 mm and
   roughly doubles the allowed q''' at 10 mm. So the large-g numbers here are pessimistic by up to about 2×, and the small-g numbers are accurate.
2. **Coolant velocity vs coolant temperature rise.** At the swept velocities, the core-average coolant bulk rise across a
   1.63 m flow path is only **0.3–3 K** (column `coolant_rise_across_core_K`). A design that uses the 550→650 °C (about 100 K) coolant
   rise implied by the budgets needs velocities of only **a few cm/s**. That is deep laminar flow, so use the laminar h
   (h = 8.235 k/(2 g_c)), which favours a narrow g_c, and ignore the turbulent rows. Otherwise the coolant ΔT across the core is
   small and the 50 K local budget applies almost everywhere.
3. **Fully developed h.** Entrance-region h is higher, so this is conservative near the inlet. Gnielinski is a tube
   correlation applied to parallel plates with D_h = 2 g_c, which is typically within ±20%. The 2300–3000 blend is arbitrary.
4. **Pads** are ignored in the heat path. Their pressure drop is only covered by the rough ×2 factor. Pads in the coolant gap
   will also trip and mix the flow, which could raise h (and dP) above these laminar values.
5. **Local budget.** B is applied at the peak-power location with the local coolant bulk temperature. P uses a single
   peaking factor, with no axial or radial coupling to the coolant temperature profile.
6. **No neutronics coupling.** The fuel fraction needed for criticality depends on g, g_c and the pitch (self-shielding and
   moderation). The reference 0.172 comes from the earlier slab-core result and is only a marker.
7. Radiation across the gaps, plate thermal stress, bowing, and gap-thickness tolerance are all ignored. At g = 1–2 mm, tolerance
   (for example ±0.2 mm) matters, because ΔT_fuel ∝ g².
8. Properties are fixed at 600 °C. FLiBe μ roughly doubles between 650 °C and 500 °C, which changes Re but not the laminar h.
