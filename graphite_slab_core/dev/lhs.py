import json, math, os, shutil
import numpy as np
import pandas as pd
from scipy.stats import qmc

PARAM_NAMES = ("slot_depth", "flat_width", "web_thickness", "wall_thickness")   # sampled geometry (cm)
# slot width is DERIVED: w = 2*slot_depth + flat_width (default round_location='both_sides')
OPTIONAL_PARAMS = ("enrichment", "core_radius")   # U-235 wt% of U, cylinder radius R (cm): sampled only if in bounds
ALL_PARAMS = PARAM_NAMES + OPTIONAL_PARAMS

def lhs_samples(bounds, n_samples, seed=None, fixed=None, geometry_opts=None, defaults=None, optimization=None):
    """Latin-hypercube samples.

    bounds   : {name: (low, high)} for every parameter to vary. The four geometry parameters must each be
               in `bounds` or `fixed`; 'enrichment' (U-235 wt%) and 'core_radius' (cm, H = 2R) are varied only
               if they appear in `bounds`, otherwise they take fixed[...] or defaults[...].
    fixed    : {name: value} parameters held constant (excluded from the hypercube)
    Returns a DataFrame with a 'feasible' flag and the reason for infeasible rows.
    """
    fixed = dict(fixed or {}); geometry_opts = dict(geometry_opts or {})
    defaults = {"enrichment": None, "core_radius": 70.0, **(defaults or {})}
    unknown = (set(bounds) | set(fixed)) - set(ALL_PARAMS)
    if unknown:
        raise ValueError(f"unknown parameter(s): {unknown}; allowed: {ALL_PARAMS}")
    free = [k for k in ALL_PARAMS if k in bounds and k not in fixed]
    missing = [k for k in PARAM_NAMES if k not in bounds and k not in fixed]
    if missing:
        raise ValueError(f"need bounds (or a fixed value) for {missing}")
    for k in free:
        lo, hi = bounds[k]
        if not ((0 <= lo < hi) if k == "flat_width" else (0 < lo < hi)):
            raise ValueError(f"bad bounds for {k}: {bounds[k]}")
    if "enrichment" in free and bounds["enrichment"][1] > 100:
        raise ValueError("enrichment is U-235 wt% of U: upper bound must be <= 100")
    if free:
        try:
            sampler = qmc.LatinHypercube(d=len(free), rng=seed, optimization=optimization)
        except TypeError:                                     # older SciPy
            sampler = qmc.LatinHypercube(d=len(free), seed=seed, optimization=optimization)
        unit = sampler.random(n_samples)
        vals = qmc.scale(unit, [bounds[k][0] for k in free], [bounds[k][1] for k in free])
    else:
        vals = np.empty((n_samples, 0))
    df = pd.DataFrame(vals, columns=free)
    for k, v in fixed.items():
        df[k] = float(v)
    for k in OPTIONAL_PARAMS:
        if k not in df:
            df[k] = defaults[k]
    df = df[list(ALL_PARAMS)]
    feas = [is_feasible(**{k: row[k] for k in PARAM_NAMES}, **geometry_opts) for row in df.to_dict("records")]
    df.insert(2, "slot_width", [resolve_params(**{k: r[k] for k in PARAM_NAMES}, **geometry_opts)["slot_width"]
                                if f else np.nan for r, (f, _) in zip(df.to_dict("records"), feas)])
    df["feasible"] = [f for f, _ in feas]
    df["reason"] = [r for _, r in feas]
    df.index.name = "sample"
    return df

def run_lhs(samples, work_dir, geometry_opts=None, run_opts=None, run_transport=False,
            threads=None, keep_going=True):
    """Build (and optionally run) one model per feasible sample in work_dir/sample_XXX."""
    geometry_opts = dict(geometry_opts or {}); run_opts = dict(run_opts or {})
    os.makedirs(work_dir, exist_ok=True)
    rows = []
    for idx, s in samples[samples["feasible"]].iterrows():
        geo = {k: float(s[k]) for k in PARAM_NAMES}
        enr = None if pd.isna(s.get("enrichment", np.nan)) else float(s["enrichment"])
        rad = float(s["core_radius"])
        d = os.path.join(work_dir, f"sample_{idx:03d}")
        os.makedirs(d, exist_ok=True)
        row = {"sample": idx, **geo, "enrichment": enr if enr is not None else MSRE_U235_WT_PCT,
               "core_radius": rad, **analytic_metrics(**geo, core_radius=rad, **geometry_opts), "keff": np.nan, "keff_std": np.nan,
               "runtime_s": np.nan, "status": "built"}
        try:
            model = build_model(**geo, enrichment=enr, core_radius=rad, **geometry_opts, **run_opts)
            model.export_to_model_xml(os.path.join(d, "model.xml"))
            with open(os.path.join(d, "params.json"), "w") as f:
                json.dump({**geo, "enrichment_u235_wt_pct": row["enrichment"], "core_radius": rad,
                       **geometry_opts, **run_opts}, f, indent=1)
            if run_transport:
                t0 = time.time()
                sp_path = model.run(cwd=d, output=False, threads=threads)
                row["runtime_s"] = time.time() - t0
                with openmc.StatePoint(sp_path) as sp:
                    row["keff"], row["keff_std"] = sp.keff.nominal_value, sp.keff.std_dev
                row["status"] = "ran"
        except Exception as e:
            row["status"] = f"error: {e}"
            if not keep_going:
                raise
        rows.append(row)
        print(f"sample {idx:3d}: " + ", ".join(f"{k}={geo[k]:.3f}" for k in PARAM_NAMES) + f", w={row['slot_width']:.3f}"
              + f", enr={row['enrichment']:.2f}, R={rad:.1f} | fuel_vf={row['fuel_vf']:.3f} C/F={row['graphite_to_fuel']:.2f}"
              + (f" | k={row['keff']:.5f}+/-{row['keff_std']:.5f} ({row['runtime_s']:.0f}s)" if row['status'] == 'ran' else f" | {row['status']}"))
    return pd.DataFrame(rows).set_index("sample")
