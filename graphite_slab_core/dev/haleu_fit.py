"""Critical-radius fit for the HALEU R scan + phone-friendly plot.
One-group bare-core form with H = 2R:  1/k_eff = a + b/(R + delta)^2   (a ~ 1/k_inf,eff, b ~ M^2 B^2 R^2 / k_inf).
Rc from 1/k = 1; uncertainty by parametric Monte Carlo on the k-eff statistical errors.
Writes results/haleu_critical.csv and figures/haleu_R_scan.png."""
import os, sys, json
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use("Agg")
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

if __name__ == "__main__":
    W = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
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
    print(crit.T.to_string())

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
        fig.savefig(os.path.join(W, "figures", "haleu_R_scan.png"), dpi=150)
