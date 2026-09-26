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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
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
        fig.savefig(os.path.join(W, "figures", "reflector_scan.png"), dpi=150)
    print(sc[["reflector_radial", "keff", "keff_std", "leakage_fraction", "dk_pcm", "drho_pcm", "marginal_pcm_per_cm",
              "R_equiv_bare", "equiv_radius_gain_cm"]].round(4).to_string())
    print(direct.round(3).to_string() if len(direct) else "no direct critical search yet")
