"""Analysis + phone-friendly plot for the fixed-width (1.0 cm) fuel-slot depth scan and the optional web scan.
Optimum = maximum unit-cell k-inf; located with a parabola through the best point and its neighbours (if interior).
Writes results/fuel_w1_summary.csv and figures/fuel_depth_w1_scan.png."""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = os.environ.get("MSR_WORK_DIR", "/workspace/msr_slab")
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
        fig.savefig(os.path.join(W, "figures", "fuel_depth_w1_scan.png"), dpi=150)
    print(fd[["fuel_depth", "fuel_flat", "kinf", "kinf_std", "graphite_to_fuel", "fuel_vf", "coolant_vf", "keff_bare", "keff_bare_std",
              "leakage_bare", "critical_bare_R85", "fuel_salt_volume_m3", "fuel_salt_mass_kg", "u235_mass_kg", "rel_conduction_power_limit"]].round(4).to_string())
    if wb is not None:
        print(wb[["web_thickness", "pitch_y", "pitch_z", "kinf", "kinf_std", "graphite_to_fuel", "fuel_vf", "coolant_vf"]].round(4).to_string())
    print(summ.round(4).to_string())
