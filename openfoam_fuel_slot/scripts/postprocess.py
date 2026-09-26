#!/usr/bin/env python3
"""Post-process all finished fuel-slot cases.

Outputs
  results/results.csv            one row per case (peak T, fuel-to-wall dT, max velocity, convergence, ...)
  results/allowed_power.csv      q'''_avg allowed for peak fuel T <= 700 C, and implied core power
  results/profiles/<case>.csv    axial profiles (z, T_fuel_max, T_centre, T_wall, T_coolant_face, Uz_centre)
  figures/peak_T_vs_power.png    peak fuel T vs q'''_avg, convection vs conduction-only, both slot depths
  figures/axial_profiles.png     axial temperature profiles, convection vs conduction (near-limit cases)
  figures/T_field_example.png    temperature / velocity field of one case
"""
import csv, glob, json, os, re, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from foamio import read_mesh_centres, read_internal, read_patch, all_times  # noqa: E402

T_LIMIT = 700.0
V_FUEL, F_RADIAL, F_AXIAL = 0.586, 2.3, 1.4    # m3 fuel salt (graphite_slab_core Case A), peaking estimates


def grid(C, vals):
    """Map a structured-mesh cell list to a 2D (nz, nx) array."""
    xs = np.unique(np.round(C[:, 0], 9)); zs = np.unique(np.round(C[:, 2], 9))
    ix = np.searchsorted(xs, np.round(C[:, 0], 9)); iz = np.searchsorted(zs, np.round(C[:, 2], 9))
    A = np.full((len(zs), len(xs)) + vals.shape[1:], np.nan)
    A[iz, ix] = vals
    return xs, zs, A


def log_history(case):
    """Fuel-region max T per iteration from the solver log (first 'Min/max T' line of each iteration)."""
    logs = glob.glob(f"{case}/log.chtMultiRegion*Foam")
    if not logs:
        return np.array([])
    tmax, it_lines = [], open(logs[0]).read().split("\nTime = ")[1:]
    for blk in it_lines:
        m = re.search(r"Min/max T:\s*([-\d.eE+]+)\s+([-\d.eE+]+)", blk)
        if m:
            tmax.append(float(m.group(2)))
    return np.array(tmax)


def analyse(case):
    p = json.load(open(f"{case}/caseParams.json"))
    times = all_times(case)
    if not times:
        return None
    t = times[-1]
    C, P = read_mesh_centres(case, "fuel")
    T = read_internal(f"{case}/{t}/fuel/T", len(C))
    U = read_internal(f"{case}/{t}/fuel/U", len(C))
    Tw = read_patch(f"{case}/{t}/fuel/T", "fuel_to_graphite")
    Pw = P["fuel_to_graphite"]
    Cg, Pg = read_mesh_centres(case, "graphite")
    Tg = read_internal(f"{case}/{t}/graphite/T", len(Cg))
    xs, zs, TT = grid(C, T)
    _, _, UU = grid(C, U)
    # wall temperature per row: left and right interface
    zr = np.round(Pw[:, 2], 9)
    left = Pw[:, 0] < xs.mean()
    TwL = np.array([Tw[(zr == z) & left].mean() for z in zs])
    TwR = np.array([Tw[(zr == z) & ~left].mean() for z in zs])
    Tmax_row = np.nanmax(TT, axis=1)
    Tcen = TT[:, len(xs) // 2]
    dT_row = Tmax_row - np.maximum(TwL, TwR)
    # energy balance: heat leaving through coolant faces vs source
    xg = np.unique(np.round(Cg[:, 0], 9))
    Qout = 0.0
    for patch, xc in (("coolantL", xg[0]), ("coolantR", xg[-1])):
        Tf = read_patch(f"{case}/{t}/graphite/T", patch)
        zf = np.round(Pg[patch][:, 2], 9)
        sel = np.isclose(Cg[:, 0], xc)
        Tc_cell = dict(zip(np.round(Cg[sel, 2], 9), Tg[sel]))
        dn = abs(xc - (0.0 if patch == "coolantL" else 2 * p["w_graphite"] + p["d_f"]))
        dz = p["H"] / p["nz"]
        Tf = np.broadcast_to(np.atleast_1d(Tf), zf.shape)     # 'uniform' patch values are written as a scalar
        Qout += sum(p["k_graphite"] * (Tc_cell[z] - T_f) / dn * dz for z, T_f in zip(zf, Tf))
    Qin = p["q_avg"] * p["d_f"] * p["H"]           # W per metre of slot width (y)
    hist = log_history(case)
    n = len(hist)
    last = hist[-min(1000, n // 2):] if n else np.array([np.nan])
    last2k = hist[-min(2000, n // 2):] - 273.15 if n else np.array([np.nan])
    tprev = times[-2] if len(times) > 1 else None
    dU = np.nan
    if tprev:
        Up = read_internal(f"{case}/{tprev}/fuel/U", len(C))
        um = np.linalg.norm(U, axis=1).max()
        dU = abs(um - np.linalg.norm(Up, axis=1).max()) / max(um, 1e-12)
    ipk = int(np.argmax(T))
    umag = np.linalg.norm(U, axis=1)
    name = os.path.basename(case)
    Tcool = (873.15 if p["coolant"] == "uniform" else 823.15 + 100 * zs / p["H"])
    prof = dict(z_m=zs, T_fuel_max_C=Tmax_row - 273.15, T_fuel_centre_C=Tcen - 273.15,
                T_wall_left_C=TwL - 273.15, T_wall_right_C=TwR - 273.15,
                T_coolant_C=np.broadcast_to(Tcool, zs.shape) - 273.15, Uz_centre_mm_s=UU[:, len(xs) // 2, 2] * 1e3)
    os.makedirs(f"{ROOT}/results/profiles", exist_ok=True)
    with open(f"{ROOT}/results/profiles/{name}.csv", "w", newline="") as f:
        wr = csv.writer(f); wr.writerow(prof.keys())
        for row in zip(*prof.values()):
            wr.writerow([f"{v:.5g}" for v in row])
    mesh = "coarse" if "meshCoarse" in name else ("fine" if "meshFine" in name else "base")
    return dict(case=name, d_f_cm=p["d_f"] * 100, q_avg_MW_m3=p["q_avg"] / 1e6, q_peak_MW_m3=p["q_peak"] / 1e6,
                mode="convection" if p["gravity"] else "conduction", coolant=p["coolant"],
                h_W_m2K=p["h"] if p["h"] > 0 else "fixedT", mesh=mesh, cells_fuel=len(C), iterations=n,
                T_peak_C=T.max() - 273.15, z_peak_m=C[ipk, 2],
                dT_fuel_wall_max_K=dT_row.max(), z_dTmax_m=zs[np.argmax(dT_row)],
                T_wall_max_C=max(TwL.max(), TwR.max()) - 273.15, T_graphite_max_C=Tg.max() - 273.15,
                wall_LR_asym_max_K=np.abs(TwL - TwR).max(),
                U_max_mm_s=umag.max() * 1e3, Uz_up_max_mm_s=U[:, 2].max() * 1e3, Uz_down_max_mm_s=-U[:, 2].min() * 1e3,
                Tpeak_iter_range_K=float(np.ptp(last)), Umax_rel_change_last1000=dU,
                Tpeak_iter_mean_C=float(last2k.mean()), Tpeak_iter_min_C=float(last2k.min()),
                Tpeak_iter_max_C=float(last2k.max()),
                energy_balance_out_over_in=Qout / Qin,
                steady_converged=bool(abs(Qout / Qin - 1) < 0.005 and (not n or float(np.ptp(last)) < 0.5)),
                _prof=prof, _fields=(xs, zs, TT, UU, Cg, Tg))


def interp_limit(q, T, lim=T_LIMIT):
    q, T = np.asarray(q, float), np.asarray(T, float)
    o = np.argsort(q); q, T = q[o], T[o]
    if T.max() < lim or T.min() > lim:
        return np.nan
    i = np.searchsorted(T, lim)       # T increases with q
    # interpolate linearly in (q, T)
    return q[i - 1] + (lim - T[i - 1]) * (q[i] - q[i - 1]) / (T[i] - T[i - 1])


def main():
    rows = []
    for case in sorted(glob.glob(f"{ROOT}/cases/*")):
        if not os.path.isfile(f"{case}/caseParams.json"):
            continue
        try:
            r = analyse(case)
        except Exception as e:   # unfinished case
            print("skip", case, e)
            continue
        if r:
            rows.append(r)
    # transient restarts: use the time-mean peak temperature as the best estimate for that case
    for r in rows:
        r["T_peak_best_C"] = r["T_peak_C"]
        r["T_best_lo_C"] = r["T_best_hi_C"] = r["T_peak_C"]
        if not r["steady_converged"] and r["mode"] == "convection":
            # wandering steady iterates: best estimate = mean of the fuel peak T over the last 2000 iterations
            r["T_peak_best_C"] = r["Tpeak_iter_mean_C"]
            r["T_best_lo_C"], r["T_best_hi_C"] = r["Tpeak_iter_min_C"], r["Tpeak_iter_max_C"]
        ts = transient_stats(f"{ROOT}/cases/{r['case']}_transient")
        if ts:
            r.update(T_peak_transient_mean_C=ts["T_mean_2nd_half"], T_peak_transient_min_C=ts["T_min_2nd_half"],
                     T_peak_transient_max_C=ts["T_max_2nd_half"], transient_t_end_s=ts["t_end"])
            r["T_peak_best_C"] = ts["T_mean_2nd_half"]
            r["T_best_lo_C"], r["T_best_hi_C"] = ts["T_min_2nd_half"], ts["T_max_2nd_half"]
    # convection enhancement per case: conduction/convection ratio of max fuel-to-wall dT and of peak rise
    idx = {(r["d_f_cm"], r["q_avg_MW_m3"], r["mode"], r["coolant"], r["h_W_m2K"], r["mesh"]): r for r in rows}
    for r in rows:
        c = idx.get((r["d_f_cm"], r["q_avg_MW_m3"], "conduction", r["coolant"], r["h_W_m2K"], "base"))
        if r["mode"] == "convection" and c:
            r["dT_ratio_cond_over_conv"] = c["dT_fuel_wall_max_K"] / r["dT_fuel_wall_max_K"]
            r["T_peak_cond_C"] = c["T_peak_C"]
            r["T_peak_reduction_K"] = c["T_peak_C"] - r["T_peak_best_C"]
    keys = [k for k in rows[0] if not k.startswith("_")] + ["T_peak_transient_mean_C", "T_peak_transient_min_C",
            "T_peak_transient_max_C", "transient_t_end_s", "T_peak_best_C", "T_best_lo_C", "T_best_hi_C", "dT_ratio_cond_over_conv", "T_peak_cond_C", "T_peak_reduction_K"]
    keys = list(dict.fromkeys(keys))
    with open(f"{ROOT}/results/results.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, keys, extrasaction="ignore"); wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items() if not k.startswith("_")})
    # allowed power
    lim_rows = []
    for d in sorted({r["d_f_cm"] for r in rows}):
        for cool in ("linear", "uniform"):
            out = {}
            for mode in ("convection", "conduction"):
                sel = [r for r in rows if r["d_f_cm"] == d and r["mode"] == mode and r["coolant"] == cool
                       and r["h_W_m2K"] == "fixedT" and r["mesh"] == "base"]
                if len(sel) >= 2:
                    out[mode] = interp_limit([r["q_avg_MW_m3"] for r in sel], [r["T_peak_best_C"] for r in sel])
                    out[mode + "_note"] = "interpolated" if not np.isnan(out[mode]) else ""
                    if np.isnan(out[mode]):     # limit outside the simulated range: extrapolate the last two points
                        sel = sorted(sel, key=lambda r: r["q_avg_MW_m3"])[-2:]
                        (q0, T0), (q1, T1) = [(r["q_avg_MW_m3"], r["T_peak_best_C"]) for r in sel]
                        out[mode] = q1 + (T_LIMIT - T1) * (q1 - q0) / (T1 - T0)
                        out[mode + "_note"] = "linear extrapolation"
                    if any(not r["steady_converged"] and "T_peak_transient_mean_C" not in r for r in sel):
                        out[mode + "_note"] += "; uses non-converged steady points"
            if out:
                qa = out.get("convection", np.nan); qc = out.get("conduction", np.nan)
                selc = [r for r in rows if r["d_f_cm"] == d and r["mode"] == "convection" and r["coolant"] == cool
                        and r["h_W_m2K"] == "fixedT" and r["mesh"] == "base"]
                qlo = interp_limit([r["q_avg_MW_m3"] for r in selc], [r["T_best_hi_C"] for r in selc])
                qhi = interp_limit([r["q_avg_MW_m3"] for r in selc], [r["T_best_lo_C"] for r in selc])
                lim_rows.append(dict(d_f_cm=d, coolant=cool, q_avg_allowed_conv_MW_m3=qa, q_avg_allowed_cond_MW_m3=qc,
                                     q_conv_range_lo=qlo, q_conv_range_hi=qhi, note_conv=out.get("convection_note", ""), note_cond=out.get("conduction_note", ""),
                                     power_factor_conv_over_cond=qa / qc,
                                     P_core_conv_MW=qa * V_FUEL / F_RADIAL, P_core_cond_MW=qc * V_FUEL / F_RADIAL))
    with open(f"{ROOT}/results/allowed_power.csv", "w", newline="") as f:
        wr = csv.DictWriter(f, lim_rows[0].keys()); wr.writeheader()
        for r in lim_rows:
            wr.writerow({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in r.items()})
    for r in lim_rows:
        print(r)
    plots(rows, lim_rows)


def plots(rows, lim_rows):
    plt.rcParams.update({"font.size": 15, "axes.titlesize": 16, "legend.fontsize": 12})
    os.makedirs(f"{ROOT}/figures", exist_ok=True)
    # ---------------- peak T vs power
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 11.5))
    for ax, cool, title in zip(axs, ("linear", "uniform"), ("Coolant 550→650 °C (bottom→top)", "Coolant 600 °C uniform")):
        for d, col in ((0.5, "tab:blue"), (1.0, "tab:red")):
            for mode, ls, mk in (("convection", "-", "o"), ("conduction", "--", "s")):
                sel = sorted([r for r in rows if abs(r["d_f_cm"] - d) < 1e-9 and r["mode"] == mode and r["coolant"] == cool
                              and r["h_W_m2K"] == "fixedT" and r["mesh"] == "base"], key=lambda r: r["q_avg_MW_m3"])
                if sel:
                    ax.plot([r["q_avg_MW_m3"] for r in sel], [r["T_peak_best_C"] for r in sel], ls, marker=mk, color=col,
                            lw=2.2, ms=7, label=f"d = {d:g} cm, {'natural conv.' if mode == 'convection' else 'conduction only'}")
                    nc = [r for r in sel if not r["steady_converged"]]
                    if nc:
                        ax.plot([r["q_avg_MW_m3"] for r in nc], [r["T_peak_best_C"] for r in nc], "o", mfc="white",
                                mec=col, ms=9, mew=2, label=None)
                    for r in sel:
                        if r["T_best_hi_C"] - r["T_best_lo_C"] > 0.5:
                            ax.errorbar(r["q_avg_MW_m3"], r["T_peak_best_C"], yerr=[[r["T_peak_best_C"] - r["T_best_lo_C"]],
                                        [r["T_best_hi_C"] - r["T_peak_best_C"]]], color=col, capsize=5, lw=1.5)
        ax.axhline(T_LIMIT, color="k", lw=1.2, ls=":")
        ax.text(0.98, T_LIMIT - 3, "700 °C", ha="right", va="top", transform=ax.get_yaxis_transform())
        ax.text(0.98, 0.03, "open marker = unsteady flow (steady solver not converged):\nmean of iterates or transient; bar = min/max", ha="right",
                va="bottom", fontsize=10, transform=ax.transAxes)
        ax.set_xscale("log"); ax.grid(True, which="both", alpha=0.3)
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}" if v in (2, 3, 5, 20, 30) else ""))
        ax.set_xlabel("average fuel q''' (MW/m³)"); ax.set_ylabel("peak fuel T (°C)"); ax.set_title(title)
        ax.legend(loc="upper left")
    fig.suptitle("2D OpenFOAM fuel slot, 1.63 m tall, cooled both faces\nthrough 0.5 cm graphite (axial peak/avg 1.4)",
                 fontsize=14)
    fig.tight_layout()
    fig.savefig(f"{ROOT}/figures/peak_T_vs_power.png", dpi=130)
    plt.close(fig)

    # ---------------- axial profiles (near-limit cases)
    by = {r["case"]: r for r in rows}
    pairs = [("d05mm_q20_conv", "d05mm_q20_cond", "d = 0.5 cm, q''' = 20 MW/m³"),
             ("d10mm_q5_conv", "d10mm_q5_cond", "d = 1.0 cm, q''' = 5 MW/m³ (unsteady: snapshot)")]
    fig, axs = plt.subplots(1, 2, figsize=(9, 10), sharey=True)
    for ax, (cv, cd, title) in zip(axs, pairs):
        if cv not in by or cd not in by:
            continue
        a, b = by[cv]["_prof"], by[cd]["_prof"]
        ax.plot(a["T_fuel_max_C"], a["z_m"], "-", color="tab:red", lw=2.2, label="fuel max, natural conv.")
        ax.plot(b["T_fuel_max_C"], b["z_m"], "--", color="tab:red", lw=2, label="fuel max, conduction")
        ax.plot(np.maximum(a["T_wall_left_C"], a["T_wall_right_C"]), a["z_m"], "-", color="tab:gray", lw=1.8, label="fuel/graphite wall, conv.")
        ax.plot(np.maximum(b["T_wall_left_C"], b["T_wall_right_C"]), b["z_m"], "--", color="tab:gray", lw=1.5, label="wall, conduction")
        ax.plot(a["T_coolant_C"], a["z_m"], ":", color="tab:blue", lw=2, label="coolant-side T")
        ax.axvline(T_LIMIT, color="k", lw=1, ls=":")
        ax.set_title(title, fontsize=13); ax.set_xlabel("T (°C)"); ax.grid(alpha=0.3)
    axs[0].set_ylabel("height z (m)")
    axs[0].legend(loc="lower right", fontsize=10)
    fig.suptitle("Axial temperature profiles\n(convection moves heat to the top of the slot)", fontsize=14)
    fig.tight_layout()
    fig.savefig(f"{ROOT}/figures/axial_profiles.png", dpi=130)
    plt.close(fig)

    # ---------------- field example
    ex = by.get("d05mm_q20_conv")
    if ex:
        xs, zs, TT, UU, Cg, Tg = ex["_fields"]
        gx, gz, GT = grid(Cg, Tg)
        fig = plt.figure(figsize=(8, 11))
        gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], width_ratios=[1, 1, 1])
        axT = fig.add_subplot(gs[:, 0]); axU = fig.add_subplot(gs[:, 1])
        axTop = fig.add_subplot(gs[0, 2]); axMid = fig.add_subplot(gs[1, 2])
        vmin, vmax = np.nanmin(GT) - 273.15, np.nanmax(TT) - 273.15
        for X, Z, A in ((gx, gz, GT), (xs, zs, TT)):
            A2 = A.copy()
            m = axT.pcolormesh(X * 1e3, Z, A2 - 273.15, shading="nearest", cmap="inferno", vmin=vmin, vmax=vmax)
        axT.axvline(5, color="w", lw=0.6); axT.axvline(5 + ex["d_f_cm"] * 10, color="w", lw=0.6)
        axT.set_xlabel("x (mm)"); axT.set_ylabel("z (m)"); axT.set_title("T (°C)")
        fig.colorbar(m, ax=axT, orientation="horizontal", pad=0.06)
        uz = UU[:, :, 2] * 1e3
        lim = np.nanmax(np.abs(uz))
        mu = axU.pcolormesh(xs * 1e3, zs, uz, shading="nearest", cmap="RdBu_r", vmin=-lim, vmax=lim)
        axU.set_xlabel("x (mm)"); axU.set_title("Uz (mm/s)"); axU.set_yticklabels([])
        fig.colorbar(mu, ax=axU, orientation="horizontal", pad=0.06)
        for ax, (z0, z1), ttl in ((axTop, (zs[-1] - 0.06, zs[-1] + 0.001), "top 6 cm"), (axMid, (0.80, 0.86), "z = 0.80–0.86 m")):
            sel = (zs >= z0) & (zs <= z1)
            ax.pcolormesh(xs * 1e3, zs[sel], TT[sel] - 273.15, shading="nearest", cmap="inferno", vmin=vmin, vmax=vmax)
            ax.pcolormesh(gx * 1e3, gz[(gz >= z0) & (gz <= z1)], GT[(gz >= z0) & (gz <= z1)] - 273.15, shading="nearest",
                          cmap="inferno", vmin=vmin, vmax=vmax)
            X, Z = np.meshgrid(xs * 1e3, zs[sel])
            st = max(1, len(xs) // 8)
            ax.quiver(X[::3, ::st], Z[::3, ::st], UU[sel][::3, ::st, 0], UU[sel][::3, ::st, 2], color="c", scale=None)
            ax.set_title(ttl, fontsize=12); ax.set_xlabel("x (mm)")
        fig.suptitle(f"d = {ex['d_f_cm']:g} cm, q'''avg = {ex['q_avg_MW_m3']:g} MW/m³, natural convection\n"
                     f"peak {ex['T_peak_C']:.0f} °C at z = {ex['z_peak_m']:.2f} m (x stretched)", fontsize=13)
        fig.tight_layout()
        fig.savefig(f"{ROOT}/figures/T_field_example.png", dpi=130)
        plt.close(fig)


def transient_stats(case):
    log = f"{case}/log.chtMultiRegionFoam"
    if not os.path.isfile(log):
        return None
    t, T = [], []
    for blk in open(log).read().split("\nTime = ")[1:]:
        m = re.search(r"Min/max T:\s*([-\d.eE+]+)\s+([-\d.eE+]+)", blk)
        if m:
            t.append(float(blk.split()[0])); T.append(float(m.group(2)) - 273.15)
    t, T = np.array(t), np.array(T)
    half = t > t[-1] / 2
    return dict(t=t, T=T, t_end=t[-1], T_mean_2nd_half=T[half].mean(), T_min_2nd_half=T[half].min(),
                T_max_2nd_half=T[half].max())


def transient_plot(case=f"{ROOT}/cases/d10mm_q5_conv_transient"):
    """Peak fuel T vs physical time for the transient restart (steadiness check)."""
    log = f"{case}/log.chtMultiRegionFoam"
    if not os.path.isfile(log):
        return None
    t, T = [], []
    for blk in open(log).read().split("\nTime = ")[1:]:
        m = re.search(r"Min/max T:\s*([-\d.eE+]+)\s+([-\d.eE+]+)", blk)
        if m:
            t.append(float(blk.split()[0])); T.append(float(m.group(2)) - 273.15)
    t, T = np.array(t), np.array(T)
    ref = [r for r in csv.DictReader(open(f"{ROOT}/results/results.csv")) if r["case"] in ("d10mm_q5_cond", "d10mm_q5_conv")]
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.plot(t, T, lw=1.5, color="tab:red", label="transient, natural convection")
    for r in ref:
        ax.axhline(float(r["T_peak_C"]), ls="--" if "cond" in r["case"] else ":", color="k",
                   label=("steady conduction only" if "cond" in r["case"] else "steady solver (last iterate)"))
    ax.set_xlabel("time (s)"); ax.set_ylabel("peak fuel T (°C)"); ax.grid(alpha=0.3); ax.legend(fontsize=11)
    ax.set_title("d = 1.0 cm, q" + "'" * 3 + " = 5 MW/m³: transient check", fontsize=14)
    fig.tight_layout(); fig.savefig(f"{ROOT}/figures/transient_check.png", dpi=130); plt.close(fig)
    half = t > t[-1] / 2
    return dict(t_end=t[-1], T_mean_2nd_half=T[half].mean(), T_min_2nd_half=T[half].min(), T_max_2nd_half=T[half].max())


if __name__ == "__main__":
    main()
    print("transient:", transient_plot())
