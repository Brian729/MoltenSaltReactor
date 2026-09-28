#!/usr/bin/env python3
"""1D steady heat-path model for a pad-supported flat-plate MSR core.

Repeat unit (normal to the plates):
    fuel gap g | graphite plate t_w | coolant gap g_c | graphite plate t_w | (next fuel gap) ...
Each fuel gap is cooled through both faces; each coolant gap is heated from both faces.

Series temperature rise from coolant bulk to fuel centreline (local, away from pads):
    q''        = q''' g / 2                       heat flux per fuel-gap face
    dT_fuel    = q''' g^2 / (8 k_f)               conduction in a slab heated uniformly, both faces cooled
    dT_graph   = q'' t_w / k_g                    conduction through one plate
    dT_film    = q'' / h_c                        coolant film (both coolant walls heated)
    dT_total   = q''' [ g^2/(8 k_f) + (g/2)(t_w/k_g + 1/h_c) ]  <= budget
=> q'''_allowed(g) = budget / [ g^2/(8 k_f) + (g/2)(t_w/k_g + 1/h_c) ]

Core power: P = q'''_allowed / F_peak * V_core * phi_fuel,  phi_fuel = g (1 - f_pad) / (g + g_c + 2 t_w)

Run: python3 heat_path_1d.py   -> results/*.csv, figures/*.png
"""
import csv, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- inputs
G_LIST_MM = [1, 1.5, 2, 3, 4, 5, 6, 8, 10]      # fuel gap thickness sweep
GC_LIST_MM = [3, 5, 10]                          # coolant gap thickness
V_LIST = [0.5, 1.0, 2.0]                         # coolant velocity, m/s
KG_LIST = [30.0, 15.0]                           # graphite k: baseline / irradiated, W/m/K
BUDGETS = [50.0, 100.0]                          # K, coolant bulk -> fuel centreline
PEAKING = [3.2, 1.5]                             # total peaking: bare (2.3 x 1.4) / flattened
T_W = 0.005                                      # graphite plate thickness, m
K_F = 1.05                                       # fuel salt k, W/m/K (openfoam_fuel_slot README)
V_CORE = 3.4                                     # m3, R = 81.5 cm, H = 2R (graphite_slab_core)
R_CORE = 0.815
PAD_D, PAD_P = 0.010, 0.035                      # pad diameter and hex pitch, m
F_PAD = (math.pi / 4 * PAD_D**2) / (math.sqrt(3) / 2 * PAD_P**2)   # = 0.074 area fraction
FLOW_LENGTH = 2 * R_CORE                         # coolant path across the core, m (for dP and coolant rise)

# FLiBe (2LiF-BeF2) at T_COOL, ORNL/TM-2006/12 (Williams) and Romatoski & Hu, Ann. Nucl. Energy 109 (2017)
T_COOL = 600.0 + 273.15
RHO_C = 2413.0 - 0.488 * T_COOL                  # kg/m3
MU_C = 1.16e-4 * math.exp(3755.0 / T_COOL)       # Pa s  (0.116 exp(3755/T) cP)
K_C = 1.1                                        # W/m/K
CP_C = 2386.0                                    # J/kg/K
PR_C = MU_C * CP_C / K_C

NU_LAM = 8.235     # fully developed laminar, parallel plates, both walls uniform heat flux (Shah & London)
RE_LAM, RE_TURB = 2300.0, 3000.0


def coolant(gc, v):
    """Re, Nu, h, regime, Darcy f, dP/L for a parallel-plate channel of gap gc (Dh = 2 gc)."""
    dh = 2 * gc
    re = RHO_C * v * dh / MU_C
    def gniel(r):
        f = (0.790 * math.log(r) - 1.64) ** -2          # Petukhov smooth-tube friction factor
        nu = (f / 8) * (r - 1000) * PR_C / (1 + 12.7 * math.sqrt(f / 8) * (PR_C ** (2 / 3) - 1))
        return nu, f
    if re < RE_LAM:
        nu, f, reg = NU_LAM, 96.0 / re, "laminar"
    elif re >= RE_TURB:
        nu, f = gniel(re); reg = "turbulent (Gnielinski)"
    else:   # transition: linear blend between laminar at 2300 and Gnielinski at 3000
        w = (re - RE_LAM) / (RE_TURB - RE_LAM)
        nt, ft = gniel(RE_TURB)
        nu = (1 - w) * NU_LAM + w * nt; f = (1 - w) * 96.0 / RE_LAM + w * ft; reg = "transitional (blend)"
    h = nu * K_C / dh
    dpdl = f / dh * 0.5 * RHO_C * v ** 2
    return dict(Re=re, Nu=nu, h=h, regime=reg, f_darcy=f, dP_per_m_Pa=dpdl)


def resistances(g, gc, v, kg):
    """Temperature rise per unit q''' (K per W/m3) for each term."""
    c = coolant(gc, v)
    return dict(fuel=g**2 / (8 * K_F), graphite=(g / 2) * T_W / kg, film=(g / 2) / c["h"]), c


def run():
    rows = []
    for gc_mm in GC_LIST_MM:
        gc = gc_mm / 1000
        for v in V_LIST:
            for kg in KG_LIST:
                for g_mm in G_LIST_MM:
                    g = g_mm / 1000
                    r, c = resistances(g, gc, v, kg)
                    rt = sum(r.values())
                    phi_f = g * (1 - F_PAD) / (g + gc + 2 * T_W)
                    phi_c = gc * (1 - F_PAD) / (g + gc + 2 * T_W)
                    dom = max(r, key=r.get)
                    for bud in BUDGETS:
                        q = bud / rt
                        for fp in PEAKING:
                            P = q / fp * V_CORE * phi_f
                            # coolant bulk rise across the core if all coolant gaps flow across the core diameter
                            a_flow = phi_c * (2 * R_CORE) * (2 * R_CORE)   # gap area normal to flow, m2
                            mdot = RHO_C * v * a_flow
                            rows.append(dict(
                                g_mm=g_mm, gc_mm=gc_mm, v_m_s=v, k_graphite=kg, budget_K=bud, peaking=fp,
                                q_allowed_MW_m3=q / 1e6, fuel_fraction=phi_f, V_fuel_m3=V_CORE * phi_f, P_core_MW=P / 1e6,
                                heat_flux_face_kW_m2=q * g / 2 / 1e3,
                                dT_fuel_K=q * r["fuel"], dT_graphite_K=q * r["graphite"], dT_film_K=q * r["film"],
                                frac_fuel=r["fuel"] / rt, frac_graphite=r["graphite"] / rt, frac_film=r["film"] / rt,
                                dominant=dom, Re=c["Re"], Nu=c["Nu"], h_W_m2K=c["h"], regime=c["regime"],
                                dP_per_m_kPa=c["dP_per_m_Pa"] / 1e3, dP_per_m_pads_x2_kPa=2 * c["dP_per_m_Pa"] / 1e3,
                                dP_core_kPa=c["dP_per_m_Pa"] * FLOW_LENGTH / 1e3,
                                coolant_rise_across_core_K=P / (mdot * CP_C)))
    return rows


PHI_REF = 0.586 / 3.4    # fuel volume fraction of the HALEU critical reference core (graphite_slab_core, Case A)


def design_points(gc, v, kg, bud, fp):
    """P(g) falls monotonically with g for fixed g_c (see README), so there is no interior optimum.
    Report the g -> 0 limit, the g where P drops to 90 % / 50 % of it, and the g that gives the reference
    fuel fraction PHI_REF (a stand-in for the neutronic constraint), with the power there."""
    gcm, S = gc / 1000, gc / 1000 + 2 * T_W
    def P_of(g):
        r, _ = resistances(g, gcm, v, kg)
        return bud / sum(r.values()) / fp * V_CORE * g * (1 - F_PAD) / (g + S), r
    c = coolant(gcm, v)
    R_lin = T_W / kg + 1 / c["h"]
    P0 = 2 * bud * (1 - F_PAD) / (R_lin * S) / fp * V_CORE          # analytic g -> 0 limit
    gs = np.linspace(1e-5, 20e-3, 20001)
    Ps = np.array([P_of(g)[0] for g in gs])
    g90 = gs[np.argmax(Ps < 0.9 * P0)]; g50 = gs[np.argmax(Ps < 0.5 * P0)]
    g_ref = PHI_REF * S / ((1 - F_PAD) - PHI_REF)
    P_ref, r = P_of(g_ref)
    rt = sum(r.values())
    return dict(gc_mm=gc, v_m_s=v, k_graphite=kg, budget_K=bud, peaking=fp, h_W_m2K=c["h"],
                P_limit_g0_MW=P0 / 1e6, g_at_90pct_mm=g90 * 1e3, g_at_50pct_mm=g50 * 1e3,
                g_ref_fraction_mm=g_ref * 1e3, P_at_g_ref_MW=P_ref / 1e6, q_at_g_ref_MW_m3=bud / rt / 1e6,
                frac_fuel=r["fuel"] / rt, frac_graphite=r["graphite"] / rt, frac_film=r["film"] / rt,
                dominant_at_g_ref=max(r, key=r.get))


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, rows[0].keys()); w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})


def plots(rows):
    plt.rcParams.update({"font.size": 14, "legend.fontsize": 10.5})
    cols = {3: "tab:blue", 5: "tab:green", 10: "tab:red"}
    lss = {0.5: ":", 1.0: "--", 2.0: "-"}
    for bud in BUDGETS:
        fig, axs = plt.subplots(3, 1, figsize=(7.5, 13.5), sharex=True, gridspec_kw=dict(height_ratios=[1.1, 1.1, 0.6]))
        for gc in GC_LIST_MM:
            for v in V_LIST:
                sel = [r for r in rows if r["gc_mm"] == gc and r["v_m_s"] == v and r["k_graphite"] == 30
                       and r["budget_K"] == bud and r["peaking"] == 3.2]
                g = [r["g_mm"] for r in sel]
                lab = f"g_c = {gc} mm, v = {v:g} m/s (Re {sel[0]['Re']:.0f})"
                axs[0].plot(g, [r["q_allowed_MW_m3"] for r in sel], lss[v], color=cols[gc], lw=2, marker="o", ms=4, label=lab)
                axs[1].plot(g, [r["P_core_MW"] for r in sel], lss[v], color=cols[gc], lw=2, marker="o", ms=4)
        axs[0].set_yscale("log"); axs[0].set_ylabel("allowed peak q''' (MW/m³)")
        axs[0].set_title(f"Budget {bud:.0f} K coolant bulk → fuel centreline\nk_graphite = 30 W/m·K, t_w = 5 mm", fontsize=13)
        axs[0].legend(loc="upper right", ncol=1); axs[0].grid(True, which="both", alpha=0.3)
        axs[1].set_ylabel("core power (MW), peaking 3.2")
        for gc in GC_LIST_MM:
            sel = [r for r in rows if r["gc_mm"] == gc and r["v_m_s"] == 1.0 and r["k_graphite"] == 30
                   and r["budget_K"] == bud and r["peaking"] == 3.2]
            axs[2].plot([r["g_mm"] for r in sel], [r["fuel_fraction"] for r in sel], color=cols[gc], lw=2, label=f"g_c = {gc} mm")
        axs[2].axhline(PHI_REF, color="k", ls=":", lw=1.2)
        axs[2].text(10, PHI_REF + 0.01, "reference critical core 0.172", ha="right", va="bottom", fontsize=10)
        axs[2].set_ylabel("fuel volume fraction"); axs[2].set_xlabel("fuel gap g (mm)"); axs[2].grid(alpha=0.3)
        axs[2].legend(loc="lower right", fontsize=10)
        axs[1].grid(alpha=0.3)
        ax2 = axs[1].twinx(); lo, hi = axs[1].get_ylim(); ax2.set_ylim(lo * 3.2 / 1.5, hi * 3.2 / 1.5)
        ax2.set_ylabel("core power (MW), peaking 1.5")
        axs[1].set_title("V_core = 3.4 m³; P = q'''/peaking × V_core × fuel fraction", fontsize=12)
        fig.tight_layout(); fig.savefig(f"{HERE}/figures/q_allowed_and_power_vs_g_{bud:.0f}K.png", dpi=130); plt.close(fig)

    # stacked bar breakdown, representative case
    gc, v = 5, 1.0
    fig, axs = plt.subplots(2, 1, figsize=(7.5, 11))
    for ax, kg in zip(axs, KG_LIST):
        sel = [r for r in rows if r["gc_mm"] == gc and r["v_m_s"] == v and r["k_graphite"] == kg
               and r["budget_K"] == 50 and r["peaking"] == 3.2]
        x = np.arange(len(sel))
        b1 = np.array([r["dT_film_K"] for r in sel]); b2 = np.array([r["dT_graphite_K"] for r in sel])
        b3 = np.array([r["dT_fuel_K"] for r in sel])
        ax.bar(x, b1, color="tab:blue", label="coolant film")
        ax.bar(x, b2, bottom=b1, color="tab:gray", label="graphite plate")
        ax.bar(x, b3, bottom=b1 + b2, color="tab:red", label="fuel conduction")
        for xi, r in zip(x, sel):
            ax.text(xi, 51, f"{r['q_allowed_MW_m3']:.0f}", ha="center", va="bottom", fontsize=10)
        ax.set_xticks(x, [f"{r['g_mm']:g}" for r in sel]); ax.set_ylim(0, 60)
        ax.set_xlabel("fuel gap g (mm)"); ax.set_ylabel("ΔT share of 50 K budget (K)")
        ax.set_title(f"g_c = {gc} mm, v = {v:g} m/s (h = {sel[0]['h_W_m2K']:.0f} W/m²K), k_g = {kg:g} W/m·K\n"
                     "numbers = allowed peak q''' (MW/m³)", fontsize=12)
        ax.legend(loc="lower left", fontsize=10.5, framealpha=0.9)
    fig.tight_layout(); fig.savefig(f"{HERE}/figures/dT_budget_breakdown.png", dpi=130); plt.close(fig)


def main():
    rows = run()
    write_csv(f"{HERE}/results/heat_path_sweep.csv", rows)
    hyd = []
    for gc in GC_LIST_MM:
        for v in V_LIST:
            c = coolant(gc / 1000, v)
            hyd.append(dict(gc_mm=gc, v_m_s=v, Dh_mm=2 * gc, **c, dP_per_m_kPa=c["dP_per_m_Pa"] / 1e3,
                            dP_per_m_pads_x2_kPa=2 * c["dP_per_m_Pa"] / 1e3, dP_core_kPa=c["dP_per_m_Pa"] * FLOW_LENGTH / 1e3))
    write_csv(f"{HERE}/results/coolant_hydraulics.csv", hyd)
    opt = [design_points(gc, v, kg, b, fp) for gc in GC_LIST_MM for v in V_LIST for kg in KG_LIST for b in BUDGETS for fp in PEAKING]
    write_csv(f"{HERE}/results/design_points.csv", opt)
    plots(rows)
    print(f"FLiBe @600C: rho={RHO_C:.1f} mu={MU_C*1e3:.3f} mPa.s Pr={PR_C:.2f}; pad fraction {F_PAD:.4f}")
    for h in hyd:
        print(h["gc_mm"], h["v_m_s"], f"Re={h['Re']:.0f} Nu={h['Nu']:.2f} h={h['h']:.0f} {h['regime']} dP/L={h['dP_per_m_kPa']:.3f} kPa/m")
    for o in opt:
        if o["peaking"] == 3.2 and o["budget_K"] == 50:
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in o.items()})


if __name__ == "__main__":
    main()
