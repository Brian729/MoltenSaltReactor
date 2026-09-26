#!/usr/bin/env python3
"""Metrics for one fuel-slot case at a given (default latest) time.

peak fuel T, location; max over z of (T_fuel_max(z) - T_interface(z)); max |U|; coolant-to-peak dT.
"""
import json, sys
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__file__))
from foamio import read_mesh_centres, read_internal, read_patch, latest_time


def metrics(case, t=None):
    t = t or latest_time(case)
    p = json.load(open(f"{case}/caseParams.json"))
    C, P = read_mesh_centres(case, "fuel")
    T = read_internal(f"{case}/{t}/fuel/T", len(C))
    U = read_internal(f"{case}/{t}/fuel/U", len(C))
    Tw = read_patch(f"{case}/{t}/fuel/T", "fuel_to_graphite")
    zf = P["fuel_to_graphite"][:, 2]
    Cg, Pg = read_mesh_centres(case, "graphite")
    Tg = read_internal(f"{case}/{t}/graphite/T", len(Cg))
    umag = np.linalg.norm(U, axis=1)
    # per-row statistics (structured mesh: group by z)
    zr = np.round(C[:, 2], 7)
    rows = np.unique(zr)
    Tmax_row = np.array([T[zr == z].max() for z in rows])
    zw = np.round(zf, 7)
    Tw_row = np.array([Tw[zw == z].max() for z in rows])   # both walls; use hotter wall (conservative dT uses min? see below)
    Tw_row_min = np.array([Tw[zw == z].min() for z in rows])
    dT_row = Tmax_row - Tw_row
    ipk = int(np.argmax(T))
    return dict(time=float(t), T_peak_C=float(T.max() - 273.15), z_peak=float(C[ipk, 2]), x_peak=float(C[ipk, 0]),
                dT_fuel_wall_max=float(dT_row.max()), z_dTmax=float(rows[np.argmax(dT_row)]),
                wall_asym_max=float((Tw_row - Tw_row_min).max()),
                T_wall_max_C=float(Tw.max() - 273.15), T_graphite_max_C=float(Tg.max() - 273.15),
                U_max=float(umag.max()), Uz_max=float(U[:, 2].max()), Uz_min=float(U[:, 2].min()),
                T_fuel_mean_C=float(T.mean() - 273.15))


if __name__ == "__main__":
    for c in sys.argv[1:]:
        print(c, json.dumps(metrics(c)))
