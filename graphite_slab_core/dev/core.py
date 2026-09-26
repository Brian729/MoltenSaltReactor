# ---- core code (will be pasted into notebook cells) ----
import math, os, json, shutil, warnings
from dataclasses import dataclass, field, asdict
import numpy as np
import openmc

# ============================================================================
# 1. Shared machined-slot profile (used for BOTH fuel and coolant slots)
# ============================================================================
#   local 2-D coords of the slot cross-section:
#     u = across the slot width (0 .. w),  v = depth measured from the mouth (0 .. d)
#   round_location:
#     "both_sides" : (DEFAULT) both lateral side walls are quarter-rounds of radius R = d, centred on
#                    the MOUTH plane at u = d and u = w-d, with a flat bottom of width `flat` between
#                    them  ->  w = 2d + flat  (flat >= 0; flat = 0 gives a half-round).
#     "side"       : (legacy) only one side wall is a quarter-round R = d  ->  w = d + flat.
#     "bottom"     : (legacy) flat parallel sides, arc bottom R = d centred at mid-width on the mouth
#                    plane (U-groove); needs an explicit slot_width <= 2d.
#     "none"       : (legacy) plain rectangle, w = flat.
#   round_side ("side" only): "+" -> rounded wall at the high-u side, "-" -> low-u side.
PROFILES = ("both_sides", "side", "bottom", "none")

def slot_width_from(d, flat, round_location="both_sides"):
    """Total slot width at the mouth from depth and flat-bottom width."""
    if round_location == "both_sides":
        return 2.0 * d + flat
    if round_location == "side":
        return d + flat
    if round_location == "none":
        return flat
    raise ValueError("round_location='bottom' has no flat bottom: give slot_width explicitly")

def profile_area(w, d, round_location="both_sides"):
    if round_location == "both_sides":
        return (w - 2 * d) * d + math.pi * d**2 / 2.0
    if round_location == "side":
        return (w - d) * d + math.pi * d**2 / 4.0
    if round_location == "bottom":
        a = w / 2.0
        return a * math.sqrt(d**2 - a**2) + d**2 * math.asin(a / d)
    if round_location == "none":
        return w * d
    raise ValueError(round_location)

def profile_perimeter(w, d, round_location="both_sides"):
    """Full wetted perimeter (mouth included: the mouth is closed by graphite too)."""
    if round_location == "both_sides":
        return w + (w - 2 * d) + math.pi * d
    if round_location == "side":
        return w + (w - d) + d + math.pi * d / 2.0
    if round_location == "bottom":
        a = w / 2.0
        return w + 2.0 * math.sqrt(d**2 - a**2) + 2.0 * d * math.asin(a / d)
    if round_location == "none":
        return 2.0 * (w + d)
    raise ValueError(round_location)

def profile_feasible(w, d, round_location="both_sides"):
    tol = 1e-9
    if round_location not in PROFILES:
        return False, f"round_location must be one of {PROFILES}"
    if w <= 0 or d <= 0:
        return False, "slot width and slot_depth must be > 0"
    if round_location == "both_sides" and w < 2 * d - tol:
        return False, f"round_location='both_sides' needs slot_width >= 2*slot_depth (w={w:.3g} < 2d={2*d:.3g})"
    if round_location == "side" and w < d - tol:
        return False, f"round_location='side' needs slot_width >= slot_depth (w={w:.3g} < d={d:.3g})"
    if round_location == "bottom" and w > 2 * d + tol:
        return False, f"round_location='bottom' needs slot_width <= 2*slot_depth (w={w:.3g} > 2d={2*d:.3g})"
    return True, ""

def profile_outline(w, d, round_location="both_sides", round_side="+", n=60):
    """Closed (u, v) polyline of the profile - for sketches only."""
    if round_location == "none":
        pts = [(0, 0), (w, 0), (w, d), (0, d)]
    elif round_location == "both_sides":
        th = np.linspace(0, np.pi / 2, n)
        right = [(w - d + d * np.cos(t), d * np.sin(t)) for t in th]          # (w,0) -> (w-d,d)
        left = [(d - d * np.sin(t), d * np.cos(t)) for t in th]               # (d,d) -> (0,0)
        pts = [(0, 0)] + right + left
    elif round_location == "side":
        th = np.linspace(0, np.pi / 2, n)
        arc = [(w - d + d * np.cos(t), d * np.sin(t)) for t in th]            # (w,0) -> (w-d,d)
        pts = [(0, 0)] + arc + [(0, d)]
    else:  # bottom
        a = w / 2; h = math.sqrt(d**2 - a**2); t0 = math.asin(a / d)
        th = np.linspace(-t0, t0, n)
        arc = [(a + d * np.sin(t), d * np.cos(t)) for t in th]
        pts = [(0, 0)] + [(0, h)] + arc + [(w, h), (w, 0)]
    pts = np.array(pts, float)
    if round_side == "-" and round_location == "side":
        pts[:, 0] = w - pts[:, 0]
    return np.vstack([pts, pts[:1]])

_PLANE = {"x": openmc.XPlane, "y": openmc.YPlane, "z": openmc.ZPlane}

def _cylinder(axis, depth_axis, width_axis, c_depth, c_width, r):
    """Cylinder whose axis is the slot's run direction."""
    kw = {f"{depth_axis}0": c_depth, f"{width_axis}0": c_width, "r": r}
    return {"x": openmc.XCylinder, "y": openmc.YCylinder, "z": openmc.ZCylinder}[axis](**kw)

def milled_slot_region(run_axis, depth_axis, width_axis, mouth, u0, w, d,
                       round_location="both_sides", round_side="+", depth_sign=+1):
    """CSG region of ONE machined slot (infinite along run_axis).

    mouth      : coordinate of the open face along depth_axis
    depth_sign : +1 -> slot goes from `mouth` towards +depth_axis
    u0, w      : slot occupies width_axis in [u0, u0+w] at the mouth
    """
    ok, msg = profile_feasible(w, d, round_location)
    if not ok:
        raise ValueError(msg)
    P_d, P_w = _PLANE[depth_axis], _PLANE[width_axis]
    m = P_d(mouth)
    b = P_d(mouth + depth_sign * d)
    in_depth = (+m & -b) if depth_sign > 0 else (-m & +b)
    beyond_mouth = +m if depth_sign > 0 else -m        # half-space on the graphite side of the mouth
    lo, hi = u0, u0 + w
    cyl = lambda uc: _cylinder(run_axis, depth_axis, width_axis, mouth, uc, d)
    if round_location == "none":
        return in_depth & +P_w(lo) & -P_w(hi)
    if round_location == "both_sides":
        # [left quarter-disc] U [flat-bottom rectangle] U [right quarter-disc]; arc centres on the mouth plane
        pl, pr = P_w(lo + d), P_w(hi - d)
        left = -cyl(lo + d) & beyond_mouth & -pl
        right = -cyl(hi - d) & beyond_mouth & +pr
        if w - 2 * d > 1e-9:
            return left | (in_depth & +pl & -pr) | right
        return left | right                          # flat = 0: half-round
    if round_location == "side":
        if round_side == "+":
            uc = hi - d                     # arc centre (on mouth plane)
            rect = in_depth & +P_w(lo) & -P_w(uc)
            quarter = -cyl(uc) & beyond_mouth & +P_w(uc)
        else:
            uc = lo + d
            rect = in_depth & +P_w(uc) & -P_w(hi)
            quarter = -cyl(uc) & beyond_mouth & -P_w(uc)
        return rect | quarter
    if round_location == "bottom":
        uc = 0.5 * (lo + hi)
        return beyond_mouth & +P_w(lo) & -P_w(hi) & -cyl(uc)
    raise ValueError(round_location)

# ============================================================================
# 2. Materials (MSRE fuel salt per ORNL sources - see the markdown cell above)
# ============================================================================
import re
def _formula(compound):
    return {el: int(n) if n else 1 for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", compound)}

MSRE_U235_WT_PCT = 33.477      # U-235 wt% of total U at start of 235U power operation (ORNL-4658 Table 2.8, run 4-1)

DEFAULT_FUEL = dict(           # MSRE fuel salt, 235U operation (ORNL-4658, R.E. Thoma 1971)
    name="fuel salt (MSRE 7LiF-BeF2-ZrF4-UF4 65.0-29.17-5.0-0.83)",
    composition_mol={"LiF": 65.0, "BeF2": 29.17, "ZrF4": 5.0, "UF4": 0.83},   # mol%  (ORNL-4658 p.10)
    density_a=2.575, density_b=5.13e-4,   # rho[g/cc] = a - b*T[degC]  (ORNL-4658 Table 8.2) -> 2.242 at 649 C
    density=None,                         # set a number (g/cc) to override the correlation
    u235_wt_pct=MSRE_U235_WT_PCT,         # ENRICHMENT = U-235 WEIGHT % of total uranium
    u234_wt_pct=0.342, u236_wt_pct=0.141, # ORNL-4658 Table 2.8 (run 4-1 nominal); U-238 = balance (66.041)
    minor_u="scale",                      # 'scale': U-234/U-236 scale with enrichment (exact MSRE values at 33.477)
                                          # 'fixed': use the wt% above as given;  'none': U-235 + U-238 only
    li7_at_pct=99.995,                    # 7Li assay of MSRE fuel carrier salt (ORNL-4658 Table 2.11: 99.994-99.996)
)
DEFAULT_COOLANT = dict(        # MSRE coolant/flush salt  7LiF-BeF2 66-34 mol%
    name="coolant salt (MSRE 7LiF-BeF2 66-34)",
    composition_mol={"LiF": 66.0, "BeF2": 34.0},
    density_a=2.214, density_b=4.2e-4,    # ORNL-4658 Table 8.1 (Li2BeF4) -> 1.941 g/cc at 649 C
    density=None,
    li7_at_pct=99.992,                    # ORNL-4658 Table 2.1 (coolant/flush batches 99.991-99.992)
)
DEFAULT_GRAPHITE = dict(name="graphite (MSRE CGB, 1.87 g/cc)", density=1.87, sab="c_Graphite")
# 1.87 +/- 0.02 g/cc: MSRE grade-CGB graphite density used in the IRPhEP MSRE benchmark evaluation (Fratoni et al.).
# Pure carbon (no boron impurity) - placeholder for Brian's actual graphite grade.

def salt_density(spec, temperature_K):
    if spec.get("density") is not None:
        return float(spec["density"])
    return spec["density_a"] - spec["density_b"] * (temperature_K - 273.15)

def uranium_wt_fractions(u235_wt_pct, u234_wt_pct=0.0, u236_wt_pct=0.0, minor_u="scale"):
    """Return {nuclide: wt fraction of U}.  `u235_wt_pct` = enrichment in WEIGHT percent."""
    e = float(u235_wt_pct)
    if not 0.0 < e <= 100.0:
        raise ValueError(f"enrichment (U-235 wt%) must be in (0, 100], got {e}")
    if minor_u == "scale":
        w234, w236 = u234_wt_pct * e / MSRE_U235_WT_PCT, u236_wt_pct * e / MSRE_U235_WT_PCT
    elif minor_u == "fixed":
        w234, w236 = u234_wt_pct, u236_wt_pct
    elif minor_u == "none":
        w234 = w236 = 0.0
    else:
        raise ValueError("minor_u must be 'scale', 'fixed' or 'none'")
    w238 = 100.0 - e - w234 - w236
    if w238 < -1e-9:
        raise ValueError(f"U-234 + U-235 + U-236 exceed 100 wt% (enrichment {e}); use minor_u='none'")
    w = {"U234": w234, "U235": e, "U236": w236, "U238": max(w238, 0.0)}
    return {k: v / 100.0 for k, v in w.items() if v > 0}

_AMU = {"U234": 234.0409521, "U235": 235.0439299, "U236": 236.0455680, "U238": 238.0507882}

def make_salt(spec, temperature=922.0):
    """openmc.Material from mol% of compounds, isotopics and a density correlation."""
    atoms = {}
    tot = sum(spec["composition_mol"].values())
    for comp, x in spec["composition_mol"].items():
        for el, n in _formula(comp).items():
            atoms[el] = atoms.get(el, 0.0) + n * x / tot
    mat = openmc.Material(name=spec["name"], temperature=temperature)
    for el, a in atoms.items():
        if el == "Li":
            f7 = spec.get("li7_at_pct", 99.995) / 100.0
            mat.add_nuclide("Li7", a * f7); mat.add_nuclide("Li6", a * (1 - f7))
        elif el == "U":
            wf = uranium_wt_fractions(spec["u235_wt_pct"], spec.get("u234_wt_pct", 0.0),
                                      spec.get("u236_wt_pct", 0.0), spec.get("minor_u", "scale"))
            moles = {k: v / _AMU[k] for k, v in wf.items()}                 # wt -> atom fractions
            s = sum(moles.values())
            for k, m in moles.items():
                mat.add_nuclide(k, a * m / s)
        elif el in ("Be", "F", "Na"):                                     # mono-isotopic
            mat.add_nuclide({"Be": "Be9", "F": "F19", "Na": "Na23"}[el], a)
        else:
            mat.add_element(el, a)
    mat.set_density("g/cm3", salt_density(spec, temperature))
    return mat

def make_materials(fuel=None, coolant=None, graphite=None, temperature=922.0, enrichment=None):
    """enrichment: U-235 wt% of uranium (None -> MSRE value, 33.477 wt%)."""
    fuel = {**DEFAULT_FUEL, **(fuel or {})}
    if enrichment is not None:
        fuel["u235_wt_pct"] = float(enrichment)
    coolant = {**DEFAULT_COOLANT, **(coolant or {})}
    graphite = {**DEFAULT_GRAPHITE, **(graphite or {})}
    m_fuel = make_salt(fuel, temperature); m_fuel.depletable = True
    m_cool = make_salt(coolant, temperature)
    m_gr = openmc.Material(name=graphite["name"], temperature=temperature)
    m_gr.add_element("C", 1.0)
    m_gr.set_density("g/cm3", graphite["density"])
    if graphite.get("sab"):
        m_gr.add_s_alpha_beta(graphite["sab"])
    return {"fuel": m_fuel, "coolant": m_cool, "graphite": m_gr}

# ============================================================================
# 3. Layer stack / unit-cell dimensions (pure python - no OpenMC objects)
# ============================================================================
STACKINGS = ("plates", "interleaved")

def resolve_params(slot_depth=1.0, flat_width=0.5, web_thickness=1.5, wall_thickness=0.5, *,
                   round_location="both_sides", slot_width=None, stacking="plates", n_slot_pairs=2,
                   round_side="+", **_ignored):
    """Validate the parameters and return the unit-cell layout.

    Primary parameters: slot_depth d, flat_width (flat bottom between the two radii), web_thickness,
    wall_thickness.  Derived: slot_width w = 2d + flat (round_location='both_sides').  All slots
    (fuel and coolant) have the same size and profile; ONE web thickness is used everywhere.

    stacking="plates" (default, Brian's design): each graphite plate has one row of slots machined into
        one face (depth d) and a solid backing (wall_thickness).  Plates are stacked along x, every other
        plate rotated 90 deg, so fuel rows (slots along z) and coolant rows (slots along y) alternate and
        each row is closed by the backing of the next plate -> ONE shared wall between every fuel row and
        coolant row.  Layers along x: F | wall | C | wall.  The web is the land between neighbouring slots
        of a row (in-row).
    stacking="interleaved" (legacy): a slab holds n_slot_pairs x (F, web, C) rows separated by webs,
        slabs separated by one wall: F web C [web F web C]*(n-1) wall.

    `layers` entries are (kind, thickness), kind in {'F','C','web','wall'}; F/C layers are d thick.
    """
    d, t_web, t_wall = map(float, (slot_depth, web_thickness, wall_thickness))
    errs = []
    for k, v in dict(slot_depth=d, web_thickness=t_web, wall_thickness=t_wall).items():
        if not (v > 0 and math.isfinite(v)):
            errs.append(f"{k} must be a positive finite number (got {v})")
    if round_location not in PROFILES:
        errs.append(f"round_location must be one of {PROFILES}")
    if stacking not in STACKINGS:
        errs.append(f"stacking must be one of {STACKINGS}")
    if stacking == "interleaved" and (int(n_slot_pairs) != n_slot_pairs or n_slot_pairs < 1):
        errs.append("n_slot_pairs must be an integer >= 1")
    if round_side not in ("+", "-"):
        errs.append("round_side must be '+' or '-'")
    w = None
    if not errs:
        if slot_width is not None:
            w = float(slot_width)
        elif round_location == "bottom":
            errs.append("round_location='bottom' needs an explicit slot_width")
        else:
            f = float(flat_width)
            if not (f >= 0 and math.isfinite(f)):
                errs.append(f"flat_width must be >= 0 (got {flat_width})")
            else:
                w = slot_width_from(d, f, round_location)
    if not errs:
        ok, msg = profile_feasible(w, d, round_location)
        if not ok:
            errs.append(msg)
    if errs:
        raise ValueError("; ".join(errs))
    flat = {"both_sides": w - 2 * d, "side": w - d, "none": w, "bottom": 0.0}[round_location]

    if stacking == "plates":
        layers = [("F", d), ("wall", t_wall), ("C", d), ("wall", t_wall)]
    else:
        layers = [("F", d), ("web", t_web), ("C", d)]
        for _ in range(int(n_slot_pairs) - 1):
            layers += [("web", t_web), ("F", d), ("web", t_web), ("C", d)]
        layers += [("wall", t_wall)]
    return dict(slot_depth=d, flat_width=flat, slot_width=w, web_thickness=t_web, wall_thickness=t_wall,
                round_location=round_location, round_side=round_side, stacking=stacking,
                n_slot_pairs=int(n_slot_pairs) if stacking == "interleaved" else None, layers=layers,
                pitch_x=sum(t for _, t in layers), pitch_y=w + t_web, pitch_z=w + t_web)

def is_feasible(**kw):
    try:
        resolve_params(**kw); return True, ""
    except ValueError as e:
        return False, str(e)

def layer_start(p, kind, which=0):
    """x coordinate (unit cell centred on 0) where the `which`-th layer of `kind` starts (= slot mouth)."""
    x, n = -p["pitch_x"] / 2.0, 0
    for k, t in p["layers"]:
        if k == kind:
            if n == which:
                return x
            n += 1
        x += t
    raise ValueError(kind)

def analytic_metrics(core_radius=None, **kw):
    """Volume fractions etc. of the infinite periodic unit cell (exact, no transport).
    If core_radius is given, also the (approximate: fraction x volume) salt volumes in the H=2R cylinder."""
    p = resolve_params(**kw)
    w, d, loc = p["slot_width"], p["slot_depth"], p["round_location"]
    A, Pm = profile_area(w, d, loc), profile_perimeter(w, d, loc)
    nF = sum(k == "F" for k, _ in p["layers"]); nC = sum(k == "C" for k, _ in p["layers"])
    Px, Py, Pz = p["pitch_x"], p["pitch_y"], p["pitch_z"]
    V = Px * Py * Pz
    Vf = nF * A * Pz            # fuel slots run along z
    Vc = nC * A * Py            # coolant slots run along y
    Vg = V - Vf - Vc
    return dict(slot_width=w, pitch_x=Px, pitch_y=Py, pitch_z=Pz, slot_area=A, slot_perimeter=Pm,
                fuel_vf=Vf / V, coolant_vf=Vc / V, graphite_vf=Vg / V,
                graphite_to_fuel=Vg / Vf, coolant_to_fuel=Vc / Vf,
                fuel_hydraulic_diam=4 * A / Pm,
                fuel_wetted_area_per_fuel_vol=Pm / A,      # cm^2 / cm^3 (same for coolant)
                plate_thickness=(d + p["wall_thickness"]) if p["stacking"] == "plates" else Px - p["wall_thickness"],
                **({} if core_radius is None else dict(
                    core_volume_l=2 * math.pi * core_radius**3 / 1000.0,
                    core_fuel_volume_l=Vf / V * 2 * math.pi * core_radius**3 / 1000.0,
                    core_coolant_volume_l=Vc / V * 2 * math.pi * core_radius**3 / 1000.0)))

# ============================================================================
# 4. OpenMC geometry
# ============================================================================
def build_unit_universe(p, mats):
    """Universe of one unit cell, centred on the origin, with UNBOUNDED cells so it can be
    used both as a periodic unit cell and as a lattice element."""
    w, d, loc, web = p["slot_width"], p["slot_depth"], p["round_location"], p["web_thickness"]
    Px, Py, Pz = p["pitch_x"], p["pitch_y"], p["pitch_z"]
    x = -Px / 2.0
    cells, slots = [], []
    for i, (kind, t) in enumerate(p["layers"]):
        if kind == "F":
            # vertical fuel slot: runs along z, width along y, mouth at the layer's -x face
            r = milled_slot_region("z", "x", "y", mouth=x, u0=-Py / 2 + web / 2, w=w, d=d,
                                   round_location=loc, round_side=p["round_side"])
            c = openmc.Cell(name=f"fuel slot (layer {i})", fill=mats["fuel"], region=r)
        elif kind == "C":
            # horizontal coolant slot: runs along y, width along z, mouth at the layer's -x face
            r = milled_slot_region("y", "x", "z", mouth=x, u0=-Pz / 2 + web / 2, w=w, d=d,
                                   round_location=loc, round_side=p["round_side"])
            c = openmc.Cell(name=f"coolant slot (layer {i})", fill=mats["coolant"], region=r)
        else:
            c = None
        if c is not None:
            cells.append(c); slots.append(r)
        x += t
    gr_region = ~openmc.Union(slots)
    cells.append(openmc.Cell(name="graphite (webs + walls)", fill=mats["graphite"], region=gr_region))
    return openmc.Universe(name="slab unit cell", cells=cells)

def _box(Lx, Ly, Lz, bc="transmission", center=(0, 0, 0)):
    cx, cy, cz = center
    s = [openmc.XPlane(cx - Lx / 2, boundary_type=bc), openmc.XPlane(cx + Lx / 2, boundary_type=bc),
         openmc.YPlane(cy - Ly / 2, boundary_type=bc), openmc.YPlane(cy + Ly / 2, boundary_type=bc),
         openmc.ZPlane(cz - Lz / 2, boundary_type=bc), openmc.ZPlane(cz + Lz / 2, boundary_type=bc)]
    return s, (+s[0] & -s[1] & +s[2] & -s[3] & +s[4] & -s[5])

def build_model(slot_depth=1.0, flat_width=0.5, web_thickness=1.5, wall_thickness=0.5, *,
                round_location="both_sides", slot_width=None, stacking="plates", n_slot_pairs=2, round_side="+",
                mode="cylinder", core_radius=70.0, reflector_thickness=0.0,
                enrichment=None, fuel=None, coolant=None, graphite=None, temperature=922.0,
                particles=10000, batches=100, inactive=40, seed=1):
    """Return an openmc.Model of the slotted-graphite-plate core.

    Geometry parameters (cm): slot_depth d, flat_width (slot width w = 2d + flat_width), web_thickness,
    wall_thickness (shared plate backing between fuel and coolant rows).  See resolve_params().
    mode="cylinder"  : (default) the plate stack (RectLattice of unit cells) fills a right cylinder of radius
                       R = core_radius and height H = 2R (z in [-R, R]); bare by default (reflector_thickness=0);
                       vacuum outside  -> finite k-eff.
    mode="unit_cell" : one unit cell with PERIODIC boundaries on all six faces -> k-infinity (quick mode).
    enrichment       : U-235 wt% of total uranium in the fuel salt (None -> MSRE 33.477 wt%).
    """
    p = resolve_params(slot_depth, flat_width, web_thickness, wall_thickness,
                       round_location=round_location, slot_width=slot_width, stacking=stacking,
                       n_slot_pairs=n_slot_pairs, round_side=round_side)
    if mode not in ("cylinder", "unit_cell"):
        raise ValueError("mode must be 'cylinder' or 'unit_cell'")
    mats = make_materials(fuel, coolant, graphite, temperature, enrichment)
    univ = build_unit_universe(p, mats)
    Px, Py, Pz = p["pitch_x"], p["pitch_y"], p["pitch_z"]

    if mode == "unit_cell":
        s, region = _box(Px, Py, Pz, bc="periodic")
        s[0].periodic_surface = s[1]; s[2].periodic_surface = s[3]; s[4].periodic_surface = s[5]
        root = openmc.Universe(cells=[openmc.Cell(name="unit cell", fill=univ, region=region)])
        source = openmc.IndependentSource(space=openmc.stats.Box((-Px / 2, -Py / 2, -Pz / 2), (Px / 2, Py / 2, Pz / 2)),
                                          constraints={"fissionable": True})
    else:
        R, t = float(core_radius), float(reflector_thickness)
        if not R > 0 or t < 0:
            raise ValueError("core_radius must be > 0 and reflector_thickness >= 0")
        if R < 2 * max(Px, Py, Pz):
            raise ValueError(f"core_radius {R} too small for the unit-cell pitch ({Px:.2f} x {Py:.2f} x {Pz:.2f} cm)")
        # lattice just large enough to cover the cylinder; an odd count centres a unit cell on the axis
        n = [2 * int(math.ceil(R / P - 0.5)) + 1 for P in (Px, Py, Pz)]
        lat = openmc.RectLattice(name="slab stack")
        lat.pitch = (Px, Py, Pz)
        lat.lower_left = (-n[0] * Px / 2, -n[1] * Py / 2, -n[2] * Pz / 2)
        lat.universes = np.full((n[2], n[1], n[0]), univ)          # [z][y][x]
        lat.outer = openmc.Universe(cells=[openmc.Cell(fill=mats["graphite"])])
        cyl = openmc.ZCylinder(r=R); zlo = openmc.ZPlane(-R); zhi = openmc.ZPlane(R)
        core = -cyl & +zlo & -zhi
        cells = [openmc.Cell(name="slab-stack core", fill=lat, region=core)]
        if t > 0:
            cyl_o = openmc.ZCylinder(r=R + t, boundary_type="vacuum")
            zlo_o = openmc.ZPlane(-R - t, boundary_type="vacuum"); zhi_o = openmc.ZPlane(R + t, boundary_type="vacuum")
            cells.append(openmc.Cell(name="graphite reflector", fill=mats["graphite"],
                                     region=-cyl_o & +zlo_o & -zhi_o & ~core))
        else:
            for srf in (cyl, zlo, zhi):
                srf.boundary_type = "vacuum"
        root = openmc.Universe(cells=cells)
        p.update(core_radius=R, core_height=2 * R, reflector_thickness=t, lattice_shape=tuple(n))
        source = openmc.IndependentSource(
            space=openmc.stats.CylindricalIndependent(r=openmc.stats.PowerLaw(0.0, R, 1.0),
                                                      phi=openmc.stats.Uniform(0.0, 2 * math.pi),
                                                      z=openmc.stats.Uniform(-R, R)),
            constraints={"fissionable": True})

    geometry = openmc.Geometry(root)
    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles, settings.batches, settings.inactive = int(particles), int(batches), int(inactive)
    settings.seed = int(seed)
    settings.temperature = {"method": "nearest", "tolerance": 150.0, "default": temperature}
    settings.source = source
    settings.output = {"tallies": False}
    model = openmc.Model(geometry=geometry, materials=openmc.Materials(mats.values()), settings=settings)
    model.params = {**p, "mode": mode}   # handy for bookkeeping
    model.mats = mats
    return model
