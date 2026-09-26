import os, warnings, math, textwrap
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import openmc

MAT_COLORS = {"fuel": (230, 120, 20), "coolant": (40, 110, 220), "graphite": (150, 150, 150)}

def have_cross_sections():
    xs = openmc.config.get("cross_sections")
    return xs is not None and os.path.exists(str(xs))

def _color_map(model):
    return {model.mats[k]: MAT_COLORS[k] for k in model.mats}

def plot_slice_python(model, basis="xy", origin=(0, 0, 0), width=(10, 10), pixels=(300, 300), ax=None):
    """Nuclear-data-free fallback: colour each pixel by geometry.find() (pure-Python CSG)."""
    ax = ax or plt.subplots(figsize=(7, 7 * width[1] / width[0]))[1]
    i, j = {"xy": (0, 1), "xz": (0, 2), "yz": (1, 2)}[basis]
    nx, ny = pixels
    cmap = {m.id: np.array(c) / 255 for m, c in _color_map(model).items()}
    img = np.ones((ny, nx, 3))
    us = origin[i] - width[0] / 2 + (np.arange(nx) + 0.5) * width[0] / nx
    vs = origin[j] + width[1] / 2 - (np.arange(ny) + 0.5) * width[1] / ny
    for r, v in enumerate(vs):
        for c, u in enumerate(us):
            pt = list(origin); pt[i] = u; pt[j] = v
            path = model.geometry.find(pt)
            if path and isinstance(path[-1], openmc.Cell) and isinstance(path[-1].fill, openmc.Material):
                img[r, c] = cmap[path[-1].fill.id]
    ext = [origin[i] - width[0] / 2, origin[i] + width[0] / 2, origin[j] - width[1] / 2, origin[j] + width[1] / 2]
    ax.imshow(img, extent=ext, origin="upper", interpolation="nearest")
    ax.set_xlabel(f"{basis[0]} [cm]"); ax.set_ylabel(f"{basis[1]} [cm]")
    return ax

def plot_slice(model, basis="xy", origin=(0, 0, 0), width=(10, 10), pixels=(600, 600),
               title=None, filename=None, force_python=False):
    """Slice plot coloured by material. Uses openmc's plotter (needs cross_sections.xml);
    falls back to pure-Python point sampling if nuclear data are not configured."""
    aspect = width[1] / width[0]
    fig, ax = plt.subplots(figsize=(9, max(3.2, 9 * aspect) + 1.2))
    used = "openmc"
    if have_cross_sections() and not force_python:
        try:
            model.plot(basis=basis, origin=origin, width=width, pixels=pixels, color_by="material",
                       colors=_color_map(model), axes=ax)
        except Exception as e:                        # pragma: no cover
            warnings.warn(f"openmc plot failed ({e}); using python fallback"); used = "python"
    else:
        used = "python"
    if used == "python":
        px = (min(pixels[0], 400), max(1, int(min(pixels[0], 400) * width[1] / width[0])))
        plot_slice_python(model, basis, origin, width, px, ax=ax)
    ax.set_xlabel(f"{basis[0]} [cm]"); ax.set_ylabel(f"{basis[1]} [cm]")
    ax.legend(handles=[Patch(color=np.array(c) / 255, label=model.mats[k].name) for k, c in MAT_COLORS.items()],
              loc="upper center", bbox_to_anchor=(0.5, -0.1 if aspect < 0.8 else -0.07), ncol=2, fontsize=8, frameon=False)
    ax.set_title("\n".join(textwrap.wrap((title or f"{basis} slice at {tuple(round(o, 3) for o in origin)}")
                                         + f"  [{used} renderer]", 90)), fontsize=9)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=150)
    return fig, ax

def csg_volume_check(model, n=20000, seed=0):
    """Monte-Carlo point sampling of the (periodic unit-cell) model with pure-Python
    geometry.find(); returns sampled volume fractions of fuel / coolant / graphite."""
    rng = np.random.default_rng(seed)
    p = model.params
    L = np.array([p["pitch_x"], p["pitch_y"], p["pitch_z"]])
    pts = (rng.random((n, 3)) - 0.5) * L
    ids = {m.id: k for k, m in model.mats.items()}
    counts = dict.fromkeys(model.mats, 0)
    for pt in pts:
        cell = model.geometry.find(pt)[-1]
        counts[ids[cell.fill.id]] += 1
    return {f"{k}_vf": (c / n, math.sqrt(c / n * (1 - c / n) / n)) for k, c in counts.items()}

def sketch_profiles(w=2.0, d=1.5, round_side="+", filename=None):
    """Annotated cross-sections of the shared milled-slot profile for each round_location."""
    locs = [("side", "DEFAULT  round_location='side'\n(one side wall = quarter-round R=d)"),
            ("bottom", "ALT  round_location='bottom'\n(U-groove, arc bottom R=d, needs w<=2d)"),
            ("none", "REF  round_location='none'\n(square slot)")]
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, (loc, ttl) in zip(axs, locs):
        ok, msg = profile_feasible(w, d, loc)
        ax.add_patch(plt.Rectangle((-0.8, 0), w + 1.6, d + 0.8, color="0.75", zorder=0))
        if ok:
            o = profile_outline(w, d, loc, round_side)
            ax.fill(o[:, 0], o[:, 1], color="#e67814", zorder=1)
            ax.plot(o[:, 0], o[:, 1], "k", lw=1)
            if loc != "none":
                uc = (w - d if round_side == "+" else d) if loc == "side" else w / 2
                ax.plot([uc], [0], "k+", ms=10)
                ang = np.deg2rad(55) if round_side == "+" or loc == "bottom" else np.deg2rad(125)
                if loc == "bottom": ang = np.deg2rad(90 + 20)
                ax.annotate("", xy=(uc + d * np.cos(ang), d * np.sin(ang)), xytext=(uc, 0),
                            arrowprops=dict(arrowstyle="->"))
                ax.text(uc + 0.5 * d * np.cos(ang) + 0.05, 0.5 * d * np.sin(ang), "R = d", fontsize=9)
        else:
            ax.text(w / 2, d / 2, "infeasible:\n" + msg, ha="center", va="center", fontsize=8, wrap=True)
        ax.axhline(0, color="k", lw=2.5)
        ax.text(w / 2, -0.08, "mouth (open face, closed by neighbouring graphite)", ha="center", va="bottom", fontsize=8)
        ax.annotate("", xy=(0, d + 0.35), xytext=(w, d + 0.35), arrowprops=dict(arrowstyle="<->"))
        ax.text(w / 2, d + 0.42, "w = slot_width", ha="center", fontsize=8)
        ax.annotate("", xy=(w + 0.45, 0), xytext=(w + 0.45, d), arrowprops=dict(arrowstyle="<->"))
        ax.text(w + 0.5, d / 2, "d = slot_depth", rotation=90, va="center", fontsize=8)
        ax.set_xlim(-0.8, w + 0.8); ax.set_ylim(d + 0.8, -0.6); ax.set_aspect("equal")
        ax.set_xlabel("u  (across width: y for fuel slots, z for coolant slots) [cm]", fontsize=8)
        ax.set_ylabel("v = depth from mouth (+x) [cm]", fontsize=8)
        ax.set_title(ttl, fontsize=9)
    fig.suptitle(f"Shared milled-slot profile, w={w} cm, d={d} cm, round_side='{round_side}' "
                 "(fuel slot: extruded along z; coolant slot: extruded along y)", fontsize=10)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=150)
    return fig
