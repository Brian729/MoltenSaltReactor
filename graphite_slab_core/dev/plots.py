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

def sketch_profiles(d=1.0, flat=0.5, filename=None):
    """Annotated cross-sections of the shared machined-slot profile: the default double-rounded slot
    (large panel) and the legacy options (small panels)."""
    w = 2 * d + flat
    fig = plt.figure(figsize=(14, 5.2))
    gs = fig.add_gridspec(2, 4, width_ratios=[2.2, 2.2, 1, 1], height_ratios=[1, 1])
    ax0 = fig.add_subplot(gs[:, :2])
    small = [fig.add_subplot(gs[0, 2]), fig.add_subplot(gs[0, 3]), fig.add_subplot(gs[1, 2]), fig.add_subplot(gs[1, 3])]

    def draw(ax, loc, ww, title, big=False):
        ax.add_patch(plt.Rectangle((-0.6, 0), ww + 1.2, d + 0.7, color="0.75", zorder=0))
        o = profile_outline(ww, d, loc)
        ax.fill(o[:, 0], o[:, 1], color="#e67814", zorder=1)
        ax.plot(o[:, 0], o[:, 1], "k", lw=1)
        ax.axhline(0, color="k", lw=2.5)
        ax.set_xlim(-0.6, ww + 0.6); ax.set_ylim(d + 0.7, -0.55); ax.set_aspect("equal")
        ax.set_title(title, fontsize=10 if big else 8)
        if not big:
            ax.set_xticks([]); ax.set_yticks([])

    draw(ax0, "both_sides", w, f"DEFAULT  round_location='both_sides':  w = 2d + flat = {w:g} cm\n"
         f"(d = {d:g} cm, flat = {flat:g} cm; both side walls are quarter-rounds R = d centred on the mouth plane)", big=True)
    for uc, ang in [(d, 125), (w - d, 55)]:
        ax0.plot([uc], [0], "k+", ms=12)
        a = np.deg2rad(ang)
        ax0.annotate("", xy=(uc + d * np.cos(a), d * np.sin(a)), xytext=(uc, 0), arrowprops=dict(arrowstyle="->"))
        ax0.text(uc + 0.55 * d * np.cos(a), 0.5 * d * np.sin(a), "R = d", fontsize=9, ha="center",
                 bbox=dict(fc="#e67814", ec="none", pad=0.5))
    ax0.text(w / 2, -0.08, "mouth (open face, closed by the backing wall of the next plate)", ha="center", va="bottom", fontsize=8)
    ax0.annotate("", xy=(0, -0.38), xytext=(w, -0.38), arrowprops=dict(arrowstyle="<->"))
    ax0.text(w / 2, -0.42, "w = slot width (derived)", ha="center", va="bottom", fontsize=8)
    if flat > 0:
        ax0.annotate("", xy=(d, d + 0.2), xytext=(w - d, d + 0.2), arrowprops=dict(arrowstyle="<->"))
        ax0.text(w / 2, d + 0.25, "flat_width", ha="center", va="top", fontsize=8)
    ax0.annotate("", xy=(w + 0.35, 0), xytext=(w + 0.35, d), arrowprops=dict(arrowstyle="<->"))
    ax0.text(w + 0.4, d / 2, "d = slot_depth", rotation=90, va="center", fontsize=8)
    ax0.set_xlabel("u  (across width: y for fuel slots, z for coolant slots) [cm]", fontsize=8)
    ax0.set_ylabel("v = depth from mouth (+x) [cm]", fontsize=8)

    draw(small[0], "both_sides", 2 * d, "both_sides, flat = 0\n(half-round)")
    draw(small[1], "side", d + flat + d / 2, "legacy 'side'\n(one side rounded)")
    draw(small[2], "bottom", 2 * d * 0.9, "legacy 'bottom'\n(U-groove, w <= 2d)")
    draw(small[3], "none", w, "legacy 'none'\n(square)")
    fig.suptitle("Shared machined-slot profile (identical for fuel slots, extruded along z, and coolant slots, extruded along y)",
                 fontsize=10)
    fig.tight_layout()
    if filename:
        fig.savefig(filename, dpi=150)
    return fig
