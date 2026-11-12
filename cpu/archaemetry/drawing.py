"""Stage 4: the publication drawing.

Follows the usual conventions for pottery profiles (see e.g. the user's reference
figures in docs/): the vessel stands on its base with the rim at the top; a thin
vertical line marks the axis; on one side the wall section is drawn and filled solid
black; on the other side the exterior outline is drawn over the preserved height; a
horizontal rim line runs across the full diameter; the sherd ID sits above; a block
scale bar in centimetres sits below. Everything is at 1:1 in millimetres, so the SVG
and PDF print at true size.
"""
import numpy as np

FONT = "Helvetica, Arial, sans-serif"


def layout(section, kind, ext_h, ext_r, side="left"):
    """Turn (r, h) millimetre geometry into drawing coordinates (x right, y down,
    origin at the axis on the rim line). Returns a dict of polylines."""
    s = -1.0 if side == "left" else 1.0
    top = max(section[:, 1].max(), ext_h.max())
    sec = np.column_stack([s * section[:, 0], top - section[:, 1]])
    ext = np.column_stack([-s * ext_r, top - ext_h])
    # complete vessel (section reaches the axis at the base): close the exterior view
    # along the underside to the axis, as in a whole-vessel drawing
    bottom = section[:, 1].min()
    if kind == "section" and section[np.abs(section[:, 1] - bottom) < 1.0, 0].min() < 1.0:
        ext = np.vstack([ext, [ext[-1, 0], top - bottom], [0.0, top - bottom]])
    # rim radius: outermost section radius within 1.5 mm of the top
    near_top = section[section[:, 1] >= top - 1.5]
    r_rim = float(near_top[:, 0].max()) if len(near_top) else float(section[:, 0].max())
    r_ext_top = float(ext_r[np.argmax(ext_h)])
    depth = float(max(sec[:, 1].max(), ext[:, 1].max()))
    return dict(section=sec, kind=kind, exterior=ext, side=s, rim_y=0.0,
                rim_x=(s * r_rim, -s * r_ext_top), r_rim=r_rim, depth=depth)


def _scale_bar_len(width_mm):
    for cm in (1, 2, 3, 5, 10, 15, 20, 30):
        if cm * 10 >= width_mm / 3:
            return cm
    return 30


def to_svg(G, path, sherd_id="", caption="", fill="solid", line_mm=0.3, scale_bar=True):
    sec, ext = G["section"], G["exterior"]
    xs = np.concatenate([sec[:, 0], ext[:, 0], list(G["rim_x"])])
    x0, x1 = xs.min() - 8, xs.max() + 8
    y0, y1 = -14.0, G["depth"] + (24 if scale_bar else 10) + (6 if caption else 0)
    W, H = x1 - x0, y1 - y0
    pts = lambda P: " ".join(f"{x:.3f},{y:.3f}" for x, y in P)
    el = []
    el.append('<defs><pattern id="hatch" patternUnits="userSpaceOnUse" width="1.2" height="1.2" '
              'patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="1.2" stroke="black" '
              f'stroke-width="{line_mm * 0.6:.3f}"/></pattern></defs>')
    if G["kind"] == "section":
        f = {"solid": "black", "hatch": "url(#hatch)", "none": "none"}[fill]
        el.append(f'<polygon id="section" fill="{f}" stroke="black" stroke-width="{line_mm}" '
                  f'stroke-linejoin="round" points="{pts(sec)}"/>')
    else:
        el.append(f'<polyline id="section" fill="none" stroke="black" stroke-width="{line_mm}" '
                  f'points="{pts(sec)}"/>')
    el.append(f'<polyline id="exterior" fill="none" stroke="black" stroke-width="{line_mm}" '
              f'points="{pts(ext)}"/>')
    rx0, rx1 = G["rim_x"]
    el.append(f'<line id="rim" x1="{rx0:.3f}" y1="0" x2="{rx1:.3f}" y2="0" stroke="black" '
              f'stroke-width="{line_mm}"/>')
    el.append(f'<line id="axis" x1="0" y1="0" x2="0" y2="{G["depth"] + 3:.3f}" stroke="black" '
              f'stroke-width="{line_mm * 0.6:.3f}"/>')
    if sherd_id:
        el.append(f'<text x="{min(rx0, rx1):.2f}" y="-4" font-family="{FONT}" font-size="3.5">{_esc(sherd_id)}</text>')
    yb = G["depth"] + 12
    if scale_bar:
        cm = _scale_bar_len(W)
        bx = x0 + 4
        for i in range(cm):
            el.append(f'<rect x="{bx + 10 * i:.2f}" y="{yb:.2f}" width="10" height="2" '
                      f'fill="{"black" if i % 2 == 0 else "white"}" stroke="black" stroke-width="0.2"/>')
        for i in range(cm + 1):
            if cm <= 5 or i % (cm // 5 if cm // 5 else 1) == 0 or i == cm:
                el.append(f'<text x="{bx + 10 * i:.2f}" y="{yb + 6:.2f}" text-anchor="middle" '
                          f'font-family="{FONT}" font-size="2.6">{i}</text>')
        el.append(f'<text x="{bx + 10 * cm + 3:.2f}" y="{yb + 2:.2f}" font-family="{FONT}" font-size="2.6">cm</text>')
    if caption:
        el.append(f'<text x="{x0 + 4:.2f}" y="{yb + (12 if scale_bar else 2):.2f}" font-family="{FONT}" '
                  f'font-size="2.4" fill="#555">{_esc(caption)}</text>')
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W:.2f}mm" height="{H:.2f}mm" '
           f'viewBox="{x0:.3f} {y0:.3f} {W:.3f} {H:.3f}">\n<title>{_esc(sherd_id or "profile")}</title>\n'
           + "\n".join(el) + "\n</svg>\n")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(svg)


def to_raster(G, paths, sherd_id="", caption="", fill="solid", line_mm=0.3, scale_bar=True, dpi=300):
    """PNG/PDF at 1:1 (figure size in inches = mm / 25.4)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle
    sec, ext = G["section"], G["exterior"]
    xs = np.concatenate([sec[:, 0], ext[:, 0], list(G["rim_x"])])
    x0, x1 = xs.min() - 8, xs.max() + 8
    y0, y1 = -14.0, G["depth"] + (24 if scale_bar else 10) + (6 if caption else 0)
    fig = plt.figure(figsize=((x1 - x0) / 25.4, (y1 - y0) / 25.4))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(x0, x1); ax.set_ylim(y1, y0); ax.axis("off")
    lw = line_mm / 25.4 * 72
    if G["kind"] == "section":
        ax.add_patch(Polygon(sec, closed=True, facecolor="black" if fill == "solid" else "white",
                             hatch="////" if fill == "hatch" else None, edgecolor="black", lw=lw))
    else:
        ax.plot(sec[:, 0], sec[:, 1], "k-", lw=lw)
    ax.plot(ext[:, 0], ext[:, 1], "k-", lw=lw)
    ax.plot(G["rim_x"], [0, 0], "k-", lw=lw)
    ax.plot([0, 0], [0, G["depth"] + 3], "k-", lw=lw * 0.6)
    pt = 72 / 25.4
    if sherd_id:
        ax.text(min(G["rim_x"]), -4, sherd_id, fontsize=3.5 * pt)
    yb = G["depth"] + 12
    if scale_bar:
        cm = _scale_bar_len(x1 - x0); bx = x0 + 4
        for i in range(cm):
            ax.add_patch(Rectangle((bx + 10 * i, yb), 10, 2, facecolor="black" if i % 2 == 0 else "white",
                                   edgecolor="black", lw=0.2 * pt))
        for i in range(cm + 1):
            if cm <= 5 or i % (cm // 5 if cm // 5 else 1) == 0 or i == cm:
                ax.text(bx + 10 * i, yb + 6, str(i), ha="center", fontsize=2.6 * pt)
        ax.text(bx + 10 * cm + 3, yb + 2, "cm", fontsize=2.6 * pt)
    if caption:
        ax.text(x0 + 4, yb + (12 if scale_bar else 2), caption, fontsize=2.4 * pt, color="#555")
    for p in paths:
        fig.savefig(p, dpi=dpi, facecolor="white")
    plt.close(fig)


def qa_figure(path, h, r, bp, G_section_rh, thickness, title=""):
    """Diagnostic figure for the user: every vertex in (r, h), the chosen cut, the robust
    exterior, and thickness down the wall."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    top = max(h.max(), G_section_rh[:, 1].max())
    fig, ax = plt.subplots(1, 2, figsize=(10, 6), gridspec_kw={"width_ratios": [2, 1]})
    k = max(1, len(h) // 40000)
    ax[0].scatter(r[::k], top - h[::k], s=0.3, c="#c9b79c", label="all vertices (r, depth)")
    ax[0].plot(bp["r_out"], top - bp["h"], "-", c="#3b2a1e", lw=1, label="robust exterior (median)")
    P = G_section_rh
    ax[0].plot(np.r_[P[:, 0], P[:1, 0]], np.r_[top - P[:, 1], top - P[:1, 1]], "-", c="#c2561d", lw=1,
               label="drawn section cut")
    ax[0].invert_yaxis(); ax[0].set_aspect("equal"); ax[0].set_xlabel("radius (mm)")
    ax[0].set_ylabel("depth below rim (mm)"); ax[0].legend(fontsize=7, loc="lower left"); ax[0].set_title(title)
    if thickness:
        dd, tt = zip(*thickness)
        ax[1].plot(tt, dd, "o-", c="#c2561d", ms=3)
        ax[1].set_xlim(0, max(tt) * 1.3)
    ax[1].invert_yaxis(); ax[1].set_xlabel("wall thickness (mm)"); ax[1].set_title("perpendicular thickness")
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
