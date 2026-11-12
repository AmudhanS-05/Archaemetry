"""Stage 3: from a mesh and its axis to the drawn profile.

Two complementary views, both in (r, h) millimetres with h measured down from the rim:

band profile   every vertex binned by height; the median radius per band is the
               exterior line. Uses the whole surface, so it is robust to noise and
               ignores handles (a minority of the circumference).
section        the mesh cut by a plane containing the axis, as a draughtsperson or
               GigaMesh does. For a mesh with both wall surfaces this is a closed
               polygon (outer wall, lip, inner wall, break) that is drawn filled; wall
               thickness is measured on it perpendicular to the wall, like a calliper.
"""
import numpy as np
from .axis import basis_from_dir, fit_planar_circle, border_vertices


# ---------------------------------------------------------------- frame
def frame(V, p0, d):
    """Cylindrical coordinates: h along d, r from axis, angle in the e1/e2 plane."""
    e1, e2, d = basis_from_dir(d)
    W = V - p0
    h = W @ d
    u, v = W @ e1, W @ e2
    return h, np.hypot(u, v), np.arctan2(v, u), (e1, e2, d)


def orient_axis(V, T, p0, d, mode="auto"):
    """Choose the sign of d so +d points toward the rim. Returns (d, how).

    auto: if the mesh has open borders (a sherd, or a vessel with an open mouth), the
    rim is the end whose border is the cleanest planar circle (a rim is round, a break
    is ragged). A closed mesh has no such evidence: keep +d closest to world +Y, the
    usual 'up' after orienting in MeshLab. 'flip' reverses whatever auto decided."""
    d = d / np.linalg.norm(d)
    how = "world_up"
    bidx = border_vertices(T) if len(T) else np.array([], int)
    if len(bidx) >= 24:
        B = V[bidx]
        h = (B - p0) @ d
        top, bot = B[h >= np.quantile(h, 0.85)], B[h <= np.quantile(h, 0.15)]
        if len(top) >= 12 and len(bot) >= 12:
            nt = fit_planar_circle(top)[3]
            nb = fit_planar_circle(bot)[3]
            if abs(nt - nb) > 0.2 * max(nt, nb, 1e-9):
                how = "rim_border"
                if nb < nt:
                    d = -d
    if how == "world_up":
        # closed mesh: a rim is level all round the arc, a break is ragged in height
        e1, e2, _ = basis_from_dir(d)
        W = V - p0
        h = W @ d
        a = np.arctan2(W @ e2, W @ e1)
        c = np.angle(np.mean(np.exp(1j * a)))
        rel = (a - c + np.pi) % (2 * np.pi) - np.pi
        sec = np.floor((rel - rel.min()) / (np.ptp(rel) + 1e-12) * 12).clip(0, 11)
        tops = [h[sec == k].max() for k in range(12) if (sec == k).sum() > 20]
        bots = [h[sec == k].min() for k in range(12) if (sec == k).sum() > 20]
        if len(tops) >= 6:
            st, sb = np.std(tops), np.std(bots)
            if max(st, sb) > 2.5 * max(min(st, sb), 0.05):
                how = "level_end"
                if st > sb:
                    d = -d
    if how == "world_up":
        up = 1 if abs(d[1]) >= abs(d[2]) else 2      # +Y up (MeshLab), else +Z up
        if d[up] < 0:
            d = -d
    if mode == "flip":
        d = -d
        how += "+flipped"
    return d, how


# ---------------------------------------------------------------- band profile
def split_walls(r, min_sep):
    """Split a band's radii into inner/outer wall by a 1-D Otsu threshold.
    Returns (r_inner, r_outer, two_walled)."""
    if len(r) < 12:
        m = float(np.median(r)); return m, m, False
    lo, hi = np.percentile(r, [2, 98])
    if hi - lo < min_sep:
        m = float(np.median(r)); return m, m, False
    bins = np.linspace(lo, hi, 33)
    hist, _ = np.histogram(r, bins=bins)
    centers = 0.5 * (bins[:-1] + bins[1:])
    w0 = np.cumsum(hist); w1 = hist.sum() - w0
    mu = np.cumsum(hist * centers)
    with np.errstate(divide="ignore", invalid="ignore"):
        between = np.nan_to_num((mu[-1] * w0 - mu) ** 2 / (w0 * w1))
    t = centers[int(np.argmax(between[:-1]))]
    inner, outer = r[r <= t], r[r > t]
    if len(inner) < 3 or len(outer) < 3:
        m = float(np.median(r)); return m, m, False
    ri, ro = float(np.median(inner)), float(np.median(outer))
    if ro - ri < min_sep:
        m = float(np.median(r)); return m, m, False
    return ri, ro, True


def band_profile(h, r, ang, n_bands=300, min_sep=0.4):
    """Per height band: centre h, median r, 98th-pct r, angular coverage, and the
    Otsu inner/outer split. Returns a dict of equal-length arrays."""
    edges = np.linspace(h.min(), h.max(), n_bands + 1)
    out = {k: [] for k in ("h", "r_med", "r_max", "cover", "r_in", "r_out", "two")}
    which = np.clip(np.searchsorted(edges, h, side="right") - 1, 0, n_bands - 1)
    order = np.argsort(which)
    bounds = np.searchsorted(which[order], np.arange(n_bands + 1))
    for i in range(n_bands):
        sel = order[bounds[i]:bounds[i + 1]]
        if len(sel) < 10:
            continue
        rb = r[sel]
        filled = len(np.unique(np.floor((ang[sel] + np.pi) / (2 * np.pi) * 24).astype(int)))
        ri, ro, two = split_walls(rb, min_sep)
        out["h"].append(0.5 * (edges[i] + edges[i + 1]))
        out["r_med"].append(float(np.median(rb)))
        out["r_max"].append(float(np.percentile(rb, 98)))
        out["cover"].append(filled / 24.0)
        out["r_in"].append(ri); out["r_out"].append(ro); out["two"].append(two)
    return {k: np.asarray(v) for k, v in out.items()}


# ---------------------------------------------------------------- plane section
def _plane_segments(V, T, p0, d, phi):
    """Intersect the mesh with the half-plane through the axis at angle phi.
    Returns (segments as pairs of edge keys, dict edge_key -> (r, h))."""
    e1, e2, d = basis_from_dir(d)
    u = np.cos(phi) * e1 + np.sin(phi) * e2          # in-plane radial direction
    n = -np.sin(phi) * e1 + np.cos(phi) * e2         # plane normal
    W = V - p0
    s = W @ n
    s[s == 0] = 1e-12
    sg = s > 0
    tri_s = sg[T]
    cross = tri_s.any(1) & ~tri_s.all(1)
    TT = T[cross]
    pts = {}
    segs = []
    for tri in TT:
        keys = []
        for a, b in ((tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])):
            if sg[a] != sg[b]:
                k = (a, b) if a < b else (b, a)
                if k not in pts:
                    t = s[a] / (s[a] - s[b])
                    P = W[a] + t * (W[b] - W[a])
                    pts[k] = (float(P @ u), float(P @ d))
                keys.append(k)
        if len(keys) == 2 and (pts[keys[0]][0] > 0 or pts[keys[1]][0] > 0):
            segs.append(tuple(keys))
    return segs, pts


def _chain(segs):
    """Join segments sharing endpoints into polylines. Returns list of (keys, closed)."""
    from collections import defaultdict
    adj = defaultdict(list)
    for i, (a, b) in enumerate(segs):
        adj[a].append(i); adj[b].append(i)
    used = np.zeros(len(segs), bool)
    chains = []

    def walk(start_key, seg_i):
        keys = [start_key]
        cur_key, cur = start_key, seg_i
        while True:
            used[cur] = True
            a, b = segs[cur]
            nxt_key = b if a == cur_key else a
            keys.append(nxt_key)
            cand = [j for j in adj[nxt_key] if not used[j]]
            if not cand:
                return keys
            cur_key, cur = nxt_key, cand[0]

    # open chains first (start at degree-1 endpoints), then closed loops
    for k, lst in adj.items():
        if len(lst) == 1 and not used[lst[0]]:
            chains.append((walk(k, lst[0]), False))
    for i in range(len(segs)):
        if not used[i]:
            keys = walk(segs[i][0], i)
            chains.append((keys, keys[0] == keys[-1]))
    return chains


def plane_section(V, T, p0, d, phi, axis_tol=1.0):
    """Best polyline of the half-plane cut at angle phi, in (r, h) mesh units.
    Returns (pts (K,2), kind) with kind 'section' (closed, or both ends on the axis:
    a two-walled wall) or 'line' (open: a single surface) or (None, None)."""
    segs, pts = _plane_segments(V, T, p0, d, phi)
    if not segs:
        return None, None
    best, best_ext = None, -1
    for keys, closed in _chain(segs):
        P = np.array([pts[k] for k in keys])
        if len(P) < 4:
            continue
        ext = np.ptp(P[:, 1]) + 0.25 * np.ptp(P[:, 0])
        if ext > best_ext:
            best, best_ext, kind = P, ext, ("section" if closed else "line")
    if best is None:
        return None, None
    if kind == "line" and best[0, 0] < axis_tol and best[-1, 0] < axis_tol and np.ptp(best[:, 0]) > 3 * axis_tol:
        kind = "section"           # complete vessel: the cut runs from axis to axis
    best[:, 0] = np.maximum(best[:, 0], 0.0)
    return best, kind


def best_section(V, T, p0, d, h, ang, exterior, n_candidates=15):
    """Cut at several angles across the sherd and keep the most complete profile
    (largest height extent), preferring among near-complete cuts the one that best
    matches the robust exterior (so a cut through a handle or a chip loses)."""
    # angular span actually present: circular percentiles around the mean direction
    c = np.angle(np.mean(np.exp(1j * ang)))
    rel = (ang - c + np.pi) % (2 * np.pi) - np.pi
    lo, hi = np.percentile(rel, [10, 90])
    cands = []
    for phi in c + np.linspace(lo, hi, n_candidates):
        P, kind = plane_section(V, T, p0, d, phi)
        if P is None:
            continue
        cands.append((phi, P, kind))
    if not cands:
        return None, None, None, []
    ext = np.array([np.ptp(P[:, 1]) for _, P, _ in cands])
    good = [cands[i] for i in np.where(ext >= 0.95 * ext.max())[0]]
    eh, er = exterior

    def mismatch(P):
        # compare the outermost cut radius per height with the robust exterior
        hb = np.round(P[:, 1], 0)
        o = [P[hb == v, 0].max() for v in np.unique(hb)]
        return float(np.median(np.abs(np.array(o) - np.interp(np.unique(hb), eh, er))))

    phi, P, kind = min(good, key=lambda g: mismatch(g[1]))
    return phi, P, kind, good


# ---------------------------------------------------------------- walls + thickness
def split_section(P):
    """Split a section polygon into the outer and inner wall at the lip (top-most
    point) and the lowest point. Returns (outer, inner), each ordered top -> bottom."""
    top = int(np.argmax(P[:, 1]))
    Q = np.roll(P, -top, axis=0)                    # starts at the lip
    bot = int(np.argmin(Q[:, 1]))
    a, b = Q[:bot + 1], np.vstack([Q[bot:], Q[:1]])[::-1]
    return (a, b) if np.median(a[:, 0]) >= np.median(b[:, 0]) else (b, a)


def _point_to_polyline(X, L):
    """Min distance from each point in X (n,2) to polyline L (m,2)."""
    A, B = L[:-1], L[1:]
    AB = B - A
    L2 = np.maximum((AB ** 2).sum(1), 1e-18)
    t = np.clip(((X[:, None, :] - A[None]) * AB[None]).sum(2) / L2[None], 0, 1)
    proj = A[None] + t[..., None] * AB[None]
    return np.sqrt(((X[:, None, :] - proj) ** 2).sum(2)).min(1)


def wall_thickness(outer, inner, interval, top_h, skip_ends=2.0):
    """Perpendicular wall thickness (calliper-style: shortest distance from the outer
    wall to the inner wall), sampled every `interval` along the height below the rim.
    Sampling stops within one wall thickness (at least `skip_ends`) of the lowest point,
    because there the shortest path runs across the break face, not through the wall.
    Returns list of (depth, thickness)."""
    bottom = max(outer[:, 1].min(), inner[:, 1].min())
    out = []
    depth, last_t = interval, 0.0
    while True:
        hq = top_h - depth
        if hq - bottom < max(skip_ends, 1.2 * last_t):
            break
        o = outer
        k = np.where(np.diff(np.sign(o[:, 1] - hq)) != 0)[0]
        if len(k):
            k = k[0]
            t = (hq - o[k, 1]) / (o[k + 1, 1] - o[k, 1])
            x = o[k] + t * (o[k + 1] - o[k])
            last_t = float(_point_to_polyline(x[None], inner)[0])
            out.append((round(float(depth), 3), round(last_t, 3)))
        depth += interval
    return out


def thickness_over_cuts(cuts, interval):
    """Median perpendicular thickness per depth over several section cuts (each cut is
    noisy on its own; the median of near-complete cuts is not). Returns [(depth, t)]."""
    per = {}
    for _, P, kind in cuts:
        if kind != "section":
            continue
        o, i = split_section(P)
        for dep, t in wall_thickness(o, i, interval, P[:, 1].max()):
            per.setdefault(dep, []).append(t)
    return [(dep, round(float(np.median(v)), 3)) for dep, v in sorted(per.items()) if len(v) >= 1]


def polyline_rms_distance(A, B):
    """Symmetric RMS of point-to-polyline distances between two polylines/polygons."""
    da = _point_to_polyline(A, np.vstack([B, B[:1]]))
    db = _point_to_polyline(B, np.vstack([A, A[:1]]))
    return float(np.sqrt(np.mean(np.concatenate([da, db]) ** 2))), float(max(da.max(), db.max()))
