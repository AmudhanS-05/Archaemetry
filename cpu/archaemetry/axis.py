"""Symmetry-axis estimators.

pca            longest extent as a seed, refined through circle centres of slices
               (complete or large vessels)
rim_arc        plane + circle fitted to the rim edge; refuses short/ragged arcs
normals_refine EXPERIMENTAL: surface-normal line fit + profile-consistency refinement,
               multi-start with an ambiguity check (small rim and body sherds)
manual         a known axis (e.g. the mesh was oriented in MeshLab with axis = Y)

All return (point_on_axis, unit_direction, info_dict).
"""
import numpy as np


# ---------- geometry helpers ----------
def robust_circle_fit(P2, n_iter=3, k=2.5):
    """Algebraic (Kasa) circle fit with iterative outlier rejection.
    P2: (N,2). Returns center(2,), radius, inlier_fraction."""
    idx = np.ones(len(P2), dtype=bool)
    cx = cy = r = 0.0
    for _ in range(n_iter):
        p = P2[idx]
        if len(p) < 8:
            break
        u, v = p[:, 0], p[:, 1]
        A = np.column_stack([u, v, np.ones_like(u)])
        b = -(u**2 + v**2)
        sol, *_ = np.linalg.lstsq(A, b, rcond=None)
        D, E, F = sol
        cx, cy = -D / 2, -E / 2
        r = np.sqrt(max(cx**2 + cy**2 - F, 1e-12))
        d = np.abs(np.sqrt((P2[:, 0] - cx)**2 + (P2[:, 1] - cy)**2) - r)
        mad = np.median(d) + 1e-9
        idx = d < k * mad
    return np.array([cx, cy]), r, float(idx.mean())


def basis_from_dir(d):
    """Orthonormal (e1,e2) spanning the plane perpendicular to unit vector d."""
    d = d / np.linalg.norm(d)
    a = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(d, a); e1 /= np.linalg.norm(e1)
    e2 = np.cross(d, e1)
    return e1, e2, d


def fit_axis(V, n_bands=40, n_refine=4, seed_dir=None):
    """PCA + circle-centre refinement. Without a seed, refine from each of the three
    principal directions and keep the most rotationally consistent result: the
    longest extent is the axis of a tall vessel but the diameter of a squat bowl."""
    if seed_dir is not None:
        return _fit_axis_seeded(V, n_bands, n_refine, seed_dir)
    c = V.mean(0)
    _, evecs = np.linalg.eigh(np.cov((V - c).T))
    best = None
    for k in (2, 1, 0):
        p, d, resid = _fit_axis_seeded(V, n_bands, n_refine, evecs[:, k])
        cost = consistency_cost(V[_subsample(len(V), MAX_FIT_POINTS)], p, d)
        if best is None or cost < best[0]:
            best = (cost, p, d, resid)
    return best[1], best[2], best[3]


def _fit_axis_seeded(V, n_bands=40, n_refine=4, seed_dir=None):
    """Return (point_on_axis, unit_dir, per-slice residual). Circle-centre
    refinement seeded by PCA major axis, or by `seed_dir` when given (used by
    rim_arc to avoid PCA's seed flipping on wider-than-tall rim bands)."""
    c = V.mean(0)
    if seed_dir is not None:
        d = np.asarray(seed_dir, float)
        d = d / np.linalg.norm(d)
    else:
        cov = np.cov((V - c).T)
        evals, evecs = np.linalg.eigh(cov)
        d = evecs[:, np.argmax(evals)]      # PCA seed = major axis
    p0 = c.copy()
    resid = np.nan
    for _ in range(n_refine):
        e1, e2, d = basis_from_dir(d)
        h = (V - p0) @ d
        lo, hi = np.percentile(h, [5, 95])   # trim caps
        edges = np.linspace(lo, hi, n_bands + 1)
        centers3d, resids = [], []
        for i in range(n_bands):
            m = (h >= edges[i]) & (h < edges[i + 1])
            if m.sum() < 30:
                continue
            Pw = V[m] - p0
            P2 = np.column_stack([Pw @ e1, Pw @ e2])
            ctr2, r, inl = robust_circle_fit(P2)
            hb = h[m].mean()
            c3 = p0 + ctr2[0] * e1 + ctr2[1] * e2 + hb * d
            centers3d.append(c3)
            resids.append(r)
        centers3d = np.array(centers3d)
        if len(centers3d) < 4:
            break
        # refit axis as best line through the circle centres
        cc = centers3d.mean(0)
        _, _, Vt = np.linalg.svd(centers3d - cc)
        d = Vt[0]
        p0 = cc
        # residual = how far the circle centres deviate from the fitted line
        proj = (centers3d - cc) @ d
        line_pts = cc + np.outer(proj, d)
        resid = float(np.sqrt(((centers3d - line_pts)**2).sum(1)).mean())
    return p0, d / np.linalg.norm(d), resid


def fit_planar_circle(P3):
    """Fit a plane + circle to 3D points. Return (center3d, normal, radius,
    norm_residual, arc_span_deg, planarity). planarity = out-of-plane extent /
    in-plane extent; small => a thin planar band (a clean rim), large => a
    non-planar band (e.g. a flared mouth) whose normal is unreliable."""
    c = P3.mean(0)
    _, S, Vt = np.linalg.svd(P3 - c)
    normal = Vt[2]                       # smallest-variance dir = plane normal
    e1, e2 = Vt[0], Vt[1]
    planarity = float(S[2] / (S[0] + 1e-12))
    P2 = np.column_stack([(P3 - c) @ e1, (P3 - c) @ e2])
    ctr2, r, inl = robust_circle_fit(P2)
    d = np.abs(np.sqrt((P2[:, 0] - ctr2[0])**2 + (P2[:, 1] - ctr2[1])**2) - r)
    norm_resid = float(np.median(d) / (r + 1e-9))
    ang = np.arctan2(P2[:, 1] - ctr2[1], P2[:, 0] - ctr2[0])
    arc_span = float(np.degrees(ang.max() - ang.min()))
    center3d = c + ctr2[0] * e1 + ctr2[1] * e2
    return center3d, normal, r, norm_resid, arc_span, planarity


def border_edges(T):
    """(M,2) vertex pairs of edges used by exactly one triangle (vectorised)."""
    E = np.sort(np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]]), axis=1)
    uniq, counts = np.unique(E, axis=0, return_counts=True)
    return uniq[counts == 1]


def border_vertices(T):
    """Indices of vertices on a border edge (edge used by exactly one triangle)."""
    return np.unique(border_edges(T).ravel())


def fit_axis_rim_arc(V, T, max_resid=0.05, min_span_deg=70.0, band_frac=0.15):
    """Axis from the horizontal circular arcs on the mesh border (the rim, and
    any clean horizontal break). A rim sherd is topologically a disk, so its
    border is one mixed loop; we isolate horizontal arcs by taking border
    vertices in several height bands (height along a PCA-seeded axis) and
    robustly circle-fitting each.

    The best arc's plane normal SEEDS a circle-centre refinement over the whole
    fragment (fit_axis with seed_dir) -- the rim normal is a robust direction for
    a long arc, and refinement then pins the position. We pick the arc by
    LARGEST span among those with a clean circle fit (a longer arc constrains the
    normal far better than a short one), and gate on a minimum span. Falls back
    to PCA when no clean, long-enough arc exists (a short-arc fragment cannot be
    auto-axised reliably -- it should be flagged for the manual/GigaMesh route).
    Mirrors the archaeological rim-chart method."""
    c = V.mean(0)
    evals, evecs = np.linalg.eigh(np.cov((V - c).T))
    d0 = evecs[:, np.argmax(evals)]

    bidx = border_vertices(T)
    if len(bidx) < 16:
        p0, d, resid = fit_axis(V)
        return p0, d, dict(method="pca_fallback", n_arcs=0, resid=resid)
    B = V[bidx]
    h = (B - c) @ d0
    top = B[h >= np.quantile(h, 1 - band_frac)]
    bot = B[h <= np.quantile(h, band_frac)]

    arcs = []
    for name, band in (("top", top), ("bot", bot)):
        if len(band) < 12:
            continue
        center, normal, r, nres, span, planarity = fit_planar_circle(band)
        # accept only a clean, planar AND long-enough arc: a short or non-planar
        # arc's plane normal is unreliable, so we would rather refuse (fall back)
        # than guess an axis.
        if (r > 0 and nres < max_resid and span >= min_span_deg
                and planarity < 0.12):
            arcs.append(dict(name=name, center=center, normal=normal, r=float(r),
                             nres=float(nres), span=float(span),
                             planarity=float(planarity)))
    if len(arcs) == 0:
        p0, d, resid = fit_axis(V)
        return p0, d, dict(method="pca_fallback", n_arcs=0, resid=resid,
                           note="no clean rim arc >= min_span_deg; flag for manual")

    # Axis = the best-fit (lowest-residual) clean arc's plane normal; position =
    # its centre. Circle-centre refinement was tried and made it worse (it tilts
    # on partial wedges), so we trust the arc normal directly.
    best = min(arcs, key=lambda a: a["nres"])
    d = best["normal"] / np.linalg.norm(best["normal"])
    p0 = best["center"]
    info = dict(method="rim_arc", n_arcs=len(arcs),
                rim_radius=round(best["r"], 3), rim_arc_deg=round(best["span"], 1),
                arc_resid=round(best["nres"], 4))
    return p0, d, info

def vertex_normals(V, T):
    f = np.cross(V[T[:, 1]] - V[T[:, 0]], V[T[:, 2]] - V[T[:, 0]])  # area-weighted
    N = np.zeros_like(V)
    for k in range(3):
        np.add.at(N, T[:, k], f)
    L = np.linalg.norm(N, axis=1)
    ok = L > 1e-12
    N[ok] /= L[ok, None]
    return N, ok


def _basis(d):
    a = np.array([1.0, 0, 0]) if abs(d[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(d, a); e1 /= np.linalg.norm(e1)
    return e1, np.cross(d, e1)


MAX_FIT_POINTS = 80_000          # a 4-6 parameter axis needs far fewer points than a scan has


def _subsample(n, k, seed=0):
    """Deterministic random subset of range(n) of size <= k (sorted)."""
    if n <= k:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, k, replace=False))


def fit_axis_normals(V, T, n_iter=12, ridge=1e-6, mask=None):
    """Return (p0, d, info). info carries conditioning/residual signals for gating.
    `mask` restricts the fit to a subset of vertices (e.g. wall surfaces only)."""
    N, ok = vertex_normals(V, T)
    if mask is not None:
        ok = ok & mask
    V, N = V[ok], N[ok]
    keep = _subsample(len(V), 3 * MAX_FIT_POINTS)
    V, N = V[keep], N[keep]
    c = V.mean(0)
    s = float(np.mean(np.linalg.norm(V - c, axis=1))) + 1e-9
    P = (V - c) / s
    Q = np.cross(P, N)                      # moment of each normal line
    A_rows = np.hstack([Q, N])              # residual = A_rows @ [d, m]
    w = np.ones(len(P))
    for it in range(n_iter):
        A = (A_rows * w[:, None]).T @ A_rows
        Add, Adm, Amm = A[:3, :3], A[:3, 3:], A[3:, 3:]
        Amm_r = Amm + ridge * np.trace(Amm) / 3 * np.eye(3)
        S = Add - Adm @ np.linalg.solve(Amm_r, Adm.T)        # Schur complement
        S = 0.5 * (S + S.T)
        ev, evec = np.linalg.eigh(S)
        d = evec[:, 0]
        m = -np.linalg.solve(Amm_r, Adm.T @ d)
        m = m - (d @ m) * d                                   # enforce d.m = 0
        res = A_rows @ np.concatenate([d, m])
        sig = 1.4826 * np.median(np.abs(res)) + 1e-9
        w = 1.0 / (1.0 + (res / (2.5 * sig)) ** 2)            # Cauchy IRLS
    p0 = c + s * np.cross(d, m)
    cond = float(ev[0] / (ev[1] + 1e-15))                     # ->0: axis well pinned; ->1: ambiguous
    info = dict(method="normals", n_used=int(len(P)),
                resid_mm=round(float(sig * s), 4), cond=round(cond, 4),
                inlier_frac=round(float((np.abs(res) < 3 * sig).mean()), 3))
    return p0, d / np.linalg.norm(d), info


def consistency_cost(V, p0, d, band_pts=60):
    """RMS within-band spread of radius (mm), equal-count height bands."""
    W = V - p0
    h = W @ d
    r = np.linalg.norm(W - np.outer(h, d), axis=1)
    order = np.argsort(h)
    r = r[order]
    nb = len(r) // band_pts
    if nb < 4:
        return np.inf
    B = r[:nb * band_pts].reshape(nb, band_pts)
    dev = B - np.median(B, axis=1, keepdims=True)
    return float(np.sqrt(np.mean(dev ** 2)))


def refine_consistency(V, p0, d, band_pts=60, tilt0=0.03, off0=1.5, n_round=60):
    """Pattern search over (2 tilt angles, 2 lateral offsets). Returns (p0, d, info)."""
    d = d / np.linalg.norm(d)
    e1, e2 = _basis(d)
    x = np.zeros(4)
    steps = np.array([tilt0, tilt0, off0, off0])

    def build(x):
        dd = d + x[0] * e1 + x[1] * e2
        dd /= np.linalg.norm(dd)
        return p0 + x[2] * e1 + x[3] * e2, dd

    best = consistency_cost(V, *build(x))
    start = best
    for _ in range(n_round):
        improved = False
        for i in range(4):
            for sgn in (+1, -1):
                y = x.copy(); y[i] += sgn * steps[i]
                cst = consistency_cost(V, *build(y))
                if cst < best - 1e-9:
                    x, best, improved = y, cst, True
        if not improved:
            steps *= 0.5
            if steps[0] < 1e-5:
                break
    p1, d1 = build(x)
    return p1, d1, dict(refined=True, spread_before_mm=round(start, 4),
                        spread_after_mm=round(best, 4),
                        shift_tilt_deg=round(float(np.degrees(np.hypot(x[0], x[1]))), 3),
                        shift_offset_mm=round(float(np.hypot(x[2], x[3])), 3))


def axis_uncertainty(V, p0, d, band_pts=60, n_mc=40, seed=0):
    """Indicative 1-sigma radial uncertainty (mm) of the drawn profile caused by axis
    error. Takes the curvature (Hessian) of the consistency cost around the fitted
    axis: a flat direction (e.g. a sphere cap, or a thin horizontal strip) means the
    axis is not identifiable and the uncertainty blows up. Covariance = 2*c0^2/N_eff
    * H^-1 with N_eff = number of height bands (neighbouring vertices are correlated,
    so this is deliberately conservative). Parameters match refine_consistency."""
    d = d / np.linalg.norm(d)
    e1, e2 = _basis(d)
    steps = np.array([0.01, 0.01, 0.5, 0.5])

    def f(x):
        dd = d + x[0] * e1 + x[1] * e2; dd /= np.linalg.norm(dd)
        return consistency_cost(V, p0 + x[2] * e1 + x[3] * e2, dd, band_pts) ** 2

    f0 = f(np.zeros(4))
    H = np.zeros((4, 4))
    for i in range(4):
        ei = np.zeros(4); ei[i] = steps[i]
        H[i, i] = (f(ei) - 2 * f0 + f(-ei)) / steps[i] ** 2
        for j in range(i + 1, 4):
            ej = np.zeros(4); ej[j] = steps[j]
            H[i, j] = H[j, i] = (f(ei + ej) - f(ei - ej) - f(-ei + ej) + f(-ei - ej)) / (4 * steps[i] * steps[j])
    n_eff = max(len(V) // band_pts, 4)
    ev, U = np.linalg.eigh(0.5 * (H + H.T))
    if ev.min() <= 1e-12 * max(ev.max(), 1e-30):
        return 99.0                                  # not identifiable
    cov = 2 * f0 / n_eff * (U @ np.diag(1 / ev) @ U.T)
    cov = 0.5 * (cov + cov.T)
    rng = np.random.default_rng(seed)
    L = np.linalg.cholesky(cov + 1e-18 * np.eye(4))
    r0 = np.linalg.norm((V - p0) - np.outer((V - p0) @ d, d), axis=1)
    acc = []
    for _ in range(n_mc):
        x = L @ rng.standard_normal(4)
        dd = d + x[0] * e1 + x[1] * e2; dd /= np.linalg.norm(dd)
        pp = p0 + x[2] * e1 + x[3] * e2
        W = V - pp
        r = np.linalg.norm(W - np.outer(W @ dd, dd), axis=1)
        acc.append(np.mean((r - r0) ** 2))
    return float(np.sqrt(np.mean(acc)))


def fit_axis_multistart(V, T, extra_seeds=(), n_jitter=4, jitter_deg=8.0, tol=1.15, seed=0):
    """Option 1 + 2 combined, with an ambiguity check.

    Seeds: surface-normal fit, PCA, any caller seeds, plus random tilts of the normals
    seed. Each is refined by profile consistency. The best-cost axis wins. Every
    refined axis whose cost is within `tol` x best is 'near-optimal'; if those
    disagree (a sphere cap, a cone, a too-small fragment: a continuum of equally good
    axes) the axis is ambiguous and the caller should flag the sherd for review.
    Returns (p0, d, info)."""
    rng = np.random.default_rng(seed)
    pn, dn, inf_n = fit_axis_normals(V, T)
    # A real sherd mesh is closed: outer wall, inner wall, and break faces. Break-face
    # normals do not meet the axis, so refit on wall vertices only (normal roughly
    # radial); and use only the OUTER wall for the consistency cost, otherwise the wall
    # thickness itself looks like radius spread.
    Nv, okn = vertex_normals(V, T)
    for _ in range(2):
        W = V - pn
        rad = W - np.outer(W @ dn, dn)
        rad /= np.maximum(np.linalg.norm(rad, axis=1), 1e-12)[:, None]
        cosr = (Nv * rad).sum(1)
        wall = okn & (np.abs(cosr) > 0.6)
        if wall.sum() < 200:
            break
        pn, dn, inf_n = fit_axis_normals(V, T, mask=wall)
    outer = _outer_wall(V, pn, dn, cosr, wall)
    Vfull, V = V, V[outer]
    V = V[_subsample(len(V), MAX_FIT_POINTS, seed)]
    c = Vfull.mean(0)
    cov = np.cov((V - c).T)
    ev, evec = np.linalg.eigh(cov)
    seeds = [(pn, dn), (c, evec[:, -1])] + list(extra_seeds)
    e1, e2 = _basis(dn)
    for _ in range(n_jitter):
        a = np.radians(jitter_deg) * rng.uniform(0.4, 1.0)
        ph = rng.uniform(0, 2 * np.pi)
        dj = dn + np.tan(a) * (np.cos(ph) * e1 + np.sin(ph) * e2)
        seeds.append((pn, dj / np.linalg.norm(dj)))
    sols = []
    for p, d in seeds:
        p1, d1, j = refine_consistency(V, p, d)
        sols.append((j["spread_after_mm"], p1, d1))
    sols.sort(key=lambda s: s[0])
    best_cost, pb, db = sols[0]
    near = [s for s in sols if s[0] <= best_cost * tol]
    r_best = _radii(V, pb, db)
    amb_deg = max(float(np.degrees(np.arccos(np.clip(abs(db @ s[2]), 0, 1)))) for s in near)
    amb_mm = max(float(np.sqrt(np.mean((_radii(V, s[1], s[2]) - r_best) ** 2))) for s in near)
    info = dict(method="normals_multistart", cond=inf_n["cond"], spread_mm=round(best_cost, 4),
                outer_wall_fraction=round(float(outer.mean()), 3),
                # a closed sherd is ~half outer wall, an open scan mostly outer wall;
                # far less means the walls could not be told apart from this seed
                wall_classification_ok=bool(outer.mean() >= 0.2),
                n_seeds=len(sols), n_near_optimal=len(near),
                ambiguity_deg=round(amb_deg, 2), ambiguity_mm=round(amb_mm, 3))
    return pb, db, info


def _outer_wall(V, p, d, cosr, wall):
    """Mask of outer-wall vertices. Wall vertices split by the sign of normal . radial;
    the group lying at larger radius is the outer wall (works whichever way the mesh is
    wound). An open, single-surface mesh has one dominant group: keep all wall points."""
    pos, neg = wall & (cosr > 0), wall & (cosr < 0)
    if min(pos.sum(), neg.sum()) < 0.15 * wall.sum():
        return wall if wall.sum() >= 200 else np.ones(len(V), bool)
    r = _radii(V, p, d)
    return pos if np.median(r[pos]) > np.median(r[neg]) else neg


def _radii(V, p, d):
    W = V - p
    return np.linalg.norm(W - np.outer(W @ d, d), axis=1)


def fit_axis_manual(direction=(0, 1, 0), point=(0, 0, 0)):
    d = np.asarray(direction, float)
    return np.asarray(point, float), d / np.linalg.norm(d), dict(method="manual")


def estimate_axis(V, T, method):
    """Dispatch by name. 'auto' tries rim_arc, then normals_refine if rim_arc refuses.
    info['method'] records what actually produced the axis."""
    if method == "pca":
        p, d, resid = fit_axis(V)
        return p, d, dict(method="pca", resid=round(float(resid), 4))
    if method == "rim_arc":
        return fit_axis_rim_arc(V, T)
    if method == "normals_refine":
        return fit_axis_multistart(V, T)
    if method == "auto":
        p, d, info = fit_axis_rim_arc(V, T)
        if info.get("method") == "rim_arc":
            return p, d, info
        p, d, info = fit_axis_multistart(V, T)
        info["note"] = "rim_arc refused; used normals_refine"
        return p, d, info
    raise ValueError(f"unknown axis method {method!r}")
