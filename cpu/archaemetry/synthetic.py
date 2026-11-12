"""Test geometry with an exactly known answer.

revolve()          exterior-only surface of revolution (the old validation vessel)
section_loop()     a 2D wall section (r, h): outer wall, rounded lip, inner wall, and
                   either a straight break (sherd) or a base reaching the axis (vessel)
sweep_section()    sweeps that section around the axis into a CLOSED, two-walled mesh;
                   for a sherd the two side break faces are capped, like a real break
"""
import numpy as np
import trimesh

# default exterior profile (h mm, r mm): flared foot, body, shoulder, neck, rim
KH = np.array([0, 4, 20, 50, 75, 85, 100, 110, 114.0])
KR = np.array([14, 17, 22, 21, 16, 7, 6, 9.5, 11.0])


def true_r(h, knots=None):
    kh, kr = (KH, KR) if knots is None else knots
    return np.interp(h, kh, kr)


def revolve(n_h=500, n_t=360, rng=None, noise=0.0, knots=None, oval=0.0):
    """Exterior-only surface of revolution about +Y. `oval` adds cos(2 theta)
    out-of-roundness (fraction of radius) to mimic hand-made ware."""
    kh, kr = (KH, KR) if knots is None else knots
    h = np.linspace(kh[0], kh[-1], n_h)
    t = np.linspace(0, 2 * np.pi, n_t, endpoint=False)
    H, TH = np.meshgrid(h, t, indexing="ij")
    R = np.interp(H, kh, kr) * (1 + oval * np.cos(2 * TH))
    V = np.stack([R * np.cos(TH), H, R * np.sin(TH)], -1).reshape(-1, 3)
    idx = np.arange(n_h * n_t).reshape(n_h, n_t)
    a, b = idx[:-1], idx[1:]
    a2, b2 = np.roll(a, -1, 1), np.roll(b, -1, 1)
    T = np.concatenate([np.stack([a, b, a2], -1).reshape(-1, 3),
                        np.stack([a2, b, b2], -1).reshape(-1, 3)])
    if noise:
        rng = rng or np.random.default_rng(0)
        n = rng.normal(0, noise, len(V))
        rad = np.stack([np.cos(TH), np.zeros_like(TH), np.sin(TH)], -1).reshape(-1, 3)
        V = V + rad * n[:, None]
    return V, T


def section_loop(outer_kh, outer_kr, thick_kh, thick_kt, h0, h1, n=160, n_lip=14,
                 base=None):
    """Closed section polygon in (r, h), counter-clockwise, no repeated end point.

    outer wall r_o(h) and horizontal wall thickness t(h) are piecewise-linear knots.
    The lip is a half-round joining outer and inner at h1. Below h0 the section ends in
    a straight break (sherd) unless `base` = floor thickness, in which case the outer
    wall turns in to the axis at h0 and the floor's top sits at h0 + base (vessel)."""
    hs = np.linspace(h0, h1, n)
    ro = np.interp(hs, outer_kh, outer_kr)
    ri = ro - np.interp(hs, thick_kh, thick_kt)
    outer = np.column_stack([ro, hs])
    inner = np.column_stack([ri, hs])[::-1]
    # half-round lip from outer top to inner top, bulging upward
    c = 0.5 * (outer[-1] + inner[0]); rad = 0.5 * (outer[-1][0] - inner[0][0])
    ang = np.linspace(0, np.pi, n_lip + 2)[1:-1]
    lip = np.column_stack([c[0] + rad * np.cos(ang), c[1] + rad * np.sin(ang)])
    if base is None:
        return np.vstack([outer, lip, inner])
    # vessel: inner wall stops at the floor, floor and underside run to the axis
    keep = inner[:, 1] >= h0 + base
    inner = inner[keep]
    floor = np.array([[0.0, h0 + base]])
    under = np.array([[0.0, h0]])
    return np.vstack([under, outer, lip, inner, floor])


def sweep_section(loop, theta_deg=360.0, n_t=240, noise=0.0, rng=None, cap=None):
    """Sweep a closed (r,h) loop about +Y. For theta < 360 the side faces are capped
    (closed, watertight solid with break faces). Returns (V, T) with consistent winding."""
    full = theta_deg >= 360.0
    t = np.radians(np.linspace(0, theta_deg, n_t, endpoint=not full))
    K = len(loop)
    R, Hh = loop[:, 0][None, :], loop[:, 1][None, :]
    C, S = np.cos(t)[:, None], np.sin(t)[:, None]
    V = np.stack([R * C, np.broadcast_to(Hh, (len(t), K)), R * S], -1).reshape(-1, 3)
    idx = np.arange(len(t) * K).reshape(len(t), K)
    a = idx[:, :]; b = np.roll(idx, -1, 1)
    if full:
        a2, b2 = np.roll(a, -1, 0), np.roll(b, -1, 0)
    else:
        a, b = a[:-1], b[:-1]
        a2, b2 = idx[1:], np.roll(idx, -1, 1)[1:]
    T = np.concatenate([np.stack([a, b, a2], -1).reshape(-1, 3),
                        np.stack([b, b2, a2], -1).reshape(-1, 3)])
    if not full:
        if cap is None:
            raise ValueError("a partial sweep needs cap triangles (sherd_cap_triangles)")
        for c0, c1 in ((idx[0][cap], idx[-1][cap][:, ::-1]), (idx[0][cap][:, ::-1], idx[-1][cap])):
            m = trimesh.Trimesh(V, np.concatenate([T, c0, c1]), process=True)
            if m.is_winding_consistent:
                break
    else:
        m = trimesh.Trimesh(V, T, process=True)   # merges axis-duplicate vertices
    if m.volume < 0:
        m.invert()
    V, T = np.asarray(m.vertices, float), np.asarray(m.faces, np.int64)
    if noise:
        rng = rng or np.random.default_rng(0)
        V = V + m.vertex_normals * rng.normal(0, noise, len(V))[:, None]
    return V, T


def sherd_cap_triangles(n, n_lip):
    """Triangulate a sherd section_loop (outer n pts, lip n_lip pts, inner n pts):
    quads between matching outer/inner heights, plus a fan over the convex lip."""
    o = np.arange(n)
    lip = n + np.arange(n_lip)
    inn = n + n_lip + np.arange(n)[::-1]          # inn[k] is at the height of o[k]
    tris = []
    for k in range(n - 1):
        tris += [[o[k], o[k + 1], inn[k + 1]], [o[k], inn[k + 1], inn[k]]]
    fan = [o[-1], *lip, inn[-1]]
    for j in range(1, len(fan) - 1):
        tris.append([fan[0], fan[j], fan[j + 1]])
    return np.array(tris, np.int64)


# a ready-made jar-like sherd with varying thickness, used by tests and the self-test
JAR_OUTER = (np.array([0, 30, 70, 95, 110.0]), np.array([60, 72, 70, 58, 64.0]))
JAR_THICK = (np.array([0, 50, 100, 110.0]), np.array([9.0, 7.0, 5.0, 4.0]))


def make_sherd(theta_deg=90.0, h0=40.0, h1=110.0, noise=0.0, rng=None, n=160, n_lip=14,
               outer=JAR_OUTER, thick=JAR_THICK, n_t=None):
    """Closed rim sherd (lip at h1, break at h0, two side breaks). Returns (V, T, loop)."""
    loop = section_loop(*outer, *thick, h0, h1, n=n, n_lip=n_lip)
    n_t = n_t or max(int(theta_deg * 1.5), 24)
    V, T = sweep_section(loop, theta_deg, n_t=n_t, noise=noise, rng=rng,
                         cap=sherd_cap_triangles(n, n_lip))
    return V, T, loop


def make_vessel(base=8.0, noise=0.0, rng=None, outer=JAR_OUTER, thick=JAR_THICK, n_t=240):
    """Complete two-walled vessel with a solid base reaching the axis. (V, T, loop)."""
    h1 = outer[0][-1]
    loop = section_loop(*outer, *thick, 0.0, h1, base=base)
    V, T = sweep_section(loop, 360.0, n_t=n_t, noise=noise, rng=rng)
    return V, T, loop
