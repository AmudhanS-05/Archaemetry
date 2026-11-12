"""Self-test: run known geometry through the real pipeline and compare the drawn
section, exterior and wall thickness with the exact answer.

`archaemetry selftest` runs the quick set (about a minute); tools/synthetic_validation.py
runs the full set. Used by CI as the accuracy regression check."""
import json
import os
import tempfile

import numpy as np

from . import synthetic as sy
from .meshio import save_mesh
from .pipeline import run
from .profile import polyline_rms_distance, split_section, wall_thickness

# pass/fail limits for the quick set (mm). Chosen with headroom over measured values.
LIMITS = {"section_rms": 0.35, "thickness_mean_err": 0.35, "exterior_rms": 0.35}


def _align_top(P, top):
    Q = P.copy(); Q[:, 1] += top - Q[:, 1].max(); return Q


def _case_closed(name, V, T, loop, cfg, outdir, interval=5.0):
    path = os.path.join(outdir, name + ".ply")
    save_mesh(path, V, T)
    cfg = dict(cfg); cfg.setdefault("sherd", {"id": name})
    cfg.setdefault("drawing", {"formats": ["svg", "png"]})
    rep = run(path, os.path.join(outdir, name), cfg, log=lambda *_: None)
    prof = json.load(open(os.path.join(outdir, name, rep["outputs"]["profile"])))
    P = _align_top(np.array(prof["section_rh"]), loop[:, 1].max())
    res = {"case": name, "verdict": rep["verdict"], "flags": rep["flags"],
           "axis_method": rep["axis"]["method"], "section_kind": prof["section_kind"]}
    if prof["section_kind"] == "section":
        rms, mx = polyline_rms_distance(P, loop)
        res.update(section_rms=round(rms, 3), section_max=round(mx, 3))
        o, i = split_section(loop)
        truth = dict(wall_thickness(o, i, interval, loop[:, 1].max()))
        got = {d["depth_below_rim"]: d["thickness"] for d in prof["thickness_mm"]}
        errs = [abs(got[k] - truth[k]) for k in got if k in truth]
        if errs:
            res.update(thickness_mean_err=round(float(np.mean(errs)), 3),
                       thickness_max_err=round(float(np.max(errs)), 3), thickness_n=len(errs))
    res["drawing"] = os.path.join(outdir, name, rep["outputs"].get("png", ""))
    return res


def _case_exterior(name, V, T, cfg, outdir, knots=None):
    path = os.path.join(outdir, name + ".ply")
    save_mesh(path, V, T)
    cfg = dict(cfg); cfg.setdefault("sherd", {"id": name})
    cfg.setdefault("drawing", {"formats": ["svg"]})
    rep = run(path, os.path.join(outdir, name), cfg, log=lambda *_: None)
    prof = json.load(open(os.path.join(outdir, name, rep["outputs"]["profile"])))
    E = np.array(prof["exterior_rh"])
    kh = sy.KH if knots is None else knots[0]
    E[:, 1] += (V[:, 1].max() if rep["axis"]["direction"][1] > 0 else -V[:, 1].min()) - E[:, 1].max()
    if rep["axis"]["direction"][1] < 0:
        E[:, 1] = -E[:, 1]
    rt = sy.true_r(E[:, 1], knots)
    inside = (E[:, 1] > kh[0] + 1) & (E[:, 1] < kh[-1] - 1)
    e = E[inside, 0] - rt[inside]
    return {"case": name, "verdict": rep["verdict"], "flags": rep["flags"],
            "axis_method": rep["axis"]["method"], "section_kind": prof["section_kind"],
            "exterior_rms": round(float(np.sqrt(np.mean(e ** 2))), 3),
            "exterior_max": round(float(np.abs(e).max()), 3)}


def cases(outdir, full=False, noise=0.05):
    rng = np.random.default_rng(0)
    out = []
    V, T, L = sy.make_vessel(noise=noise, rng=rng)
    out.append(_case_closed("vessel_two_walled", V, T, L, {"profile": {"axis_method": "pca"}}, outdir))
    for th in ((120, 90, 60, 45) if full else (120,)):
        V, T, L = sy.make_sherd(th, noise=noise, rng=rng)
        out.append(_case_closed(f"rim_sherd_{th}deg", V, T, L, {"profile": {"axis_method": "auto"}}, outdir))
    V, T = sy.revolve(rng=rng, noise=noise)
    out.append(_case_exterior("vessel_exterior_only", V, T, {"profile": {"axis_method": "pca"}}, outdir))
    return out


def check(results):
    """Return list of failure messages for the quick set (empty = pass). Cases flagged
    REVIEW are expected to be inaccurate and are judged only on being flagged."""
    fails = []
    for r in results:
        if r["verdict"] != "GO":
            continue
        for k, lim in LIMITS.items():
            if k in r and r[k] > lim:
                fails.append(f"{r['case']}: {k} = {r[k]} mm > {lim} mm")
    return fails


def main(outdir=None, full=False, log=print):
    outdir = outdir or tempfile.mkdtemp(prefix="archaemetry_selftest_")
    os.makedirs(outdir, exist_ok=True)
    log(f"Self-test ({'full' if full else 'quick'}) in {outdir}")
    res = cases(outdir, full=full)
    for r in res:
        nums = ", ".join(f"{k} {r[k]}" for k in ("section_rms", "thickness_mean_err", "exterior_rms") if k in r)
        log(f"  {r['case']:<24} {r['verdict']:<6} {nums}  [{r['axis_method']}]")
    fails = check(res)
    json.dump(res, open(os.path.join(outdir, "selftest.json"), "w"), indent=2)
    log("PASS" if not fails else "FAIL\n  " + "\n  ".join(fails))
    return not fails, res
