"""The whole CPU lane in one call: mesh -> cleaned mesh -> axis -> profile -> drawing,
with a GO / REVIEW verdict. Used by the CLI, the desktop app and the tests."""
import copy
import json
import os
import time

import numpy as np

from . import __version__
from .meshio import load_mesh, save_mesh
from .cleanup import clean, close_holes
from .axis import estimate_axis, fit_axis_manual
from . import profile as pf
from . import drawing as dr

DEFAULTS = {
    "sherd": {"id": None},                      # default: input file name
    "scale": {"mm_per_unit": 1.0,               # mesh units -> mm
              "mesh_distance": None,            # or: a distance measured on the mesh...
              "distance_mm": None},             # ...and the same distance by calliper
    "cleanup": {"keep_largest_only": True, "min_faces": None,
                "close_holes_max_edges": 0, "save_mesh": True},
    "profile": {"axis_method": "auto",          # auto | pca | rim_arc | normals_refine | manual
                "axis": {"direction": [0, 1, 0], "point": [0, 0, 0]},   # used by manual
                "orientation": "auto",          # auto | flip
                "bands": 300, "thickness_interval_mm": 5.0, "min_wall_sep_mm": 0.4},
    "drawing": {"section_side": "left", "fill": "solid", "line_weight_mm": 0.3,
                "scale_bar": True, "caption": True, "formats": ["svg", "png", "pdf"]},
    "ware": "plain",                            # plain | decorated
}
CHOICES = {
    ("profile", "axis_method"): {"auto", "pca", "rim_arc", "normals_refine", "manual"},
    ("profile", "orientation"): {"auto", "flip"},
    ("drawing", "section_side"): {"left", "right"},
    ("drawing", "fill"): {"solid", "hatch", "none"},
    ("ware",): {"plain", "decorated"},
}
# Gates. Provisional: set from synthetic tests (cpu/tools/axis_methods_eval.py);
# recheck once real sherds have been measured.
AMBIGUITY_MM_MAX = 0.3
COVERAGE_MIN = 0.25
COVERAGE_MIN_SURFACE_FIT = 0.10
KEPT_FRACTION_MIN = 0.5


def load_config(path_or_dict=None):
    """Merge a YAML file (or dict) over DEFAULTS. Returns (config, warnings)."""
    user = {}
    if isinstance(path_or_dict, dict):
        user = path_or_dict
    elif path_or_dict:
        import yaml
        with open(path_or_dict, encoding="utf-8") as fh:
            user = yaml.safe_load(fh) or {}
    cfg = copy.deepcopy(DEFAULTS)
    warnings = []
    for k, v in user.items():
        if k not in cfg:
            warnings.append(f"unknown config key '{k}' ignored")
        elif isinstance(cfg[k], dict) and isinstance(v, dict):
            for k2, v2 in v.items():
                if k2 not in cfg[k]:
                    warnings.append(f"unknown config key '{k}.{k2}' ignored")
                else:
                    cfg[k][k2] = v2
        else:
            cfg[k] = v
    for path, allowed in CHOICES.items():
        val = cfg[path[0]] if len(path) == 1 else cfg[path[0]][path[1]]
        if val not in allowed:
            raise ValueError(f"config {'.'.join(path)} = {val!r}; expected one of {sorted(allowed)}")
    return cfg, warnings


def mm_per_unit(cfg):
    s = cfg["scale"]
    if s.get("mesh_distance") and s.get("distance_mm"):
        return float(s["distance_mm"]) / float(s["mesh_distance"])
    return float(s.get("mm_per_unit") or 1.0)


def run(input_path, outdir, config=None, log=print):
    """Run the lane. Returns the report dict (also written to <id>_report.json)."""
    t0 = time.time()
    cfg, warnings = load_config(config)
    for w in warnings:
        log(f"[config] {w}")
    os.makedirs(outdir, exist_ok=True)
    sid = cfg["sherd"]["id"] or os.path.splitext(os.path.basename(input_path))[0]
    base = os.path.join(outdir, _safe(sid))

    log(f"Loading {os.path.basename(input_path)}")
    V, T = load_mesh(input_path)
    k = mm_per_unit(cfg)
    V = V * k

    log("Stage 2: cleanup")
    if cfg["cleanup"]["keep_largest_only"] or cfg["cleanup"]["min_faces"]:
        V, T, clean_rep = clean(V, T, cfg["cleanup"]["min_faces"])
    else:
        clean_rep = {"skipped": True}
    if cfg["cleanup"]["close_holes_max_edges"]:
        V, T, msg = close_holes(V, T, cfg["cleanup"]["close_holes_max_edges"])
        clean_rep["close_holes"] = msg
    if cfg["cleanup"]["save_mesh"]:
        save_mesh(base + "_clean.ply", V, T)

    log("Stage 3: axis")
    pcfg = cfg["profile"]
    if pcfg["axis_method"] == "manual":
        p0, d, axis_info = fit_axis_manual(**pcfg["axis"])
        p0 = p0 * k
    else:
        p0, d, axis_info = estimate_axis(V, T, pcfg["axis_method"])
    d, orient_how = pf.orient_axis(V, T, p0, d, pcfg["orientation"])

    log("Stage 3: profile and section")
    h, r, ang, _ = pf.frame(V, p0, d)
    bp = pf.band_profile(h, r, ang, pcfg["bands"], pcfg["min_wall_sep_mm"])
    ext_r = np.where(bp["two"], bp["r_out"], bp["r_med"])
    phi, P, kind, cuts = pf.best_section(V, T, p0, d, h, ang, (bp["h"], ext_r))
    if P is None:                                   # no usable cut: fall back to the band line
        P, kind, cuts = np.column_stack([ext_r, bp["h"]]), "line", []
    thickness = pf.thickness_over_cuts(cuts, pcfg["thickness_interval_mm"]) if kind == "section" else []

    log("Stage 4: drawing")
    if kind == "section":       # exterior view mirrors the drawn section's outer wall
        o, _ = pf.split_section(P)
        view_h, view_r = o[:, 1], o[:, 0]
    else:
        view_h, view_r = bp["h"], ext_r
    G = dr.layout(P, kind, view_h, view_r, cfg["drawing"]["section_side"])
    th_vals = [t for _, t in thickness]
    caption = ""
    if cfg["drawing"]["caption"]:
        bits = [f"rim Ø {2 * G['r_rim']:.0f} mm"]
        if th_vals:
            bits.append(f"wall {min(th_vals):.1f}–{max(th_vals):.1f} mm")
        bits.append(f"axis {axis_info.get('method')}")
        bits.append(f"archaemetry {__version__}")
        caption = " · ".join(bits)
    outputs = {}
    D = cfg["drawing"]
    if "svg" in D["formats"]:
        outputs["svg"] = base + "_drawing.svg"
        dr.to_svg(G, outputs["svg"], sid, caption, D["fill"], D["line_weight_mm"], D["scale_bar"])
    raster = [f for f in ("png", "pdf") if f in D["formats"]]
    if raster:
        for f in raster:
            outputs[f] = base + f"_drawing.{f}"
        dr.to_raster(G, [outputs[f] for f in raster], sid, caption, D["fill"], D["line_weight_mm"], D["scale_bar"])
    outputs["qa"] = base + "_qa.png"
    dr.qa_figure(outputs["qa"], h, r, {"h": bp["h"], "r_out": ext_r}, P, thickness, title=sid)

    # ---- verdict
    method = axis_info.get("method", "")
    surface_fit = method == "normals_multistart"
    coverage = float(np.mean(bp["cover"]))
    flags, notes = [], []
    if "fallback" in method:
        flags.append("rim_arc could not fix the axis from a clean rim arc (too little or ragged rim)")
    amb = axis_info.get("ambiguity_mm")
    if amb is not None and amb > AMBIGUITY_MM_MAX:
        flags.append(f"axis ambiguous: equally good axes disagree by {amb} mm (sphere/cone-like or too little surface)")
    if axis_info.get("wall_classification_ok") is False:
        flags.append("could not separate outer and inner wall surfaces; axis unreliable")
    cov_min = COVERAGE_MIN_SURFACE_FIT if surface_fit else COVERAGE_MIN
    if coverage < cov_min:
        flags.append(f"very low angular coverage ({coverage:.2f})")
    kf = clean_rep.get("kept_fraction")
    if kf is not None and kf < KEPT_FRACTION_MIN:
        flags.append(f"cleanup kept only {kf:.0%} of the mesh; check the sherd is intact")
    if cfg["ware"] == "decorated":
        flags.append("decorated ware: rollout and decoration need the GigaMesh station")
    if surface_fit:
        notes.append("axis from the experimental surface-fit method; check the drawing")
    if kind != "section":
        notes.append("only one surface in the mesh: drawn as an outline, no wall thickness")
    if orient_how.startswith("world_up"):
        notes.append("rim direction guessed from mesh orientation; check the rim is at the top (use orientation: flip)")
    verdict = "REVIEW" if flags else "GO"

    profile_data = {
        "units": "mm", "h_origin": "rim at top; h increases toward the rim",
        "section_kind": kind, "section_rh": np.round(P, 4).tolist(),
        "exterior_rh": np.round(np.column_stack([ext_r, bp["h"]]), 4).tolist(),
        "rim_radius_mm": round(G["r_rim"], 3),
        "thickness_mm": [{"depth_below_rim": dpt, "thickness": t} for dpt, t in thickness],
    }
    with open(base + "_profile.json", "w") as fh:
        json.dump(profile_data, fh)
    outputs["profile"] = base + "_profile.json"

    report = {
        "archaemetry_version": __version__,
        "sherd_id": sid, "input": os.path.basename(input_path),
        "verdict": verdict, "flags": flags, "notes": notes, "config_warnings": warnings,
        "scale_mm_per_unit": k,
        "cleanup": clean_rep,
        "axis": {"method": method, "requested": pcfg["axis_method"], "orientation": orient_how,
                 "direction": [round(float(x), 5) for x in d],
                 "point_mm": [round(float(x), 3) for x in p0],
                 **{kk: v for kk, v in axis_info.items() if kk != "method"}},
        "profile": {"section_kind": kind, "height_mm": round(float(np.ptp(bp["h"])), 3),
                    "rim_diameter_mm": round(2 * G["r_rim"], 2),
                    "max_radius_mm": round(float(ext_r.max()), 3),
                    "mean_angular_coverage": round(coverage, 3),
                    "cut_angle_deg": None if phi is None else round(float(np.degrees(phi)), 2),
                    "n_cuts_used_for_thickness": sum(1 for c in cuts if c[2] == "section")},
        "wall_thickness_mm": ({"min": min(th_vals), "max": max(th_vals),
                               "mean": round(float(np.mean(th_vals)), 3),
                               "interval_mm": pcfg["thickness_interval_mm"]} if th_vals else None),
        "outputs": {kk: os.path.basename(v) for kk, v in outputs.items()},
        "runtime_s": round(time.time() - t0, 2),
        "config": cfg,
    }
    with open(base + "_report.json", "w") as fh:
        json.dump(report, fh, indent=2, default=_json_default)
    log(f"Verdict: {verdict}" + "".join(f"\n  ! {f}" for f in flags) + "".join(f"\n  - {n}" for n in notes))
    return report


def _json_default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _safe(s):
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in str(s)) or "sherd"
