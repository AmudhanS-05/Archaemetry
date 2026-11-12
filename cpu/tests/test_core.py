import json
import re

import numpy as np
import pytest

from archaemetry import synthetic as sy
from archaemetry import profile as pf
from archaemetry.axis import fit_axis, border_edges
from archaemetry.cleanup import clean
from archaemetry.meshio import load_mesh, save_mesh
from archaemetry.pipeline import load_config, mm_per_unit, run

Y = np.array([0, 1.0, 0])


def angle(a, b):
    return np.degrees(np.arccos(min(1.0, abs(a @ b) / np.linalg.norm(a) / np.linalg.norm(b))))


@pytest.mark.parametrize("ext", ["ply", "obj", "stl"])
def test_mesh_roundtrip(tmp_path, ext):
    V, T, _ = sy.make_sherd(60, n=40, n_lip=6)
    p = tmp_path / f"m.{ext}"
    save_mesh(p, V, T)
    V2, T2 = load_mesh(p)
    assert len(T2) == len(T)
    assert np.allclose(np.sort(V2, 0)[[0, -1]], np.sort(V, 0)[[0, -1]], atol=1e-4)


def test_cleanup_keeps_sherd_drops_turntable_and_noise():
    V, T, _ = sy.make_sherd(90, n=60)
    ang = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    disc = np.vstack([[0, -30, 0], np.c_[150 * np.cos(ang), -30 + 0 * ang, 150 * np.sin(ang)]])
    Td = np.array([[0, 1 + i, 1 + (i + 1) % 24] for i in range(24)]) + len(V)
    noise_V = np.random.default_rng(0).normal(0, 50, (30, 3))
    Tn = np.arange(30).reshape(10, 3) + len(V) + len(disc)
    V2, T2, rep = clean(np.vstack([V, disc, noise_V]), np.vstack([T, Td, Tn]))
    assert rep["components_before"] == 12 and len(T2) == len(T)


def test_pca_axis_tall_and_squat_vessels():
    V, T = sy.revolve(noise=0.05)                       # tall lekythos
    assert angle(fit_axis(V)[1], Y) < 0.5
    V, T, _ = sy.make_vessel()                          # squat jar, wider than tall
    assert angle(fit_axis(V)[1], Y) < 0.5


def test_section_and_thickness_exact_with_true_axis():
    V, T, loop = sy.make_sherd(90)
    assert len(border_edges(T)) == 0                    # closed, like a real scanned sherd
    h, r, ang, _ = pf.frame(V, np.zeros(3), Y)
    bp = pf.band_profile(h, r, ang)
    phi, P, kind, cuts = pf.best_section(V, T, np.zeros(3), Y, h, ang, (bp["h"], bp["r_out"]))
    assert kind == "section"
    rms, _ = pf.polyline_rms_distance(P, loop)
    assert rms < 0.01
    o, i = pf.split_section(loop)
    truth = dict(pf.wall_thickness(o, i, 5, loop[:, 1].max()))
    got = dict(pf.thickness_over_cuts(cuts, 5))
    errs = [abs(got[k] - truth[k]) for k in got if k in truth]
    assert len(errs) >= 10 and np.mean(errs) < 0.05


def test_config_rejects_bad_choice_and_warns_unknown():
    with pytest.raises(ValueError):
        load_config({"profile": {"axis_method": "magic"}})
    _, w = load_config({"profile": {"axis_metod": "pca"}})
    assert w and "axis_metod" in w[0]


def test_scale_from_marker_distance():
    cfg, _ = load_config({"scale": {"mesh_distance": 2.0, "distance_mm": 62.4}})
    assert mm_per_unit(cfg) == pytest.approx(31.2)


@pytest.fixture(scope="module")
def sherd_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("run")
    V, T, loop = sy.make_sherd(120, noise=0.05, rng=np.random.default_rng(3))
    save_mesh(d / "s.ply", V / 10.0, T)                  # mesh in cm: exercise scaling
    rep = run(str(d / "s.ply"), str(d / "out"),
              {"sherd": {"id": "T_1.1"}, "scale": {"mm_per_unit": 10.0},
               "drawing": {"formats": ["svg", "png"]}}, log=lambda *_: None)
    return d, rep, loop


def test_pipeline_go_and_accurate(sherd_run):
    d, rep, loop = sherd_run
    assert rep["verdict"] == "GO"
    top = loop[:, 1].max()
    outer_rim_r = loop[loop[:, 1] >= top - 1.5, 0].max()      # rim diameter = outer lip edge
    assert abs(rep["profile"]["rim_diameter_mm"] - 2 * outer_rim_r) < 1.0
    t = rep["wall_thickness_mm"]
    assert 3.0 < t["min"] < t["max"] < 9.5


def test_drawing_geometry_is_the_profile(sherd_run):
    """The SVG section polygon must be exactly the profile section at 1:1 mm."""
    d, rep, _ = sherd_run
    svg = (d / "out" / rep["outputs"]["svg"]).read_text()
    prof = json.loads((d / "out" / rep["outputs"]["profile"]).read_text())
    pts = re.search(r'id="section"[^>]*points="([^"]+)"', svg).group(1).split()
    S = np.array([[float(v) for v in p.split(",")] for p in pts])
    P = np.array(prof["section_rh"])
    top = P[:, 1].max()
    assert np.allclose(np.abs(S[:, 0]), P[:, 0], atol=2e-3)
    assert np.allclose(S[:, 1], top - P[:, 1], atol=2e-3)
    m = re.search(r'width="([\d.]+)mm" height="([\d.]+)mm" viewBox="[-\d.]+ [-\d.]+ ([\d.]+) ([\d.]+)"', svg)
    assert abs(float(m.group(1)) - float(m.group(3))) < 0.01   # 1 user unit == 1 mm


def test_flip_and_decorated_flag(tmp_path):
    V, T, _ = sy.make_sherd(120, n=80)
    save_mesh(tmp_path / "s.ply", V, T)
    a = run(str(tmp_path / "s.ply"), str(tmp_path / "a"), {"drawing": {"formats": []}}, log=lambda *_: None)
    b = run(str(tmp_path / "s.ply"), str(tmp_path / "b"),
            {"profile": {"orientation": "flip"}, "ware": "decorated", "drawing": {"formats": []}},
            log=lambda *_: None)
    assert np.dot(a["axis"]["direction"], b["axis"]["direction"]) < -0.99
    assert b["verdict"] == "REVIEW" and any("decorated" in f for f in b["flags"])
