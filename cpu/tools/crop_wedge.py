"""Crop a complete vessel mesh to an angular wedge over a height band, to simulate an
open (exterior-only) rim/body fragment with ground truth.

Usage:
  python crop_wedge.py IN.ply OUT.ply --deg 70 --hmin-frac 0.6 --hmax-frac 1.0
"""
import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaemetry.axis import fit_axis  # noqa: E402
from archaemetry.meshio import load_mesh, save_mesh  # noqa: E402
from accuracy_oracle import crop_wedge_faces  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("inp"); ap.add_argument("out")
ap.add_argument("--deg", type=float, default=70.0, help="angular width of wedge")
ap.add_argument("--hmin-frac", type=float, default=0.6)
ap.add_argument("--hmax-frac", type=float, default=1.0)
ap.add_argument("--center-deg", type=float, default=0.0)
a = ap.parse_args()

V, T = load_mesh(a.inp)
p, d, _ = fit_axis(V)
Vs, Ts, _ = crop_wedge_faces(V, T, p, d, a.deg, a.hmin_frac, a.hmax_frac, a.center_deg)
save_mesh(a.out, Vs, Ts)
print(f"wedge: {a.deg} deg, height {a.hmin_frac:.2f}-{a.hmax_frac:.2f}; kept {len(Vs):,} verts / {len(Ts):,} faces")
