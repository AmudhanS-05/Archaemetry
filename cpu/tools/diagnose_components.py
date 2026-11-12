"""List the connected components of a raw SfM mesh (sherd vs turntable vs noise).

Usage: python diagnose_components.py mesh.obj
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaemetry.meshio import load_mesh, face_components  # noqa: E402

V, T = load_mesh(sys.argv[1])
lab, n = face_components(T, len(V))
counts = np.bincount(lab, minlength=n)
print(f"{len(V):,} verts / {len(T):,} faces / {n} components")
for k in np.argsort(counts)[::-1][:10]:
    P = V[np.unique(T[lab == k].ravel())]
    ext = P.max(0) - P.min(0)
    print(f"  comp {k:5d}: {counts[k]:9,} faces  bbox {np.round(ext, 3)}  faces/area-proxy "
          f"{counts[k] / max(np.prod(np.sort(ext)[1:]), 1e-12):.3g}")
