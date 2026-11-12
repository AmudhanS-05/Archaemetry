"""Stage 2: isolate the sherd from turntable, backdrop and noise.

Keeps the connected component with the most triangles (a turntable is a large area
made of few triangles; the photographed sherd is a small area made of very many),
or every component above `min_faces`. Optional capped hole fill via pymeshlab.
"""
import numpy as np
from .meshio import face_components, submesh


def clean(V, T, min_faces=None):
    lab, n = face_components(T, len(V))
    counts = np.bincount(lab, minlength=n)
    if min_faces:
        keep = np.isin(lab, np.where(counts >= min_faces)[0])
    else:
        keep = lab == int(counts.argmax())
    V2, T2 = submesh(V, T, keep)
    ext = V2.max(0) - V2.min(0)
    report = {
        "components_before": n,
        "faces_before": int(len(T)),
        "faces_after": int(len(T2)),
        "kept_fraction": round(len(T2) / max(len(T), 1), 4),
        "bbox_extent": [round(float(x), 4) for x in ext],
        "compactness_ratio": round(float(ext.min() / ext.max()), 3) if ext.max() > 0 else 0.0,
        "components_kept": int(len(np.unique(lab[keep]))),
    }
    return V2, T2, report


def close_holes(V, T, max_edges):
    """Fill holes with at most `max_edges` boundary edges. Needs pymeshlab (optional)."""
    try:
        import pymeshlab
    except ImportError:
        return V, T, "skipped: pymeshlab not installed"
    ms = pymeshlab.MeshSet()
    ms.add_mesh(pymeshlab.Mesh(V, T.astype(np.int32)))
    ms.apply_filter("meshing_close_holes", maxholesize=int(max_edges), selfintersection=False)
    m = ms.current_mesh()
    return m.vertex_matrix(), m.face_matrix().astype(np.int64), "ok"
