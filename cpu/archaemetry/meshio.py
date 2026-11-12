"""Mesh input/output and connected components (trimesh + scipy; no GPU, no Open3D)."""
import numpy as np
import trimesh
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def load_mesh(path):
    """Read OBJ/PLY/STL/OFF/GLB into (V float64 (N,3), T int64 (M,3)).
    Multi-object files are merged; vertices are not merged or reordered."""
    m = trimesh.load(path, force="mesh", process=False)
    if not isinstance(m, trimesh.Trimesh) or len(m.faces) == 0:
        raise ValueError(f"{path}: no triangle mesh found")
    return np.asarray(m.vertices, float), np.asarray(m.faces, np.int64)


def save_mesh(path, V, T):
    trimesh.Trimesh(V, T, process=False).export(path)


def face_components(T, n_verts):
    """Label triangles by connected component (sharing a vertex). Returns (labels, n)."""
    m = len(T)
    rows = np.repeat(np.arange(m), 3)
    A = coo_matrix((np.ones(3 * m), (rows, T.ravel())), shape=(m, n_verts)).tocsr()
    # faces connected through vertices: run components on the bipartite face-vertex graph
    from scipy.sparse import bmat
    G = bmat([[None, A], [A.T, None]]).tocsr()
    n, lab = connected_components(G, directed=False)
    face_lab = lab[:m]
    _, face_lab = np.unique(face_lab, return_inverse=True)
    return face_lab, int(face_lab.max() + 1) if m else 0


def submesh(V, T, face_mask):
    """Keep the masked triangles and the vertices they use (re-indexed)."""
    Tk = T[face_mask]
    used = np.unique(Tk.ravel())
    remap = -np.ones(len(V), np.int64)
    remap[used] = np.arange(len(used))
    return V[used], remap[Tk]
