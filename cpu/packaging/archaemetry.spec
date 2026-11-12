# PyInstaller spec: one-folder app (faster start, fewer antivirus false positives than one-file).
#   pip install "./cpu[build]" && pyinstaller cpu/packaging/archaemetry.spec
import os
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

here = os.path.abspath(SPECPATH)
hidden = collect_submodules("trimesh.exchange") + collect_submodules("archaemetry") + [
    "scipy.sparse.csgraph._validation", "matplotlib.backends.backend_agg",
    "matplotlib.backends.backend_pdf", "matplotlib.backends.backend_svg"]

a = Analysis(
    [os.path.join(here, "archaemetry_app.py")],
    pathex=[os.path.dirname(here)],
    datas=collect_data_files("trimesh"),
    hiddenimports=hidden,
    excludes=["open3d", "torch", "tensorflow", "IPython", "jupyter", "notebook", "PyQt5", "PySide6"],
    noarchive=False,
)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="Archaemetry",
          console=False, disable_windowed_traceback=False)
coll = COLLECT(exe, a.binaries, a.datas, name="Archaemetry")
app = BUNDLE(coll, name="Archaemetry.app", bundle_identifier="org.archaemetry.app") \
    if os.sys.platform == "darwin" else None
