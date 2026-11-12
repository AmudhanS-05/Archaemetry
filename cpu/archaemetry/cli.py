"""Command line: archaemetry run | selftest | gui | config"""
import argparse
import sys

from . import __version__


def main(argv=None):
    ap = argparse.ArgumentParser(prog="archaemetry", description="Profile drawings of pottery from 3D meshes.")
    ap.add_argument("--version", action="version", version=f"archaemetry {__version__}")
    sub = ap.add_subparsers(dest="cmd")

    r = sub.add_parser("run", help="mesh -> drawing + report")
    r.add_argument("input", help="mesh file (OBJ, PLY, STL, OFF, GLB)")
    r.add_argument("-c", "--config", help="YAML config (see `archaemetry config`)")
    r.add_argument("-o", "--outdir", default="out")
    r.add_argument("--id", help="sherd ID (default: file name)")
    r.add_argument("--axis", choices=["auto", "pca", "rim_arc", "normals_refine", "manual"])
    r.add_argument("--side", choices=["left", "right"], help="section side")
    r.add_argument("--flip", action="store_true", help="rim is at the other end")
    r.add_argument("--mm-per-unit", type=float, help="scale mesh units to mm")

    s = sub.add_parser("selftest", help="check the install against known geometry")
    s.add_argument("-o", "--outdir")
    s.add_argument("--full", action="store_true")

    sub.add_parser("gui", help="open the desktop app")
    sub.add_parser("config", help="print the default config as YAML")

    a = ap.parse_args(argv)
    if a.cmd == "run":
        from .pipeline import load_config, run
        cfg, _ = load_config(a.config)
        if a.id: cfg["sherd"]["id"] = a.id
        if a.axis: cfg["profile"]["axis_method"] = a.axis
        if a.side: cfg["drawing"]["section_side"] = a.side
        if a.flip: cfg["profile"]["orientation"] = "flip"
        if a.mm_per_unit: cfg["scale"]["mm_per_unit"] = a.mm_per_unit
        rep = run(a.input, a.outdir, cfg)
        return 0 if rep["verdict"] == "GO" else 2
    if a.cmd == "selftest":
        from .selftest import main as st
        ok, _ = st(a.outdir, a.full)
        return 0 if ok else 1
    if a.cmd == "config":
        import yaml
        from .pipeline import DEFAULTS
        print(yaml.safe_dump(DEFAULTS, sort_keys=False))
        return 0
    if a.cmd == "gui" or a.cmd is None:
        from .gui import main as gui
        gui()
        return 0


if __name__ == "__main__":
    sys.exit(main())
