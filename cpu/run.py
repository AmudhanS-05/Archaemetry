"""Backward-compatible entry point. Prefer `archaemetry run` (pip install ./cpu).

  python run.py --config config/sherd.example.yaml --input mesh.obj --outdir out/
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from archaemetry.cli import main  # noqa: E402

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config")
    ap.add_argument("--input", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    argv = ["run", a.input, "-o", a.outdir] + (["-c", a.config] if a.config else [])
    sys.exit(main(argv))
