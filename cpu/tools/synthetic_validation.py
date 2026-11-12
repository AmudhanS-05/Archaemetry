"""Full synthetic validation: every case of the self-test, plus small arcs that should
be flagged REVIEW. Writes results.json and prints a table.

Usage: python synthetic_validation.py OUTDIR
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from archaemetry.selftest import main  # noqa: E402

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "synthetic_validation_out"
    ok, res = main(out, full=True)
    json.dump(res, open(os.path.join(out, "results.json"), "w"), indent=2)
    sys.exit(0 if ok else 1)
