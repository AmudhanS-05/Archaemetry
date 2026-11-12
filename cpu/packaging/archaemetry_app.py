"""Frozen-app entry point. No arguments: open the desktop app.
With arguments: behave like the `archaemetry` command (run, selftest, config)."""
import multiprocessing
import sys

from archaemetry.cli import main

if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main(sys.argv[1:] or ["gui"]))
