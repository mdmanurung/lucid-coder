"""Render the template self-check figures into site/img/ for the docs.

Run: uv run --with plotnine-extra,pandas,numpy python site/make_figures.py
Outputs: site/img/<name>.png, one per figure in templates/plots.py's self-check.
"""

import contextlib
import runpy
import sys
import tempfile
from pathlib import Path

TEMPLATES = Path(__file__).resolve().parents[1] / "templates"
OUT = Path(__file__).resolve().parent / "img"

# ponytail: reuses the self-check by pointing its temp dir at OUT; breaks if the
# self-check stops using TemporaryDirectory, then give plots.py a demo entry point.
tempfile.TemporaryDirectory = lambda: contextlib.nullcontext(str(OUT))
sys.path.insert(0, str(TEMPLATES))
runpy.run_path(str(TEMPLATES / "plots.py"), run_name="__main__")
for path in OUT.iterdir():
    if path.suffix != ".png":
        path.unlink()
