## Python
Style: `ruff format` at 88 columns; it wraps calls, not strings, so keep messages short. f-strings, `pathlib` with `ROOT = Path(__file__).resolve().parents[1]`, `argparse` with `Path` arguments, `main()` behind `if __name__ == "__main__"`.
Errors: raise `ValueError`/`KeyError` with the fix in the message; `assert` only in tests; no bare `except:`.
Figures: `pdf.fonttype = 42` and `ps.fonttype = 42`; `figsize` in mm / 25.4; `fig.savefig(run_dir / name, bbox_inches="tight")` then `plt.close(fig)`, before any `plt.show()`; `sc.settings.figdir = run_dir` once.
Smells: an AnnData subset is a view, `.copy()` before assigning into `.obs`/`.var`. `validate=` on every merge. `random_state` as an argument. End with `sc.logging.print_versions()` written to the run dir.
