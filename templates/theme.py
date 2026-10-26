"""Figure theme and panel export, shared by every figure script."""

from pathlib import Path

import matplotlib

# Type 42 embeds editable TrueType text instead of outlines; 8 pt matches the
# R theme (theme_bw(base_size = 8)).
matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.size": 8})


def save_panel(fig, data, run_dir, name, width_mm, height_mm):
    """Write <run_dir>/<name>.pdf (vector), .png (review) and .csv (source data).

    `data` is the same frame the panel was drawn from, so the table and the
    panel cannot drift apart.
    """
    run_dir = Path(run_dir)
    if not run_dir.is_dir():
        raise FileNotFoundError(
            f"save_panel: run_dir {run_dir} does not exist; create it before saving."
        )
    fig.set_size_inches(width_mm / 25.4, height_mm / 25.4)
    fig.savefig(run_dir / f"{name}.pdf")
    fig.savefig(run_dir / f"{name}.png", dpi=300)
    data.to_csv(run_dir / f"{name}.csv", index=False)
    return run_dir / name
