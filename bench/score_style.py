#!/usr/bin/env python3
"""Plot-style score for an acumen run tree: the house-style fig.* rules, per run.

acumen grades only answer.md. This reads what each run left behind and checks the
code that drew the figure. A run keeps only answer.md, script.py and its transcript
(runner._collect_artifacts), so the code is script.py plus every Write, Edit and
Bash call in transcript.jsonl; that is where an .R file or an Rscript -e call lives.
Figures are not retained; the file check runs only when a run dir holds them.

Usage: python score_style.py [runs_dir] [--csv out.csv] [--selftest]
"""

import argparse
import csv
import json
import re
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent

TITLE = re.compile(
    r"set_title\(|suptitle\(|plt\.title\(|ggtitle\("
    r"|labs\([^)]*\b(?:title|subtitle|caption)\s*=\s*(?!\s|\"\"|'')"
)
CUTOFF = re.compile(r"axhline\(|axvline\(|geom_hline\(|geom_vline\(|plot_volcano\(")
JITTER = re.compile(r"jitter|stripplot|swarmplot|plot_box\(")
SEED = re.compile(
    r"\bseed\s*=|set\.seed\(|random\.seed\(|default_rng\(\s*\d|random_state\s*="
)
TEMPLATES = re.compile(
    r"from\s+(?:plots|theme|palettes)\s+import|import\s+plots\b"
    r"|source\([^)]*(?:plots|theme|palettes)\.R"
)
MM = re.compile(
    r"save_panel\(|width_mm|height_mm|units\s*=\s*[\"']mm|/\s*25\.4|\*\s*0\.03937"
)
LITERAL = re.compile(r"[\"']([^\"'\s]+)\.(pdf|png|csv)[\"']")
CHECKS = ["no_title", "cutoff_lines", "jitter_seed", "triplet", "templates", "mm_size"]


def transcript_code(path):
    """Text of every Write, Edit and Bash tool call in a Claude transcript."""
    if not path.is_file():
        return ""
    parts = []
    for line in path.read_text(errors="replace").splitlines():
        try:
            record = json.loads(line)
        except ValueError:
            continue
        if not isinstance(record, dict) or record.get("type") != "assistant":
            continue
        content = record.get("message", {}).get("content")
        for block in content if isinstance(content, list) else []:
            if not isinstance(block, dict) or block.get("type") != "tool_use":
                continue
            args = block.get("input") or {}
            for key in ("file_path", "content", "new_string", "command"):
                if isinstance(args.get(key), str):
                    parts.append(args[key])
    return "\n".join(parts)


def triplet(code, run_dir):
    """A pdf, png and csv sharing one stem: retained files first, else the code."""
    files = [p for p in run_dir.iterdir() if p.suffix in {".pdf", ".png", ".csv"}]
    if files:
        stems = defaultdict(set)
        for path in files:
            stems[path.stem].add(path.suffix)
        return any(len(exts) == 3 for exts in stems.values())
    if "save_panel(" in code:
        return True
    stems = defaultdict(set)
    for stem, ext in LITERAL.findall(code):
        stems[stem].add(ext)
    return any(len(exts) == 3 for exts in stems.values())


def score_run(run_dir, task_id):
    """Pass (True), fail (False) or not applicable (None) per check."""
    script = run_dir / "script.py"
    code = script.read_text(errors="replace") if script.is_file() else ""
    code += "\n" + transcript_code(run_dir / "transcript.jsonl")
    if not code.strip():
        return dict.fromkeys(CHECKS)
    return {
        "no_title": not TITLE.search(code),
        "cutoff_lines": bool(CUTOFF.search(code)) if "volcano" in task_id else None,
        "jitter_seed": bool(SEED.search(code)) if JITTER.search(code) else None,
        "triplet": triplet(code, run_dir),
        "templates": bool(TEMPLATES.search(code)),
        "mm_size": bool(MM.search(code)),
    }


def collect(runs):
    rows = []
    for run_dir in sorted(runs.glob("*/*/*/*/rep_*")):
        arm, split, model, task_id = run_dir.relative_to(runs).parts[:4]
        result = run_dir / "result.json"
        try:
            success = json.loads(result.read_text()).get("success")
        except (OSError, ValueError):
            success = None
        row = dict(arm=arm, split=split, model=model, task=task_id, rep=run_dir.name)
        rows.append({**row, "success": success, **score_run(run_dir, task_id)})
    return rows


def table(rows):
    """Share of applicable runs passing, per arm x check."""
    lines = [f"{'arm':12s}" + "".join(f"{c:>14s}" for c in CHECKS)]
    for arm in sorted({r["arm"] for r in rows}):
        cells = []
        for check in CHECKS:
            seen = [r[check] for r in rows if r["arm"] == arm and r[check] is not None]
            share = f"{sum(seen) / len(seen):.0%} ({len(seen)})" if seen else "n/a"
            cells.append(f"{share:>14s}")
        lines.append(f"{arm:12s}" + "".join(cells))
    return "\n".join(lines)


def selftest():
    good = """import sys
sys.path.insert(0, ".claude/skills/lucid-coder/templates")
from plots import plot_volcano, plot_box
from theme import save_panel
fig = plot_volcano(res, "lfc", "padj", padj_max=0.1, lfc_min=0.5)
save_panel(fig, res, ".", "volcano", 89, 70)
"""
    bad = """import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(x, y)
ax.set_title("Volcano")
fig.savefig("volcano.png")
"""
    r_write = {
        "type": "assistant",
        "message": {
            "content": [
                {
                    "type": "tool_use",
                    "name": "Write",
                    "input": {
                        "file_path": "box.R",
                        "content": 'p <- ggplot(d) + geom_jitter() + ggtitle("x")\n'
                        'ggsave("box.pdf", p, width = 89, units = "mm")\n'
                        'ggsave("box.png", p)\nwrite.csv(d, "box.csv")\n',
                    },
                }
            ]
        },
    }
    with tempfile.TemporaryDirectory() as tmp:
        runs = Path(tmp)
        for arm, task, body in [
            ("skill_v1", "volcano", good),
            ("noskill", "volcano", bad),
            ("noskill", "box_r", ""),
        ]:
            run_dir = runs / arm / "train" / "claude-sonnet-5" / task / "rep_1"
            run_dir.mkdir(parents=True)
            if body:
                (run_dir / "script.py").write_text(body)
            (run_dir / "result.json").write_text('{"success": true}')
        box = runs / "noskill/train/claude-sonnet-5/box_r/rep_1"
        (box / "transcript.jsonl").write_text(json.dumps(r_write) + "\n")
        rows = {r["arm"] + "/" + r["task"]: r for r in collect(runs)}
    assert all(rows["skill_v1/volcano"][c] for c in CHECKS if c != "jitter_seed")
    assert rows["skill_v1/volcano"]["jitter_seed"] is None
    bad_row = rows["noskill/volcano"]
    assert not any(bad_row[c] for c in ("no_title", "cutoff_lines", "triplet"))
    assert not bad_row["templates"] and not bad_row["mm_size"]
    r_row = rows["noskill/box_r"]
    assert r_row["no_title"] is False and r_row["jitter_seed"] is False
    assert r_row["triplet"] and r_row["mm_size"] and r_row["cutoff_lines"] is None
    assert not TITLE.search('labs(title = "", x = "Genes")')  # blank title is fine
    print("selftest ok")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", nargs="?", type=Path, default=HERE / "plots" / "runs")
    parser.add_argument("--csv", type=Path, help="per-run CSV (default runs/style.csv)")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        return selftest()
    if not args.runs.is_dir():
        sys.exit(f"no run tree at {args.runs}; run `acumen bench` first")
    rows = collect(args.runs)
    if not rows:
        sys.exit(f"no rep_* run dirs under {args.runs}")
    out = args.csv or args.runs / "style.csv"
    with open(out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(table(rows))
    print(f"\n{len(rows)} runs; share passing (applicable runs); per-run CSV: {out}")


if __name__ == "__main__":
    main()
