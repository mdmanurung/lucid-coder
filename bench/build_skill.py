"""Build skills/v1/ for an acumen project from the lucid-coder skill and plot templates.

Usage: python build_skill.py <project_dir>   (e.g. plots)

acumen reads a version's SKILL.md at skills/v1/SKILL.md and installs the folder's
content as .claude/skills/<skill_name>/ in each sandbox, so the files go straight
into v1/. acumen's fit skips creating v1 when it exists, so v1 is our skill.
"""

import re
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
TEMPLATES = [
    f"{stem}.{ext}" for stem in ("plots", "theme", "palettes") for ext in ("py", "R")
]
NOTE = """
## Plot templates

The plot templates are in `templates/` next to this file. In Python, add that
directory to `sys.path` and import them:
`sys.path.insert(0, ".claude/skills/lucid-coder/templates")`, then
`from plots import plot_volcano` and `from theme import save_panel`. In R,
`source()` `theme.R`, `palettes.R` and `plots.R` from
`.claude/skills/lucid-coder/templates/`, in that order.
"""


PIPELINE_NOTE = """
## Where this skill lives

This skill is installed at `.claude/skills/{name}/`. Run its script from the working
directory as `python .claude/skills/{name}/scripts/{script}`.
"""
# Pipeline skills ship their own script and vendored templates; tests stay behind.
PIPELINES = {"pseudobulk-de": "pseudobulk_de.py"}

UMBRELLA = """---
name: sc-pipelines
description: >-
  Router to tested, deterministic scripts for single-cell and bulk-omics analyses
  whose exact numeric answer depends on parameter and formula choices an ad-hoc
  script is likely to get slightly wrong: QC and doublet filtering, normalisation,
  embedding, clustering, cell-type annotation, marker genes (rank_genes_groups),
  pseudobulk differential expression, cell-type composition and differential
  abundance, group-comparison or correlation tests (Wilcoxon, Spearman, Pearson),
  longitudinal mixed models, gene-set enrichment, gene-program discovery (consensus
  NMF, cNMF), gene-signature scoring, TF or pathway activity, cell-cell
  communication, trajectory and pseudotime, RNA velocity, MOFA multi-omics factor
  analysis, mashr effect sharing, and variance partitioning. Fire on any of those
  analyses, even when the task already spells out its own step-by-step method --
  check whether a fixed script covers it before writing one by hand.
  (acumen bench 2026-10-27: this wording raised script use from 30/48 to 39/48 runs.)
---

# sc-pipelines

Each analysis below has a tested script that fixes every step. Do not write your own
analysis code when a pipeline covers the task.

1. Pick the pipeline from the table.
2. Read its `SKILL.md`, at `.claude/skills/sc-pipelines/<pipeline>/SKILL.md`.
3. Run its script as `python .claude/skills/sc-pipelines/<pipeline>/scripts/<script>`
   with every required argument. Take each value from the task; never invent one.
4. Read the answer from the script's `result.json` and `tables/`. Do not recompute it
   by hand, even for a question about one intermediate step.

| Pipeline | Script | Use it for |
|---|---|---|
{rows}
"""


def first_sentence(skill_md):
    """The frontmatter description's first sentence, for the routing table."""
    text = skill_md.read_text().split("---")[1]
    desc = " ".join(text.split("description:", 1)[1].replace(">-", "").split())
    first = re.split(r"(?<!e\.g)(?<!i\.e)\. (?=[A-Z])", desc)[0]
    return first.split(" Fire on")[0].rstrip(".") + "."


def build_umbrella(dest, slim=False):
    """dest/SKILL.md routes to every pipeline skill, copied as dest/<name>/.

    slim: no examples/ or references/ (maintainer notes), and a router row cut at
    its first colon.
    """
    rows = []
    for src in sorted((REPO / "skills").glob("*/scripts")):
        skill = src.parent
        script = f"{skill.name.replace('-', '_')}.py"
        if not (src / script).is_file() or not (skill / "SKILL.md").is_file():
            continue
        skip = (
            ("tests", "__pycache__", "examples", "references")
            if slim
            else ("tests", "__pycache__")
        )
        shutil.copytree(skill, dest / skill.name, ignore=shutil.ignore_patterns(*skip))
        use_for = first_sentence(skill / "SKILL.md")
        if slim:
            use_for = use_for.split(": ")[0].rstrip(".") + "."
        rows.append(f"| `{skill.name}` | `{script}` | {use_for} |")
    (dest / "SKILL.md").write_text(UMBRELLA.format(rows="\n".join(rows)))
    print(f"built {dest} with {len(rows)} pipelines")


def main(project, version="v1", slim=False):
    v1 = REPO / "bench" / Path(project).name / "skills" / version  # any cwd
    if v1.exists():
        sys.exit(f"{v1} exists; delete it to rebuild (acumen treats v1 as frozen)")
    name = Path(project).name
    if name in (
        "pipelines",
        "method-free",
    ):  # projects benched against the umbrella router
        v1.mkdir(parents=True)
        build_umbrella(v1, slim)
        return
    if name in PIPELINES:
        shutil.copytree(
            REPO / "skills" / name,
            v1,
            ignore=shutil.ignore_patterns("tests", "__pycache__"),
        )
        with open(v1 / "SKILL.md", "a") as skill_md:
            skill_md.write(PIPELINE_NOTE.format(name=name, script=PIPELINES[name]))
        print(f"built {v1}")
        return
    shutil.copytree(REPO / "skills" / "lucid-coder", v1)
    (v1 / "templates").mkdir()
    for name in TEMPLATES:
        shutil.copy2(REPO / "templates" / name, v1 / "templates" / name)
    with open(v1 / "SKILL.md", "a") as skill_md:
        skill_md.write(NOTE)
    print(f"built {v1}")


if __name__ == "__main__":
    # build_skill.py <project> [--slim vN]: --slim builds a lean umbrella as version vN.
    if len(sys.argv) == 4 and sys.argv[2] == "--slim":
        main(sys.argv[1], sys.argv[3], slim=True)
    elif len(sys.argv) == 2:
        main(sys.argv[1])
    else:
        sys.exit(__doc__)
