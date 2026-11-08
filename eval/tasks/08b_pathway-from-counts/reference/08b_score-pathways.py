"""Purpose: test IFN-beta stimulation in CD14+ monocyte pseudobulk counts and score
         PROGENy pathway activity on the result.
Run:     uv run --python 3.12 --with-requirements requirements.txt python \
         scripts/08b_score-pathways.py results/run1
Inputs:  data/counts.tsv    [gap: pseudobulk counts, gene x sample, given with the task]
         data/samples.tsv   [gap: sample, replicate (donor), label, given with the task]
         data/progeny.tsv   [gap: PROGENy top 100 footprints, given with the task]
Outputs: <run_dir>/de.tsv            PyDESeq2 table per gene, stim vs ctrl
         <run_dir>/pathways.tsv      ulm score and padj per pathway
         <run_dir>/session_info.txt  package versions
Notes:   [Stochastic] none. Design ~ label, as in sc-best-practices; the donor-paired
         ~ replicate + label gives JAK-STAT 14.89 instead of 16.72. The pathway input
         is the signed Wald `stat`: it carries the direction and the precision of each
         change; log2FoldChange ignores precision and p-values lose the sign. padj is
         decoupler's Benjamini-Hochberg across the 13 pathways scored.
"""

import argparse
import sys
from importlib.metadata import version
from pathlib import Path

import decoupler as dc
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

ROOT = Path(__file__).resolve().parents[1]
DESIGN = "~ label"  # as in sc-best-practices; "~ replicate + label" pairs donors
CONTRAST = ["label", "stim", "ctrl"]  # positive means higher after stimulation


def main():
    # Setup --------------------------------------------------------------------
    # Make a fresh run dir.
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path, help="new directory under results/")
    run_dir = parser.parse_args().run_dir
    if run_dir.exists():
        raise FileExistsError(f"run_dir {run_dir} already exists; pick a new name.")
    run_dir.mkdir(parents=True)

    # Load inputs --------------------------------------------------------------
    # Counts as samples x genes for PyDESeq2, samples aligned by name.
    counts = pd.read_csv(ROOT / "data" / "counts.tsv", sep="\t", index_col="gene").T
    samples = pd.read_csv(ROOT / "data" / "samples.tsv", sep="\t", index_col="sample")
    progeny = pd.read_csv(ROOT / "data" / "progeny.tsv", sep="\t")
    if not counts.index.equals(samples.index):
        raise ValueError("samples.tsv rows do not match the counts columns")

    # Test stim vs ctrl --------------------------------------------------------
    # PyDESeq2 fit and Wald test.
    dds = DeseqDataSet(counts=counts, metadata=samples, design=DESIGN, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=CONTRAST, quiet=True)
    stats.summary()
    de = stats.results_df.rename_axis("gene").sort_values("pvalue")

    # Score pathways on the Wald stat ------------------------------------------
    # One ulm score per pathway.
    wald = de[["stat"]].T.rename(index={"stat": "stim_vs_ctrl"})
    score, padj = dc.mt.ulm(data=wald, net=progeny)
    pathways = pd.DataFrame(
        {"pathway": score.columns, "score": score.iloc[0], "padj": padj.iloc[0]}
    ).sort_values("score", ascending=False)

    # Write outputs ------------------------------------------------------------
    # One write line per table.
    de.to_csv(run_dir / "de.tsv", sep="\t")
    pathways.to_csv(run_dir / "pathways.tsv", sep="\t", index=False)

    # Session info -------------------------------------------------------------
    # Versions last.
    session_info = [f"python {sys.version}"] + [
        f"{pkg} {version(pkg)}" for pkg in ["decoupler", "pydeseq2", "pandas"]
    ]
    (run_dir / "session_info.txt").write_text("\n".join(session_info) + "\n")


if __name__ == "__main__":
    main()
