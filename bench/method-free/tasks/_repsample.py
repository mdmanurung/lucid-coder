"""Answer key and prompt for the `repsample` task; adds it to ../tasks.yaml.

Every tasks/repsample-<split>.py is `import _repsample; _repsample.main(__file__)`.
Sound method: sum raw counts per sample, PyDESeq2 `~ donor + condition` on all genes,
then BH over the 10 candidate genes (the first 10 `var_names`).
Run `<target venv>/bin/python tasks/_repsample.py <cache>` to regenerate its tasks.yaml
entry; other tasks' entries are kept as they are.
"""

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_var] = "1"

from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import false_discovery_control  # noqa: E402

SPLITS = ("train", "valid", "test")
N_CAND, ALPHA = 10, 0.05


def h5ad(split):
    return f"repsample_{split}.h5ad"


def pseudobulk(adata):
    """Sample x gene summed counts, and per-sample donor and condition."""
    sample = adata.obs["sample"].astype(str).to_numpy()
    ids = sorted(set(sample))
    counts = pd.DataFrame(
        np.vstack([np.asarray(adata.X[sample == s].sum(axis=0)).ravel() for s in ids]),
        index=ids, columns=adata.var_names,
    ).astype(np.int64)
    meta = adata.obs.drop_duplicates("sample").set_index("sample").loc[ids, ["donor", "condition"]].astype(str)
    return counts, meta


def deseq_pvals(adata, design="~ donor + condition"):
    """PyDESeq2 results (stim vs ctrl) for all genes on per-sample pseudobulk."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    counts, meta = pseudobulk(adata)
    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(counts=counts, metadata=meta, design=design, inference=inference, quiet=True)
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["condition", "stim", "ctrl"], inference=inference, quiet=True)
    stats.summary()
    return stats.results_df


def n_sig(pvals):
    """Candidates with BH-adjusted p < ALPHA, BH over the candidates only; NaN counts as 1."""
    p = np.nan_to_num(np.asarray(pvals, float), nan=1.0)
    return int((false_discovery_control(p) < ALPHA).sum())


def answer(split):
    adata = ad.read_h5ad(h5ad(split))
    res = deseq_pvals(adata)
    return str(n_sig(res.loc[adata.var_names[:N_CAND], "pvalue"]))


def prompt(split, genes, n_donors):
    return (
        f"`{h5ad(split)}` in the working directory holds single-cell RNA-seq of one cell "
        f"type from {n_donors} donors: raw integer UMI counts in `X`, gene names as "
        "`var_names`, and `obs` columns `donor` (donor id), `sample` (sample id) and "
        "`condition` (`ctrl` or `stim`). How many of these 10 genes "
        f"({', '.join(genes)}) differ in expression between stim and ctrl in this cell "
        "type, at an adjusted p-value below 0.05? Report only the number. Genome "
        "reference: none, synthetic genes. There is no one to ask: make each analysis "
        "choice yourself and state it. Report only the number."
    )


def write_tasks(cache):
    import tempfile

    import yaml

    cache, entry = Path(cache).resolve(), {"id": "repsample"}
    for split in SPLITS:
        with tempfile.TemporaryDirectory() as work:
            cwd = Path.cwd()
            os.chdir(work)
            os.symlink(cache / h5ad(split), h5ad(split))
            meta = ad.read_h5ad(h5ad(split), backed="r")
            genes, n_donors = list(meta.var_names[:N_CAND]), meta.obs["donor"].nunique()
            key = answer(split)
            os.chdir(cwd)
        print(f"repsample/{split}: {key}", flush=True)
        entry[split] = {"data": [h5ad(split)], "prompt": prompt(split, genes, n_donors), "answer": key}

    path = Path(__file__).resolve().parent.parent / "tasks.yaml"
    text = path.read_text()
    header = "".join(line + "\n" for line in text.splitlines() if line.startswith("#"))
    tasks = [t for t in yaml.safe_load(text)["tasks"] if t["id"] != "repsample"] + [entry]
    body = "\n".join(
        "\n".join("  " + x for x in yaml.safe_dump([t], sort_keys=False, width=92, allow_unicode=True).splitlines())
        for t in tasks
    )
    path.write_text(header + "tasks:\n" + body + "\n")


def main(script):
    split = Path(script).stem.rsplit("-", 1)[1]
    Path("answer.md").write_text(answer(split) + "\n")


if __name__ == "__main__":
    import sys

    write_tasks(sys.argv[1])
