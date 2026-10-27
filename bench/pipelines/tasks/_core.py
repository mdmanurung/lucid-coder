"""Spec, prompts and ground truth for the `core` pipeline tasks (qc, embed, markers, annotate).

Each tasks/<id>-<split>.py reproducer is `import _core; _core.main(__file__)`. Each
procedure mirrors skills/<pipeline>/scripts/*.py but is written independently, so a bug
in a skill cannot hide in its own answer key. Corpus rows are named in task_parts/core.yaml.

    <target venv>/bin/python tasks/_core.py <cache dir>

recomputes every answer in a temp dir linking the cache and rewrites task_parts/core.yaml
and the reproducer stubs.
"""

import os
import random
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from threadpoolctl import threadpool_limits

KANG = "kang2018.h5ad"
PBMC3K = "core_pbmc3k.h5ad"
REF, QUERY, MARKERS = (
    "core_kang_stim_ref.h5ad",
    "core_kang_ctrl_query.h5ad",
    "core_pbmc_markers.tsv",
)

SPEC = {
    "qc": dict(
        skill="sc-qc",
        corpus="B-an-01 (Scrublet per sample, seeded), B-an-02 (fixed mito/gene "
        "limits, drop predicted_doublet)",
        data=[PBMC3K],
        train=dict(min_genes=200, max_genes=2500, min_counts=500, mito=5, rate=0.06, seed=0),
        valid=dict(min_genes=300, max_genes=2000, min_counts=800, mito=10, rate=0.08, seed=1),
        test=dict(min_genes=250, max_genes=3000, min_counts=600, mito=8, rate=0.05, seed=2),
    ),
    "embed": dict(
        skill="sc-embed",
        corpus="B-an-06 (normalise, HVG, scale, arpack PCA, Harmony max_iter 20, kNN, "
        "UMAP), C-an-05 (leiden igraph, n_iterations=2)",
        data=[KANG],
        train=dict(cell_type="CD8 T cells", resolution=0.4),
        valid=dict(cell_type="NK cells", resolution=0.5),
        test=dict(cell_type="B cells", resolution=0.3),
    ),
    "markers": dict(
        skill="sc-markers",
        corpus="D-an-40 (Wilcoxon, benjamini-hochberg, reference='rest', use_raw=False)",
        data=[KANG],
        train=dict(label="ctrl", cell_type="B cells", padj=0.05, lfc=1, pct=0.25),
        valid=dict(label="stim", cell_type="NK cells", padj=0.01, lfc=1.5, pct=0.3),
        test=dict(label="ctrl", cell_type="CD8 T cells", padj=0.05, lfc=0.5, pct=0.2),
    ),
    "annotate": dict(
        skill="sc-annotate",
        corpus="B-an-04 / D-an-41 (train CellTypist on a log1p-1e4 reference by a label "
        "column, annotate with majority voting toggled)",
        data=[REF, QUERY, MARKERS],
        train=dict(cell_type="CD8 T cells"),
        valid=dict(cell_type="NK cells"),
        test=dict(cell_type="Dendritic cells"),
    ),
}

KANG_TEXT = (
    f"`{KANG}` in the working directory holds Kang et al. 2018 PBMCs: raw UMI counts "
    "in `X`; `obs` columns `label` (`ctrl` or `stim`, interferon-beta), `replicate` "
    "(patient id) and `cell_type`."
)


def prompt(task_id, split_name):
    s = SPEC[task_id][split_name]
    if task_id == "qc":
        return (
            f"`{PBMC3K}` in the working directory holds the 10x PBMC 3k library: raw UMI "
            "counts in `X`, gene symbols as `var_names`, and `obs` column `sample` (one "
            "library). Run cell QC and doublet removal. Compute per-cell QC metrics on the "
            "raw counts; mitochondrial genes start with `MT-`. Score doublets on all cells "
            "before any filtering with scanpy's `sc.pp.scrublet`, per `sample`: expected "
            f"doublet rate {s['rate']}, 30 principal components, exact nearest neighbours "
            f"(no approximate search), random state {s['seed']}, and Scrublet's own "
            "automatic threshold; a cell is a doublet when its score exceeds that "
            f"threshold. Keep a cell when it has {s['min_genes']} to {s['max_genes']} "
            f"genes detected (inclusive), at least {s['min_counts']} UMI counts, at most "
            f"{s['mito']}% mitochondrial counts, and is not a predicted doublet. Save the "
            "filtered h5ad, a per-cell QC table and QC figures in the working directory. "
            "Report only the number of cells kept."
        )
    if task_id == "embed":
        return (
            f"{KANG_TEXT} Embed and cluster the {s['cell_type']}, both conditions "
            "together. Keep the raw counts in a layer; normalise each cell to 10,000 "
            "total counts, then log1p. Pick the 2,000 most highly variable genes with the "
            "`seurat` flavor on the log-normalised data, with `replicate` as the batch "
            "key. Scale those genes to zero mean and unit variance with no clipping and "
            "run PCA with 30 components and the ARPACK solver. Integrate the PCs by "
            "`replicate` with harmonypy's `run_harmony` (at most 20 iterations). Build a "
            "15-neighbour kNN graph on the first 20 Harmony-corrected components, compute "
            "a UMAP with min_dist 0.3, and run Leiden clustering (igraph flavor, 2 "
            f"iterations, undirected graph) at resolution {s['resolution']}. Use seed 0 "
            "for PCA, Harmony, neighbours, UMAP and Leiden. Save the embedded h5ad, a "
            "cluster-size table and UMAP figures in the working directory. Report only "
            "the number of Leiden clusters."
        )
    if task_id == "markers":
        return (
            f"{KANG_TEXT} Find marker genes for each `cell_type` among the "
            f"{s['label']} cells only. Normalise each cell to 10,000 total counts, then "
            "log1p. Test each cell type against all other cells with a Wilcoxon rank-sum "
            "test (scanpy `rank_genes_groups`, every gene, Benjamini-Hochberg correction, "
            "on the log-normalised values, not `.raw`). A gene is a marker of a cell type "
            f"when its adjusted p < {s['padj']}, its log2 fold change (scanpy's, up in "
            f"that type) > {s['lfc']}, and it is detected (value > 0) in at least "
            f"{s['pct']} of that type's cells. Save the full marker table, a dot plot and "
            "a heatmap of the top 5 markers per cell type in the working directory. "
            f"Report only the number of markers of {s['cell_type']}."
        )
    return (
        f"`{REF}` in the working directory holds the interferon-beta stimulated PBMCs of "
        "Kang et al. 2018: raw UMI counts in `X`, `obs` column `cell_type`. "
        f"`{QUERY}` holds the control PBMCs: raw UMI counts in `X`, a UMAP in "
        "`obsm['X_umap']`, `obs` columns `replicate` and `cell_type`. "
        f"`{MARKERS}` lists marker genes per label (columns `label`, `gene`). Transfer "
        "labels from the stimulated cells to the control cells with CellTypist: "
        "normalise copies of both to 10,000 counts per cell, then log1p; train a model "
        "on the stimulated cells by `cell_type` with `celltypist.train` (default "
        "solver, n_jobs 1, random state 0); annotate the control cells with it in "
        "best-match mode without majority voting. Tabulate the transferred labels by "
        "the control cells' `cell_type`, and check them against the marker genes in a "
        "dot plot. Save the annotated h5ad, the tables and the figures in the working "
        f"directory. Report only the number of control cells with `cell_type` "
        f"{s['cell_type']} whose transferred label is also {s['cell_type']}."
    )


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)


def answer_qc(s):
    import scanpy as sc

    adata = ad.read_h5ad(PBMC3K)
    x = sp.csr_matrix(adata.X)
    genes = np.asarray((x > 0).sum(axis=1)).ravel()
    counts = np.asarray(x.sum(axis=1)).ravel()
    mito = adata.var_names.str.startswith("MT-")
    pct_mito = 100 * np.asarray(x[:, mito].sum(axis=1)).ravel() / counts
    passed = (
        (genes >= s["min_genes"])
        & (genes <= s["max_genes"])
        & (counts >= s["min_counts"])
        & (pct_mito <= s["mito"])
    )
    cells = ad.AnnData(x.astype(np.float32), obs=adata.obs[["sample"]].astype(str))
    sc.pp.scrublet(
        cells,
        batch_key="sample",
        expected_doublet_rate=s["rate"],
        n_prin_comps=30,
        use_approx_neighbors=False,
        random_state=s["seed"],
        verbose=False,
    )
    cutoff = cells.obs["sample"].map(
        {b: v["threshold"] for b, v in cells.uns["scrublet"]["batches"].items()}
    )
    doublet = (cells.obs["doublet_score"] > cutoff.astype(float)).to_numpy()
    return str(int((passed & ~doublet).sum()))


def answer_embed(s):
    import harmonypy
    import scanpy as sc

    seed_all(0)
    kang = ad.read_h5ad(KANG)
    adata = kang[kang.obs["cell_type"].astype(str) == s["cell_type"]].copy()
    adata.obs["replicate"] = adata.obs["replicate"].astype(str).astype("category")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, n_top_genes=2000, flavor="seurat", batch_key="replicate"
    )
    hvg = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(hvg)
    sc.tl.pca(hvg, n_comps=30, svd_solver="arpack", random_state=0)
    corrected = harmonypy.run_harmony(
        hvg.obsm["X_pca"],
        adata.obs,
        "replicate",
        max_iter_harmony=20,
        random_state=0,
        ncores=1,
        verbose=False,
    ).Z_corr
    adata.obsm["X_harmony"] = np.asarray(corrected)
    sc.pp.neighbors(
        adata, n_neighbors=15, n_pcs=20, use_rep="X_harmony", random_state=0
    )
    sc.tl.leiden(
        adata,
        resolution=s["resolution"],
        flavor="igraph",
        n_iterations=2,
        directed=False,
        random_state=0,
    )
    return str(adata.obs["leiden"].nunique())


def answer_markers(s):
    import scanpy as sc

    kang = ad.read_h5ad(KANG)
    adata = kang[kang.obs["label"].astype(str) == s["label"]].copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        "cell_type",
        groups=[s["cell_type"]],
        method="wilcoxon",
        reference="rest",
        corr_method="benjamini-hochberg",
        pts=True,
        use_raw=False,
        n_genes=adata.n_vars,
    )
    res = sc.get.rank_genes_groups_df(adata, group=s["cell_type"])
    marker = (
        (res["pvals_adj"] < s["padj"])
        & (res["logfoldchanges"] > s["lfc"])
        & (res["pct_nz_group"] >= s["pct"])
    )
    return str(int(marker.sum()))


def lognorm(adata):
    import scanpy as sc

    out = adata.copy()
    out.X = sp.csr_matrix(out.X, dtype=np.float32)
    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)
    return out


def answer_annotate(s):
    import celltypist

    ref, query = ad.read_h5ad(REF), ad.read_h5ad(QUERY)
    model = celltypist.train(
        lognorm(ref), labels="cell_type", n_jobs=1, random_state=0
    )
    pred = celltypist.annotate(lognorm(query), model=model, majority_voting=False)
    label = pred.predicted_labels["predicted_labels"].astype(str).to_numpy()
    truth = query.obs["cell_type"].astype(str).to_numpy()
    return str(int(((truth == s["cell_type"]) & (label == s["cell_type"])).sum()))


ANSWER = {
    "qc": answer_qc,
    "embed": answer_embed,
    "markers": answer_markers,
    "annotate": answer_annotate,
}


def answer(task_id, split_name):
    import scanpy as sc

    sc.settings.n_jobs = 1
    with threadpool_limits(1):
        return ANSWER[task_id](SPEC[task_id][split_name])


def main(script):
    task_id, split_name = Path(script).stem.rsplit("-", 1)
    Path("answer.md").write_text(answer(task_id, split_name) + "\n")


def write_part(cache):
    import yaml

    here = Path(__file__).resolve().parent
    lines = [
        "# core group: sc-qc, sc-embed, sc-markers, sc-annotate; generated by",
        "# `<target venv>/bin/python tasks/_core.py <cache>` from tasks/_core.py SPEC.",
        "tasks:",
    ]
    for task_id, task in SPEC.items():
        entry = {"id": task_id, "data": task["data"]}
        for split_name in ("train", "valid", "test"):
            with tempfile.TemporaryDirectory() as work:
                os.chdir(work)
                for name in task["data"]:
                    os.symlink(Path(cache).resolve() / name, name)
                result = answer(task_id, split_name)
            entry[split_name] = {"prompt": prompt(task_id, split_name) + " Genome reference: GRCh37 (hg19), 10x Cell Ranger hg19 annotation.", "answer": result}
            (here / f"{task_id}-{split_name}.py").write_text(
                '"""Reproducer; the procedure lives in _core.py."""\n\n'
                "import _core\n\n_core.main(__file__)\n"
            )
            print(f"{task_id}/{split_name}: {result}", flush=True)
        body = yaml.safe_dump([entry], sort_keys=False, width=88)
        lines.append(f"  # skill: {task['skill']}; corpus: {task['corpus']}")
        lines += ["  " + line for line in body.splitlines()]
    (here.parent / "task_parts" / "core.yaml").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    write_part(sys.argv[1])
