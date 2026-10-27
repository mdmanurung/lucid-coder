"""Shared spec and ground truth for the `struct` pipeline tasks: traj, programs, mofa, pbde.

Every tasks/<id>-<split>.py is `import _structg; _structg.main(__file__)` (not `_struct`: that name is a
stdlib C module, which shadows a file of the same name). Each answer
function mirrors one pipeline skill's fixed procedure (skills/trajectory,
gene-programs, mofa, pseudobulk-de) but is written independently: nothing here imports
a skill script, so a skill bug cannot hide in its own answer key.
"""

import os

# One BLAS/numba thread before numpy loads: the same numbers on any host.
for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ[_var] = "1"

from pathlib import Path  # noqa: E402

import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import sparse  # noqa: E402

SPEC = {
    # views (struct_mofa/<view>.tsv), factors, R2 % for an active factor, seed
    "mofa": {
        "train": dict(views=["b", "mono14", "nk"], factors=5, min_r2=15, seed=0),
        "valid": dict(views=["cd4t", "cd8t", "mono16"], factors=5, min_r2=5, seed=1),
        "test": dict(views=["b", "nk"], factors=4, min_r2=15, seed=2),
    },
    # cell type, min cells per sample, min detection fraction, padj and |LFC| cut-offs
    "pbde": {
        "train": dict(cell_type="FCGR3A+ Monocytes", min_cells=20, min_frac=0.2, padj=0.01, lfc=1.5),
        "valid": dict(cell_type="CD4 T cells", min_cells=100, min_frac=0.15, padj=0.05, lfc=2.0),
        "test": dict(cell_type="Dendritic cells", min_cells=15, min_frac=0.25, padj=0.1, lfc=1.0),
    },
    # root group, n diffusion components, PAGA edge threshold, seed
    # cell type, patients kept, k, n HVGs, NMF replicates, seed
    "programs": {
        "train": dict(cell_type="B cells", patients=["patient_1015", "patient_1488"], k=3, n_hvg=500, n_replicates=12, seed=0),
        "valid": dict(cell_type="NK cells", patients=["patient_1256", "patient_1488"], k=4, n_hvg=800, n_replicates=12, seed=0),
        "test": dict(cell_type="CD14+ Monocytes", patients=["patient_1015", "patient_1244"], k=5, n_hvg=1000, n_replicates=12, seed=0),
    },
    "traj": {
        "train": dict(root="g1", n_dcs=10, threshold=0.2, seed=0),
        "valid": dict(root="g2", n_dcs=8, threshold=0.42, seed=1),
        "test": dict(root="g4", n_dcs=12, threshold=0.33, seed=2),
    },
}


def traj(s):
    import scanpy as sc

    adata = ad.read_h5ad("struct_traj.h5ad")
    sc.tl.paga(adata, groups="cluster")
    levels = list(adata.obs["cluster"].cat.categories)
    conn = adata.uns["paga"]["connectivities"].toarray()
    i, j = np.triu_indices(len(levels), k=1)
    n_edges = int((conn[i, j] >= s["threshold"]).sum())
    sc.tl.diffmap(adata, n_comps=s["n_dcs"], random_state=s["seed"])
    adata.uns["iroot"] = int(np.flatnonzero(adata.obs["cluster"] == s["root"])[0])
    sc.tl.dpt(adata, n_dcs=s["n_dcs"])
    dpt = adata.obs.loc[np.isfinite(adata.obs["dpt_pseudotime"])]
    far = dpt.groupby("cluster", observed=True)["dpt_pseudotime"].median().idxmax()
    return f"{n_edges}/{far}"


def programs(s):
    """cNMF-style consensus NMF; top-scoring gene of the program with highest stim AUC."""
    import scanpy as sc
    from sklearn.cluster import KMeans
    from sklearn.decomposition import NMF, non_negative_factorization
    from sklearn.metrics import roc_auc_score

    adata = ad.read_h5ad("kang2018.h5ad")
    keep = (adata.obs["cell_type"].astype(str) == s["cell_type"]) & adata.obs[
        "replicate"
    ].astype(str).isin(s["patients"])
    adata = adata[keep.to_numpy()]
    counts = sparse.csr_matrix(adata.X, dtype=np.float64)

    # HVGs: Seurat dispersion on normalize_total(1e4) + log1p of the expressed genes.
    expressed = np.flatnonzero(np.asarray(counts.sum(axis=0)).ravel() > 0)
    lognorm = ad.AnnData(counts[:, expressed].copy())
    sc.pp.normalize_total(lognorm, target_sum=1e4)
    sc.pp.log1p(lognorm)
    sc.pp.highly_variable_genes(lognorm, flavor="seurat", n_top_genes=s["n_hvg"])
    hvg = expressed[lognorm.var["highly_variable"].to_numpy()]
    x = counts[:, hvg].toarray()
    assert (x.sum(axis=1) > 0).all(), "a cell has zero counts over the HVGs"
    x = x / x.std(axis=0)  # unit variance per gene, not centred

    k, n_rep = s["k"], s["n_replicates"]
    seeds = np.random.default_rng(s["seed"]).integers(1, 2**31 - 1, n_rep)
    spectra = np.vstack(
        [
            NMF(k, init="random", solver="cd", max_iter=1000, random_state=int(r))
            .fit(x)
            .components_
            for r in seeds
        ]
    )
    # Consensus: unit-norm components, k-means (10 starts), median per cluster.
    # The density filter keeps everything at threshold 2 (max unit-vector distance is 2).
    unit = spectra / np.linalg.norm(spectra, axis=1, keepdims=True)
    labels = KMeans(k, n_init=10, random_state=s["seed"]).fit_predict(unit)
    centres = np.vstack([np.median(unit[labels == c], axis=0) for c in range(k)])
    centres /= centres.sum(axis=1, keepdims=True)
    usage, _, _ = non_negative_factorization(
        x, H=centres, n_components=k, update_H=False, solver="cd", max_iter=1000
    )
    usage /= usage.sum(axis=1, keepdims=True)

    # Gene scores: OLS of per-gene z-scored TPM (all expressed genes) on usage.
    tpm = counts[:, expressed].toarray()
    tpm = tpm / tpm.sum(axis=1, keepdims=True) * 1e6
    z = (tpm - tpm.mean(axis=0)) / tpm.std(axis=0)
    scores = np.linalg.pinv(usage) @ z

    stim = (adata.obs["label"].astype(str) == "stim").to_numpy()
    best = int(np.argmax([roc_auc_score(stim, usage[:, c]) for c in range(k)]))
    # Programs are named GP1..GPk by total usage over cells, largest first.
    rank = 1 + list(np.argsort(-usage.sum(axis=0), kind="stable")).index(best)
    return f"GP{rank}/{adata.var_names[expressed][np.argmax(scores[best])]}"


def mofa(s):
    """MOFA2 on Gaussian views; `<view with largest Factor1 R2>/<n active factors>`."""
    import contextlib
    import sys
    import tempfile

    import h5py
    from mofapy2.run.entry_point import entry_point

    tables = {
        v: pd.read_csv(f"struct_mofa/{v}.tsv", sep="\t", index_col=0).sort_index(axis=1)
        for v in sorted(s["views"])
    }
    samples = sorted(set().union(*(t.index for t in tables.values())))
    model = Path(tempfile.mkdtemp()) / "model.hdf5"
    with contextlib.redirect_stdout(sys.stderr):
        ent = entry_point()
        ent.set_data_options(scale_views=True)
        ent.set_data_matrix(
            [[t.reindex(samples).to_numpy(float)] for t in tables.values()],
            likelihoods=["gaussian"] * len(tables),
            views_names=list(tables),
            samples_names=[samples],
            features_names=[list(t.columns) for t in tables.values()],
        )
        ent.set_model_options(factors=s["factors"])
        ent.set_train_options(
            convergence_mode="fast", seed=s["seed"], dropR2=None, gpu_mode=False, quiet=True
        )
        ent.build()
        ent.run()
        ent.save(str(model))  # sorts factors by variance explained
    with h5py.File(model, "r") as h5:
        group = h5["groups/groups"][0].decode()
        views = [v.decode() for v in h5["views/views"][:]]
        r2 = pd.DataFrame(h5[f"variance_explained/r2_per_factor/{group}"][:], index=views)
    n_active = int((r2.max(axis=0) >= s["min_r2"]).sum())
    return f"{r2[0].idxmax()}/{n_active}"


def pbde(s):
    """Paired pseudobulk PyDESeq2; `<n significant up>/<n significant down>`."""
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.default_inference import DefaultInference
    from pydeseq2.ds import DeseqStats

    adata = ad.read_h5ad("kang2018_pb.h5ad")
    cells = adata[(adata.obs["cell_type"].astype(str) == s["cell_type"]).to_numpy()]
    obs = cells.obs[["sample", "replicate", "label"]].astype(str)
    n_cells = obs["sample"].value_counts()
    kept = sorted(n_cells.index[n_cells >= s["min_cells"]])
    in_kept = obs["sample"].isin(kept).to_numpy()
    # Detection over the cells of the kept samples.
    genes = np.asarray((cells.X[in_kept] > 0).mean(axis=0)).ravel() >= s["min_frac"]
    x, labels = cells.X[in_kept][:, genes], obs["sample"].to_numpy()[in_kept]
    counts = pd.DataFrame(
        {k: np.asarray(x[labels == k].sum(axis=0)).ravel() for k in kept},
        index=cells.var_names[genes],
    ).T.astype(np.int64)
    meta = obs.drop_duplicates("sample").set_index("sample").loc[kept]

    inference = DefaultInference(n_cpus=1)
    dds = DeseqDataSet(
        counts=counts, metadata=meta, design="~ replicate + label",
        inference=inference, quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds, contrast=["label", "stim", "ctrl"], alpha=s["padj"],
        inference=inference, quiet=True,
    )
    stats.summary()
    res = stats.results_df
    sig = res[(res["padj"] < s["padj"]) & (res["log2FoldChange"].abs() > s["lfc"])]
    return f"{(sig['log2FoldChange'] > 0).sum()}/{(sig['log2FoldChange'] < 0).sum()}"


ANSWER = {"traj": traj, "programs": programs, "mofa": mofa, "pbde": pbde}


KANG = (
    "Kang et al. 2018 PBMCs: raw UMI counts in `X`; `obs` columns `label` (`ctrl` or "
    "`stim`, interferon-beta), `replicate` (patient id) and `cell_type`"
)
SAVE = "Save the result tables and figures in the working directory."


def prompt(task_id, split):
    s = SPEC[task_id][split]
    if task_id == "traj":
        return (
            "`struct_traj.h5ad` in the working directory holds an embedded single-cell "
            "dataset: log-normalised values in `X`, a neighbour graph (`uns['neighbors']`, "
            "`obsp['connectivities']`), a UMAP in `obsm['X_umap']`, and `obs` column "
            "`cluster` (groups `g1` to `g7`). Infer the trajectory with scanpy, keeping the "
            "stored neighbour graph (do not recompute it). Run PAGA on `cluster`. Run a "
            f"diffusion map with {s['n_dcs']} components and `random_state` {s['seed']}. "
            "The root cell is the first cell of group "
            f"`{s['root']}` in obs order; run diffusion pseudotime (DPT) from it with "
            f"{s['n_dcs']} diffusion components and no branching. A pair of groups is a "
            f"PAGA edge when its connectivity is >= {s['threshold']:g}. Cells with "
            "infinite pseudotime are left out of the group medians. Save the pseudotime "
            "table, the PAGA connectivity table and figures (PAGA graph, UMAP coloured by "
            "pseudotime, pseudotime per group) in the working directory. Report only "
            "`<number of PAGA edges>/<group with the highest median pseudotime>`, for "
            "example `3/g5`."
        )
    if task_id == "programs":
        k, n, seed = s["k"], s["n_replicates"], s["seed"]
        return (
            f"`kang2018.h5ad` in the working directory holds {KANG}. Find gene programs "
            f"by consensus NMF (cNMF-style) in the {s['cell_type']} of patients "
            f"{' and '.join(s['patients'])} (`replicate`). (1) HVGs: among genes with "
            f"non-zero total counts in this subset, the top {s['n_hvg']} by scanpy "
            "`highly_variable_genes(flavor='seurat')` on a copy run through "
            "`normalize_total(target_sum=1e4)` and `log1p`. (2) NMF input: the raw counts "
            "of those HVGs, each gene divided by its standard deviation over cells "
            "(population sd, not centred). (3) Fit scikit-learn `NMF` "
            f"{n} times with {k} components, `init='random'`, `solver='cd'`, Frobenius "
            "loss, `max_iter=1000`; run i uses `random_state` equal to element i of "
            f"`numpy.random.default_rng({seed}).integers(1, 2**31 - 1, {n})`. "
            f"(4) Consensus: L2-normalise each of the {k} x {n} components and keep all of "
            "them (no density filter); cluster them with `KMeans(n_clusters="
            f"{k}, n_init=10, random_state={seed})`; the element-wise median of each "
            "cluster, scaled to sum 1, is a consensus spectrum. (5) Usage: "
            "`non_negative_factorization` of the NMF input with the consensus spectra "
            "fixed (`update_H=False`, `solver='cd'`, `max_iter=1000`), then scale each "
            f"cell's usage to sum 1. (6) Name the programs GP1 to GP{k} by total usage "
            "over cells, largest first. (7) Gene scores, for every gene expressed in the "
            "subset: z-score each gene's TPM (counts scaled to 1e6 per cell; population "
            "sd) over cells, then take the least-squares coefficients `pinv(usage) @ Z`. "
            "The program most active in stim cells is the one whose usage has the "
            "highest one-vs-rest ROC AUC for `label` stim. " + SAVE + " Report only "
            "`<program name>/<its highest-scoring gene>`, for example `GP3/CD74`."
        )
    if task_id == "mofa":
        views = ", ".join(f"`{v}`" for v in s["views"])
        return (
            "The directory `struct_mofa/` in the working directory holds pseudobulk views "
            "of Kang et al. 2018 PBMCs, 16 samples (8 patients x `ctrl`/`stim` "
            "interferon-beta): one TSV per cell type named `<view>.tsv`, samples in rows "
            "(first column `sample`), log2(CPM + 1) features named `<view>_<gene>`; and "
            "`meta.tsv` with columns `sample`, `label` and `replicate`. Run MOFA2 "
            f"(mofapy2) on the views {views}: every view Gaussian, {s['factors']} factors, "
            "views scaled to unit variance, convergence mode `fast`, seed "
            f"{s['seed']}, no factor dropped during training (`dropR2=None`), on CPU. Sort "
            "views, samples and features before building the model. Number the factors in "
            "the order of the saved model (mofapy2 sorts them by variance explained when "
            "saving), so Factor1 explains the most. A factor is active when its variance "
            f"explained (R², in %) is >= {s['min_r2']:g} in at least one view. Colour and "
            "test factors by `label`. " + SAVE + " Report only `<view in which Factor1 "
            "has the largest R²>/<number of active factors>`, for example `cd8t/3`."
        )
    return (
        "`kang2018_pb.h5ad` in the working directory holds Kang et al. 2018 PBMCs: raw "
        "UMI counts in `X`; `obs` columns `sample` (one id per patient x condition), "
        "`replicate` (patient), `label` (`ctrl` or `stim`, interferon-beta) and "
        f"`cell_type`. Pseudobulk the {s['cell_type']}: sum raw counts per `sample` over "
        f"that cell type's cells; drop samples with fewer than {s['min_cells']} cells; "
        f"keep genes detected (count > 0) in at least {s['min_frac']:g} of the cells of "
        "the kept samples. Then test stim vs ctrl with PyDESeq2: design `~ replicate + "
        "label`, no fold-change shrinkage, independent-filtering alpha "
        f"{s['padj']:g}. A gene is significant when padj < {s['padj']:g} and |log2 fold "
        f"change| > {s['lfc']:g}; a missing padj is not significant. " + SAVE + " Report "
        "only `<significant genes up in stim>/<significant genes down in stim>`, for "
        "example `120/45`."
    )


DATA = {
    "traj": ["struct_traj.h5ad"],
    "programs": ["kang2018.h5ad"],
    "mofa": ["struct_mofa"],
    "pbde": ["kang2018_pb.h5ad"],
}
CORPUS = {
    "traj": "skill: trajectory; corpus C-an-35 (sc.tl.paga on cluster groups), "
    "A-an-31 (fixed PAGA edge threshold), B-an-26 (seeded diffusion map)",
    "programs": "skill: gene-programs; corpus C-an-29 (NMF gene programs), D-an-22 "
    "(seeded NMF fit), D-an-23 (one-vs-rest AUC of program usage)",
    "mofa": "skill: mofa; corpus C-an-32 (MOFA2 on pseudobulk views)",
    "pbde": "skill: pseudobulk-de (python engine); corpus B-an-11 (sum per sample), "
    "B-an-14 (10%-style detection filter), B-an-13 (DESeq2 per cell type)",
}


def write_part(cache):
    """Compute every answer in a temp dir linking `cache`; write ../task_parts/struct.yaml."""
    import tempfile

    import yaml

    here, lines = Path(__file__).resolve().parent, []
    for task_id in SPEC:
        entry = {"id": task_id, "data": DATA[task_id]}
        for split in ("train", "valid", "test"):
            with tempfile.TemporaryDirectory() as work:
                cwd = Path.cwd()
                os.chdir(work)
                for name in DATA[task_id]:
                    os.symlink(Path(cache).resolve() / name, name)
                answer = ANSWER[task_id](SPEC[task_id][split])
                os.chdir(cwd)
            print(f"{task_id}/{split}: {answer}", flush=True)
            entry[split] = {"prompt": prompt(task_id, split)
                + (" Genome reference: none, synthetic genes." if task_id == "traj" else " Genome reference: GRCh37 (hg19), 10x Cell Ranger hg19 annotation."), "answer": answer}
        body = yaml.safe_dump([entry], sort_keys=False, width=92, allow_unicode=True)
        lines.append(f"  # {CORPUS[task_id]}\n" + "\n".join("  " + x for x in body.splitlines()))
    header = (
        "# struct group: trajectory, gene-programs, mofa, pseudobulk-de; generated by\n"
        "# `<target venv>/bin/python tasks/_structg.py <cache>` from tasks/_structg.py SPEC.\n"
    )
    (here.parent / "task_parts" / "struct.yaml").write_text(
        header + "tasks:\n" + "\n".join(lines) + "\n"
    )


def main(script):
    task_id, split = Path(script).stem.rsplit("-", 1)
    Path("answer.md").write_text(ANSWER[task_id](SPEC[task_id][split]) + "\n")


if __name__ == "__main__":
    import sys

    write_part(sys.argv[1])
