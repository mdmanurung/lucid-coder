"""Shared spec and ground truth for the interp tasks (grouptest, signature, tfact, ccc).

SPEC is the single source: `python _interp.py <cache>` prints task_parts/interp.yaml
with answers, and every tasks/<id>-<split>.py is `import _interp; _interp.main(__file__)`.
Each procedure mirrors one pipeline skill (group-tests, signature-scoring, tf-activity,
ccc) but is written independently of its script, so a skill bug cannot hide in its key.
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DATA = "kang2018_lognorm.h5ad"
NET = "collectri_human.tsv"
RESOURCE = "liana_consensus.tsv"

PANEL = (
    "ISG15 IFI6 IFIT1 IFIT3 MX1 OAS1 IFITM3 LY6E STAT1 IRF7 CD3E CD14 LYZ MS4A1 NKG7 "
    "GNLY ACTB GAPDH B2M CCL5 S100A8 S100A9 FCGR3A HLA-DRA CD74 IL7R CCR7 FTL FTH1 "
    "TMSB4X"
).split()
SETS = {
    "ifn": "ISG15 IFI6 IFIT1 IFIT3 MX1 OAS1 IFITM3 LY6E IRF7 RSAD2 ISG20 IFI44L".split(),
    "mhc2": "HLA-DRA HLA-DRB1 HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 CD74".split(),
}

SPEC = {
    "grouptest": {
        "data": [DATA],
        "train": dict(cell_type="CD4 T cells", alpha=0.05),
        "valid": dict(cell_type="NK cells", alpha=0.1),
        "test": dict(cell_type="B cells", alpha=0.1),
    },
    "signature": {
        "data": [DATA],
        "train": dict(cell_type="B cells", set="ifn", ctrl_size=50, ask="highest"),
        "valid": dict(cell_type="B cells", set="mhc2", ctrl_size=50, ask="highest"),
        "test": dict(cell_type="CD14+ Monocytes", set="mhc2", ctrl_size=100, ask="highest"),
    },
    "tfact": {
        "data": [DATA, NET],
        "train": dict(label="ctrl", cell_type="CD8 T cells"),
        "valid": dict(label="ctrl", cell_type="CD14+ Monocytes"),
        "test": dict(label="ctrl", cell_type="NK cells"),
    },
    "ccc": {
        "data": [DATA, RESOURCE],
        "train": dict(label="ctrl", expr_prop=0.1, magnitude=0.01, specificity=0.05),
        "valid": dict(label="stim", expr_prop=0.1, magnitude=0.01, specificity=0.05),
        "test": dict(label="ctrl", expr_prop=0.2, magnitude=0.05, specificity=0.01),
    },
}

DATA_TEXT = (
    f"`{DATA}` in the working directory holds Kang et al. 2018 PBMCs: "
    "log-normalised values in `X` (`normalize_total` to 10,000 then `log1p`); `obs` "
    "columns `sample` (one id per patient x condition, e.g. `patient_101_ctrl`), "
    "`replicate` (patient), `label` (`ctrl` or `stim`, interferon-beta) and "
    "`cell_type`."
)
END = "Save the results tables and figures in the working directory."


def prompt(task_id, split_name):
    s = SPEC[task_id][split_name]
    if task_id == "grouptest":
        body = (
            f"Take the cells whose `cell_type` is {s['cell_type']}. For each `sample`, compute the mean "
            f"log-normalised value over its cells of each of these {len(PANEL)} genes: "
            f"{', '.join(PANEL)}. Test each gene stim vs ctrl with a two-sided paired "
            "Wilcoxon signed-rank test, pairing the two samples of each `replicate` "
            "(ctrl is the reference, so differences are stim minus ctrl; complete "
            "pairs only; zero differences dropped; exact p-values as scipy gives for "
            "8 pairs). Correct across the genes with Benjamini-Hochberg. A gene is "
            f"significant when the adjusted p-value < {s['alpha']:g}. {END} Report "
            "only the number of significant genes."
        )
    elif task_id == "signature":
        genes = SETS[s["set"]]
        body = (
            f"Take only the cells whose `cell_type` is {s['cell_type']} (subset first, "
            "then score). Score "
            f"the gene set `{s['set']}` ({', '.join(genes)}) in every cell with scanpy "
            f"`sc.tl.score_genes`: ctrl_size {s['ctrl_size']}, n_bins 25, "
            "random_state 0, use_raw False; require at least 5 of the set's genes in "
            "the data. Average the score per `sample` (samples nest in `label`, the "
            "groups). As a control, also score 20 random gene sets of the same size "
            "drawn from the other genes (seed 0) the same way. "
            f"{END} Report only the `sample` with the {s['ask']} mean score."
        )
    elif task_id == "tfact":
        body = (
            f"Take the `{s['label']}` cells. `{NET}` is the CollecTRI network "
            "(columns `source`, `target`, `weight`). Infer TF activity per cell with "
            "decoupler's univariate linear model (`dc.mt.ulm`, tmin 5). Then rank the "
            "sources per `cell_type` against the rest with `dc.tl.rankby_group` "
            "(reference `rest`, method `t-test_overestim_var`). "
            f"{END} Report only the source with the highest statistic in "
            f"{s['cell_type']}."
        )
    else:
        body = (
            f"Take the `{s['label']}` cells. `{RESOURCE}` is liana's consensus "
            "ligand-receptor resource (columns `ligand`, `receptor`). Run liana's "
            "rank aggregate of all five methods (CellPhoneDB, Connectome, log2FC, "
            "NATMI, SingleCellSignalR) between all `cell_type` pairs with this "
            f"resource: expr_prop {s['expr_prop']:g}, min_cells 5, n_perms 100, seed "
            "0, use_raw False, one job and one BLAS thread. An interaction passes when "
            f"magnitude_rank <= {s['magnitude']:g} and specificity_rank <= "
            f"{s['specificity']:g}. {END} Report only the number of passing "
            "interactions."
        )
    return f"{DATA_TEXT} {body}"


def grouptest(adata, s):
    from scipy import stats

    cells = adata[adata.obs["cell_type"] == s["cell_type"], PANEL]
    x = pd.DataFrame(cells.X.toarray(), columns=PANEL, dtype=np.float64)
    x["sample"] = cells.obs["sample"].astype(str).to_numpy()
    means = x.groupby("sample").mean()
    pairs = cells.obs[["sample", "replicate", "label"]].astype(str).drop_duplicates()
    wide = pairs.pivot(index="replicate", columns="label", values="sample").dropna()
    d = means.loc[wide["stim"]].to_numpy() - means.loc[wide["ctrl"]].to_numpy()
    p = [stats.wilcoxon(d[:, j], zero_method="wilcox").pvalue for j in range(len(PANEL))]
    padj = stats.false_discovery_control(p, method="bh")
    return str(int((padj < s["alpha"]).sum()))


def signature(adata, s):
    import scanpy as sc

    cells = adata[adata.obs["cell_type"] == s["cell_type"]].copy()
    genes = [g for g in SETS[s["set"]] if g in cells.var_names]
    assert len(genes) >= 5
    sc.tl.score_genes(
        cells, genes, ctrl_size=s["ctrl_size"], n_bins=25, random_state=0,
        use_raw=False, score_name="score",
    )
    means = cells.obs.groupby(cells.obs["sample"].astype(str))["score"].mean()
    return str(means.idxmax() if s["ask"] == "highest" else means.idxmin())


def tfact(adata, s):
    import decoupler as dc

    cells = adata[adata.obs["label"] == s["label"]].copy()
    cells.obs["cell_type"] = cells.obs["cell_type"].cat.remove_unused_categories()
    net = pd.read_csv(NET, sep="\t")
    dc.mt.ulm(cells, net, tmin=5)
    acts = dc.pp.get_obsm(cells, "score_ulm")
    ranked = dc.tl.rankby_group(
        acts, groupby="cell_type", reference="rest", method="t-test_overestim_var"
    )
    ranked = ranked[ranked["group"].astype(str) == s["cell_type"]]
    return str(ranked.sort_values(["stat", "name"], ascending=[False, True])["name"].iloc[0])


def ccc(adata, s):
    import liana as li
    from threadpoolctl import threadpool_limits

    cells = adata[adata.obs["label"] == s["label"]].copy()
    cells.obs["cell_type"] = cells.obs["cell_type"].cat.remove_unused_categories()
    resource = pd.read_csv(RESOURCE, sep="\t")
    with threadpool_limits(1):
        res = li.mt.rank_aggregate(
            cells, groupby="cell_type", resource=resource, expr_prop=s["expr_prop"],
            min_cells=5, n_perms=100, seed=0, n_jobs=1, use_raw=False,
            inplace=False, verbose=False,
        )
    passes = (res["magnitude_rank"] <= s["magnitude"]) & (
        res["specificity_rank"] <= s["specificity"]
    )
    return str(int(passes.sum()))


def answer(task_id, split_name):
    return globals()[task_id](ad.read_h5ad(DATA), SPEC[task_id][split_name])


def main(script):
    task_id, split_name = Path(script).stem.rsplit("-", 1)
    Path("answer.md").write_text(answer(task_id, split_name) + "\n")


if __name__ == "__main__":
    # Print the part file from recorded answers: python _interp.py answers.tsv
    import yaml

    rec = pd.read_csv(sys.argv[1], sep="\t", header=None, dtype=str)
    got = {(t, sp): a for t, sp, a in rec.itertuples(index=False)}
    rows = {
        "grouptest": "group-tests; corpus C-an-47 (paired Wilcoxon of scores pre/post), C-an-20, C-an-48",
        "signature": "signature-scoring; corpus A-an-14 (score_genes ctrl_size/n_bins, per-group means)",
        "tfact": "tf-activity; corpus A-an-27 (CollecTRI + ULM, rank_sources_groups vs rest)",
        "ccc": "ccc; corpus C-an-31 (liana consensus resource, expr_prop, n_perms)",
    }
    out = ["tasks:"]
    for task_id, spec in SPEC.items():
        entry = {"id": task_id, "data": spec["data"]}
        for sp in ("train", "valid", "test"):
            entry[sp] = {"prompt": prompt(task_id, sp) + " Genome reference: GRCh37 (hg19), 10x Cell Ranger hg19 annotation.", "answer": got[(task_id, sp)]}
        body = yaml.safe_dump([entry], sort_keys=False, width=88)
        out += [f"  # skill: {rows[task_id]}"]
        out += ["  " + line if line else line for line in body.splitlines()]
    print("\n".join(out))
