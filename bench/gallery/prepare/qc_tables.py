"""Write qc_metrics.tsv: per-cell QC of three 10x PBMC libraries (qc_violin).

Kang 2018 ships no mitochondrial genes, so QC uses 10x PBMC data: pbmc3k (scanpy),
pbmc_1k_v3 and pbmc_10k_v3 (10x Genomics filtered matrices). Columns: sample, cell,
n_genes, total_counts, pct_mito (percent of counts on MT- genes). No cell is filtered.
"""

import shutil
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

URL = "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/{0}/{0}_filtered_feature_bc_matrix.h5"

out = Path(sys.argv[1])
tmp = out / ".download"
tmp.mkdir(exist_ok=True)  # pbmc3k() downloads here too
sc.settings.datasetdir = tmp
libraries = {"pbmc3k": sc.datasets.pbmc3k()}
for name in ("pbmc_1k_v3", "pbmc_10k_v3"):
    path = tmp / f"{name}.h5"
    # 10x answers 403 to urllib's default User-Agent.
    request = urllib.request.Request(URL.format(name), headers={"User-Agent": "curl/8"})
    with urllib.request.urlopen(request) as response:
        path.write_bytes(response.read())
    libraries[name] = sc.read_10x_h5(path)
shutil.rmtree(tmp)

rows = []
for name, adata in libraries.items():
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    rows.append(
        pd.DataFrame(
            {
                "sample": name,
                "cell": adata.obs_names,
                "n_genes": adata.obs["n_genes_by_counts"].to_numpy(),
                "total_counts": adata.obs["total_counts"].to_numpy().astype(int),
                "pct_mito": np.round(adata.obs["pct_counts_mt"].to_numpy(), 3),
            }
        )
    )
qc = pd.concat(rows)
qc.to_csv(out / "qc_metrics.tsv", sep="\t", index=False)
print(qc.groupby("sample").size())
