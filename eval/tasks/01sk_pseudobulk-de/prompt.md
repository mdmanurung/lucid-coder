I have single-cell RNA-seq data from six donors: three healthy and three with the disease, about 200 cells per donor, 50 genes.

- `data/counts.tsv`: raw UMI counts, genes × cells. The first column is `gene`; every other column is one cell, headed by its cell ID.
- `data/cells.tsv`: one row per cell, with `cell_id` (matches the count column names), `donor_id` and `condition` (`healthy` or `disease`).

Can you build a reusable skill that tests which genes differ between disease and healthy, in the ClawBio skill format, and run it on this data?

Write the skill at `skills/pseudobulk-de/`: a `SKILL.md`, one R script `pseudobulk_de.R` for data preparation and model fitting, and one Jupyter notebook `report.ipynb` (R kernel `ir`) for the figures and the report.

`SKILL.md` starts with YAML frontmatter holding `name: pseudobulk-de`, a `description`, and a `metadata` block that includes `trigger_keywords`. Its body has these `##` sections: Why This Exists, Core Capabilities, Input Formats, Workflow, CLI Reference, Demo, Output Structure, Dependencies, Safety, Integration with Bio Orchestrator, Citations.

The script runs as `Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output <out_dir>`. Both flags are required. It must refuse an `<out_dir>` that already exists. It writes into `<out_dir>`:

- `de.tsv` (tab-separated): one row per gene with the columns `gene`, `log2_fc` (disease vs healthy), `p_value` and `q_value` (Benjamini-Hochberg over all genes).
- `result.json`: the parameters, the counts, and the number of genes significant at q < 0.05.
- `reproducibility/commands.sh`: the exact commands that produce the outputs, including the two notebook commands below.
- any other tables the notebook needs.

The notebook reads only files in `<out_dir>`. Its first code cell is tagged `parameters` and sets `out_dir`. It renders `report.md`: the procedure, written for a reviewer. Give a data summary (counts of donors, cells and genes), the replication unit, how cells are aggregated, each filter and what it removed, the model or method and why, the contrast, the multiple-testing family, and the top results. Give any figure a caption.

Run it once:

```
Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output results/run1
papermill skills/pseudobulk-de/report.ipynb results/run1/report.ipynb -p out_dir results/run1
jupyter nbconvert --to markdown --no-input results/run1/report.ipynb
```

Use only base R and: dplyr, tidyr, readr, purrr, stringr, forcats, magrittr, ggplot2, here, limma, edgeR, jsonlite, knitr, IRdisplay.
