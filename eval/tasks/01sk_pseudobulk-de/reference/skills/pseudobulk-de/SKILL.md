---
name: pseudobulk-de
description: >-
  Disease-vs-healthy differential expression from single-cell UMI counts, with the
  donor as the replication unit: sums counts per donor, filters genes with
  filterByExpr, and tests with edgeR quasi-likelihood. Writes de.tsv, a reviewer
  report, result.json and the exact command.
version: 0.1.0
license: MIT
metadata:
  domain: transcriptomics
  tags:
    - single-cell
    - pseudobulk
    - differential expression
    - edgeR
  inputs:
    - name: input
      type: directory
      description: holds counts.tsv (genes x cells) and cells.tsv (cell_id, donor_id, condition)
      required: true
  outputs:
    - name: de
      type: file
      format: [tsv]
      description: gene, log2_fc, p_value, q_value
    - name: report
      type: file
      format: [md]
      description: the procedure and top results, for a reviewer
  dependencies:
    r: ">=4.3"
    packages: [dplyr, readr, forcats, edgeR, limma]
  trigger_keywords:
    - pseudobulk
    - pseudobulk differential expression
    - disease vs healthy single-cell DE
    - donor-level DE
    - edgeR quasi-likelihood
---

# Pseudobulk DE

Run `skills/pseudobulk-de/pseudobulk_de.R`. Do not re-implement the analysis in
session, and do not change the model or cutoff without the user.

## Why This Exists

- **Without it**: cells get tested as if each were a replicate. With ~200 cells
  per donor, a donor-to-donor difference looks like a disease effect.
- **With it**: the donor is the unit. Counts are summed per donor and tested
  3 vs 3, so only effects that hold across donors are called.

## Core Capabilities

1. Sum raw UMI counts per donor (pseudobulk).
2. Drop low-count genes with edgeR `filterByExpr` and report how many.
3. Fit a quasi-likelihood negative binomial GLM, `~ condition`, healthy as reference.
4. BH-adjust over the tested genes; write a report a reviewer can check.

## Input Formats

| File | Format | Required content |
|---|---|---|
| `<input>/counts.tsv` | TSV, genes x cells | first column `gene`, one column per cell ID, raw UMI counts |
| `<input>/cells.tsv` | TSV, one row per cell | `cell_id`, `donor_id`, `condition` (`healthy` or `disease`) |

Each donor must belong to one condition. Every count column needs a row in `cells.tsv`.

## Workflow

1. Refuse an existing output directory.
2. Read both tables; align `cells.tsv` to the count columns by `cell_id`.
3. Sum counts per donor; count cells per donor.
4. `filterByExpr` on the donor sums; TMM normalisation.
5. `estimateDisp`, `glmQLFit`, `glmQLFTest` on `conditiondisease`.
6. BH over the tested genes; significant at q < 0.05.
7. Write `de.tsv`, `report.md`, `result.json`, `reproducibility/commands.sh`.

## CLI Reference

```bash
Rscript skills/pseudobulk-de/pseudobulk_de.R --input data --output results/run1
```

Both flags are required. The output directory must not exist.

## Demo

No bundled demo. The task data (6 donors, 1200 cells, 50 genes) is the worked
example: the run calls one gene, GENE33 (log2_fc about 2).

## Output Structure

```
<output>/
├── de.tsv                       # gene, log2_fc, p_value, q_value; sorted by p_value
├── report.md                    # data summary, procedure, top 10 genes, command
├── result.json                  # parameters, counts, n_significant
├── session_info.txt
└── reproducibility/
    └── commands.sh              # the exact command
```

## Dependencies

R >= 4.3 with dplyr, readr, forcats, edgeR and limma.

## Safety

- Never overwrites: an existing output directory stops the run.
- Local files only; no network access.
- Three donors per group gives low power. A gene not called is not shown to be unchanged.

## Integration with Bio Orchestrator

Route here for "which genes differ between disease and healthy" on single-cell
counts with donor labels. Hand `de.tsv` to an enrichment skill. Cell-level marker
questions go elsewhere.

## Citations

- Robinson, McCarthy, Smyth (2010). edgeR. Bioinformatics 26:139-140.
- Lun, Chen, Smyth (2016). It's DE-licious: quasi-likelihood DE with edgeR. Methods Mol Biol 1418:391-416.
- Squair et al. (2021). Confronting false discoveries in single-cell differential expression. Nat Commun 12:5692.
