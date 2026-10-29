# Eval task catalog (2026-10-28)

Candidate tasks for the lucid-coder eval, merged from five sources. Each task gives the agent a small dataset and a prompt; `check.R` gates the result on a planted trap. Existing tasks (01, 01py, 01sk, 02, 02py, 03, 03py, 04, 05, 06) are not repeated.

Source reports, with every citation and loader: `eval/task-sources/`.
- `corpus-tasks.md`: the 17 cc_learnings repos (20 candidates, shipped-data list).
- `scbp-tasks.md`: sc-best-practices.org (15 candidates).
- `pertpy-tasks.md`: pertpy tutorials (12 candidates). Docs moved to pertpy.scverse.org; loaders are `pt.ds.*` in the tutorials and `pt.dt.*` in pertpy 1.0.2.
- `decoupler-liana-tasks.md`: decoupler v1.9.2 and latest, liana-py (10 candidates). Pin decoupler v2: liana >= 1.7.0 requires `decoupler>=2.0.0`.
- nf-core/scdownstream (read directly): MIT template scripts usable as gate references, and the `scdownstream` branch of nf-core/test-datasets, 99 files all under 1 MB, including 30 structural trap fixtures in `anndata-variations/` and three 384-cell PBMC samples in `pbmc/`.

Ranking: a trap a gate can check with one number or sign; a known answer from the source; small data available now; recurrence in the corpus. Sources agreeing independently count in favour.

## Build first

| # | Task | Trap the gate catches | Data | Known answer | Lang, packages | Sources |
|---|------|------------------------|------|--------------|----------------|---------|
| 1 | Repeated visits, mixed model | no donor random effect: visits treated as independent | real, shipped: MOCHA `SuppFig17_ZeroInflatedGLMMs/metadata.csv` (69 samples, 18 donors) + `tf_matrix.csv` (870 TF scores) | the repo's own per-TF lmer table (`tf_usage_model_results.csv`) | R: lme4, lmerTest | corpus #1 (4 repos); pertpy #3 (zhang_2021: 787 genes with cells as replicates, 0 with `(1 \| Patient)`) |
| 2 | Pathway activity from a DE contrast | scoring raw counts, log2FC or an unsigned statistic instead of the signed Wald stat | Kang 2018 CD14 monocyte pseudobulk DE table (2,345 genes) + PROGENy; or HSC TGF-beta (6 samples, 1.5 MB) | Kang: JAK-STAT 17.15, NFkB 2.78, both positive, only significant ones (ULM, deterministic); HSC: TGFb top | Python decoupler v2 (or R decoupleR) | scbp #2; decoupler #2; corpus #8 (sign loss, 3 repos) |
| 3 | Cell-type composition | per-type tests on proportions, or cells as the unit (n = 9,842 instead of 6 samples) | Haber 2017: ship the 10 samples x 8 types count table, derived once from `haber_2017_regions()` | only Enterocyte credibly up in Salmonella (scCODA, reference Goblet or Endocrine) | R crumblr, or Python scCODA (needs JAX) | scbp #1; pertpy #1, #9; corpus #5 and #2 (CLR with zeros, 4 repos) |
| 4 | Gene identifiers | Ensembl IDs or mixed symbols silently match few gene-set members, so most TFs or signatures drop with no error | HSC data indexed by Ensembl with symbols mixed in; or IFN `IFNa_IFNg_Scoring_ISGs_Mat.csv` (shipped, 272 KB); scdownstream `symbols_*`, `with_duplicates` fixtures | about 600 TFs scored (631 in v1) and SRF in the top 5 | Python decoupler v2; R for the IFN variant | decoupler #3; corpus #7 (3 repos); scdownstream |
| 5 | Contrast direction | alphabetical factor or column order picks the reference, so every sign flips | synthetic: 2 groups whose labels sort the wrong way | planted sign per gene | R base, limma | corpus #3 (4 repos) |
| 6 | Per-sample QC | one global cutoff wipes a low-depth sample; per-sample MAD keeps it | scdownstream `pbmc/` 3 x 384 cells, one sample thinned to plant the depth gap | planted cell counts kept per sample | Python scanpy or R scuttle | scbp #3 (ADT QC, per-donor 5 MAD); scdownstream `filter.py`; corpus #12 (doublets per sample, 3 repos) |
| 7 | Counts in a layer | normalising `X` that already holds log values, or reading counts from the wrong slot | scdownstream `counts_in_layer.h5ad` / `counts_in_X.h5ad` | known library sizes | Python anndata, scanpy | scdownstream |
| 8 | Interaction contrast | comparing two DE lists (significant in one cell type, not the other) instead of testing `cell_type:condition` | Kang 2018 pseudobulk, 128 samples, shipped as a table | RABGAP1L responds more in CD14 monocytes than CD4 T; ENO1 less | R edgeR or Python PyDESeq2 | scbp #4 |
| 9 | Regulator sign convention | a network without a `weight`/`mor` column scores every edge +1, so a repressor flips sign | decoupler toy data, 30 cells x 20 genes, offline | repressor T2 has mean score < 0 in group B | Python decoupler v2 | decoupler #1 |
| 10 | Visit and pair selection | keeping duplicate runs or the wrong visit, so pairs are unbalanced | synthetic: donors x visits with duplicates and one missing baseline | planted pairs | R dplyr | corpus #6 (3 repos), #14 |

## Built (2026-10-29)

All ten are built. Their gates run in `eval/tests/gates.test.js` (R) and `eval/tests/python.test.js` (Python): the reference passes and the naive answer fails under the harness itself. No arm has been run on them yet. Inputs are frozen; each task folder holds its own builder (`make-inputs.*` or `make-data.R`), not `eval/fixtures/make-fixtures.R`. Python tasks run in the env their `requirements.txt` names (`eval/envs/python-tasks.txt`, decoupler 2.3.1).

| # | Task id | Path | Gate status |
|---|---------|------|-------------|
| 1 | 07_repeated-visits | `eval/tasks/07_repeated-visits/` | pass: reference passes, naive fails (gates.test.js) |
| 2 | 08_pathway-activity | `eval/tasks/08_pathway-activity/` | pass: reference passes, naive fails, score.js scores (python.test.js) |
| 3 | 09_composition | `eval/tasks/09_composition/` | pass: reference passes, naive fails (gates.test.js) |
| 4 | 10_gene-ids | `eval/tasks/10_gene-ids/` | pass: reference passes, naive fails, score.js scores (python.test.js) |
| 5 | 11_contrast-direction | `eval/tasks/11_contrast-direction/` | pass: reference passes, naive fails (gates.test.js) |
| 6 | 12_per-sample-qc | `eval/tasks/12_per-sample-qc/` | pass: reference passes, naive fails (gates.test.js) |
| 7 | 13_counts-layer | `eval/tasks/13_counts-layer/` | pass: reference passes, naive fails, score.js scores (python.test.js) |
| 8 | 14_interaction-contrast | `eval/tasks/14_interaction-contrast/` | pass: reference passes, naive fails (gates.test.js) |
| 9 | 15_regulator-sign | `eval/tasks/15_regulator-sign/` | pass: reference passes, naive fails, score.js scores (python.test.js) |
| 10 | 16_visit-pairs | `eval/tasks/16_visit-pairs/` | pass: reference passes, naive fails (gates.test.js) |

Differences from the plan above:
- #2 (08): decoupler 2.3.1 with the current PROGENy gives JAK-STAT 16.73 on the Wald stat; the book reports 17.15 with an older PROGENy copy. The gate uses 16.73 (tolerance 0.5).
- #4 (10): the input is a DE table, so the raw-counts arm of the trap cannot occur; only the gene-ID trap is tested.
- 08, 10, 13 and 15 prompts run `uv run ... --with-requirements requirements.txt`, the same pinned lock the gate and rerun use. run.sh must copy `eval/envs/python-tasks.txt` to `requirements.txt` in each Python run repo (the task-folder file has a relative `-r` that does not resolve there).
- 07 ships `LICENSE-GPL-3.0.txt`; only that folder is GPL-3.0, the repository licence is undecided.
- Known soft spots, recorded in each `check.R`: 09 does not catch per-mouse proportion t-tests or per-mouse quasi-binomial GLMs (only 2 infected mice, too little power to call the depleted types); its up-type rule is Enterocyte log_fc > 0 and q < 0.25, so crumblr, crumblr robust, edgeR QL, t-tests, quasi-binomial and limma on logits all pass while the pooled cells-as-unit test fails; 12 and 13 do not check `qc_summary.tsv` or `combined.h5ad`, which their prompts ask for.

## Next tier

| Task | Trap | Data | Note |
|------|------|------|------|
| CLR axis for CITE-seq | CLR per protein (axis 0) instead of per cell | NeurIPS 2021 CITE-seq, 616 MB GEO; subset | scbp #5; big download |
| Guide RNA assignment | max-guide assignment drops negatives (2.9%) and multi-guide cells (8.7%) | `papalexi_2021()` gdo, 20,729 cells | pertpy #2; seconds on CPU |
| Milo design | site left out of the design, or the wrong sample column | `stephenson_2021_subsampled()`, 59,873 cells | pertpy #4; no count of significant neighbourhoods stated |
| Multiple-testing family | adjusting within each cell type vs across the family the claim spans | synthetic | corpus #4 (5 repos) |
| Age acceleration | clock minus age instead of the residual | synthetic, 80 people | corpus #17 |
| Group-aware cross-validation | folds split rows of one dataset | synthetic | corpus #15 |
| Mouse data, human resource | mouse symbols miss human CollecTRI | `dc.ds.erygast1k()`, 801 cells, v2 only | decoupler #4 |
| LIANA, too few cells | ranking interactions from a cell type with a handful of cells | `li.ds.kang_2018()` with a planted rare type | decoupler #9 |
| Pseudotime root | an arbitrary root cell | Setty 2019 bone marrow, 392 MB processed | scbp #6 |
| Dose response | unassigned cells used as the zero-dose reference | `srivatsan_2020_sciplex2()`, 145 MB | pertpy #5 |
| Negative control: cells per donor | vary cells per donor with donors fixed; the per-cell false-positive rate must not rise | synthetic, task 01 generator | design only; fold into task 01 later |
| Negative control: one subtype fixed | hold one subtype's absolute count fixed while the others grow; it must not be called changed in absolute terms | synthetic, known absolute counts | design only. Overlaps 09: both test closure. 09 has one real type rising and gates false calls on the rest; this is the mirror case with a known absolute truth. Build as a second fixture of 09, not a new task |

## Open before building

- Python env: pertpy 1.0.2 is in the conda env `pertpy`, but import did not finish in 120 s on the login node; scCODA needs JAX. decoupler, pydeseq2 and liana availability in a uv env is not yet checked.
- Downloads: figshare blocked size checks from this host. Build each fixture once, freeze it as a small table or h5ad under `eval/tasks/<id>/inputs/`, and never fetch at run time.
- Version pins: gate on the pinned version's numbers; decoupler v1 and v2 give different scores on the same data.
