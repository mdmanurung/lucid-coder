# biology-skills v0.6.0 vs our figure and identifier rules

Source: `ext/biology-skills/skills/{hadley-viz,bioinformatics,genomics}/`, paths below are relative to `skills/`.
Read in full: hadley-viz (SKILL.md plus 12 references), both other SKILL.md files, and the references on entity mapping, reference genomes, transcripts, gene expression and intervals. Skimmed by heading only: `qualifying-evidence.md`, `qualifying-variants.md`, `qc-outputs-and-sample-identity.md`, `metadata-and-provenance.md` and the file-format references (VCF, BAM, FASTA, FASTQ). They concern variant and MultiQC work we do not do.
No hadley-viz example plots or evals exist: `evals/*.jsonl` has 0 hadley entries. Every viz rule is prose with no cited figure, so none carries corpus evidence like our `fig.*` rows.
Our paths: `L` = `/exports/para-lipg-hpc/mdmanurung/lucid-coder`.
I ran no code and no colour simulator.

## 1. Adopt (ranked)

1. **Record annotation and reference identity per run.**
   - Theirs: "Genome assembly, patch release, exact reference sequence or FASTA, ... annotation release, transcript ... should not be collapsed" (`genomics/SKILL.md:12`); "Do not infer an assembly from chromosome naming alone" (`genomics/SKILL.md:34`).
   - Ours: no pipeline SKILL.md names a build or annotation release. A grep for `Ensembl|GRCh|hg38|mm10|GRCm` over `L/skills/*/SKILL.md` hits only `enrichment/SKILL.md:110` and `signature-scoring/SKILL.md:238`. I did not grep the scripts for these.
   - Put it in: `L/templates/skillkit.py` `write_result` (optional `reference` field copied from the h5ad `uns`) plus a `reproducibility/` line. It is a user question when absent, never a default (`L/docs/pipeline-skills.md`, "No defaults").
   - Why: a resource keyed to one build or species silently overlaps a data set from another. `ccc/SKILL.md:100` and `:234` already catch the species case by overlap only.
2. **Stop on duplicate or unmatched identifiers, and report the overlap.**
   - Theirs: "Treat unexpected row multiplication, row loss, duplicated identifiers ... as semantic events" (`bioinformatics/SKILL.md:26`); "compare row or entity counts before and after" (`bioinformatics/SKILL.md:56`).
   - Ours: only `enrichment.py:203` and `:217` check duplicates. sc-qc, sc-embed, sc-annotate and tf-activity have no `var_names` uniqueness check. `signature-scoring/SKILL.md:101` drops duplicate signature rows silently.
   - Put it in: an `ERRORS["duplicate_var_names"]` in each skill that reads an h5ad. Add "n genes matched / n in resource" to `result.json` for tf-activity, ccc and sc-annotate.
   - Why: a 12-skill chain is only as sound as its first join. Only `signature-scoring/SKILL.md:148` records missing genes per set.
3. **Caption carries reduction and transform facts.**
   - Theirs: "Document reduction: ... 'points downsampled by N%'" (`hadley-viz/references/aggregation-sampling-and-large-data.md:30`); "Record the random seed or parameters" (`.../high-dimensional-embeddings.md:20`); "UMAP of first 30 PCs" (`:54`).
   - Ours: the sc-embed caption template has no PC count, `n_neighbors`, `min_dist` or seed. Check `sc-embed/SKILL.md:205` (report) against the script.
   - Put it in: sc-embed and sc-markers caption strings. They stay under the 60-word `fig.text` cap.
   - Why: a reader of a UMAP panel cannot reproduce it from the caption.
4. **Define "few points" numerically in `fig.wilke`.**
   - Theirs: "If n<20, prefer dot plots"; "violin plot on n=3 yields a misleading continuous shape" (`hadley-viz/references/distributions-summaries-and-sample-size.md:16`, `:38`). Also "n<10" in `hadley-viz/SKILL.md:42`.
   - Ours: "no violin for few points" (`fig.wilke`, `L/skills/lucid-coder/references/house-style.md:126`). It gives no number. `tf-activity/scripts/tf_activity.py:389` draws violins per cell; that is fine because each group has many cells.
   - Put it in: `house-style.md` `fig.wilke` and the `plot_box` guard in `L/templates/plots.{py,R}`.
   - Why: their cut-offs conflict (10 vs 20), so adopt the principle and pick one number with the user. Neither value is evidence-backed.
5. **Name the basis of heatmap order in the caption.**
   - Theirs: "Reordered axes: Presenting clustered heatmap as if original order (no indication of reordering)" (`.../heatmaps-matrices-and-clustering.md:37`); state distance and linkage when clustering (`:17`, `:24`).
   - Ours: row and column order are the caller's factor levels (`L/templates/plots.R:92`). The caption does not say what the order means.
   - Put it in: `fig.heatmap`: "caption states the order basis (clustered with metric and linkage, or fixed by label)". tf-activity's caption at `tf_activity.py:385` already states the saturation limit.
   - Why: it costs one clause and closes the "adjacency is algorithmic" misreading (`:20`).
6. **Do not connect across missing time points.**
   - Theirs: "Do not connect lines across a missing point" (`.../missingness-censoring-and-detection-limits.md:17`); "do not drop [subjects] without note" (`.../longitudinal-paired-and-repeated-measures.md:31`).
   - Ours: no `fig.*` row covers gaps. `longitudinal/SKILL.md` was not read for this.
   - Put it in: longitudinal figure code, once its trajectory plot is checked for this.
7. **Bulk change is not cell-intrinsic change.**
   - Theirs: "A change in bulk RNA abundance can arise because cell-type proportions changed" (`genomics/references/gene-expression.md:59`).
   - Ours: pseudobulk-de and composition are separate skills. No sentence ties a DE result to a composition check.
   - Put it in: the pseudobulk-de "model will want to" section as one gotcha, with composition named as the next skill.
8. **Never resolve an ambiguous mapping by "first row".**
   - Theirs: "Do not write code that resolves an ambiguous biological mapping by selecting the first returned row" (`bioinformatics/references/entity-mapping-and-join-cardinality.md:58`).
   - Ours: no script maps IDs today, so there is nothing to fix. Add as a gotcha if a symbol-to-Ensembl step is ever added to any skill.

## 2. Conflicts

1. **Figure subtitle and in-plot n.** Ours is better supported.
   - Theirs: "Report the number of observations ... (e.g. 'n = 8' in subtitle or caption)" (`hadley-viz/references/distributions-summaries-and-sample-size.md:21`); "subtitle 'Region on chr7 (GRCh38) ...'" (`.../genomic-coordinates-variants-and-effects.md:55`); "always output figure captions or metadata strings" (`hadley-viz/SKILL.md:41`).
   - Ours: "no subtitle and no caption inside the plot; a title only when the panel stands alone" (`house-style.md:125`; `L/templates/plots.py:5`). n goes in tick labels (`fig.labels`, `:129`).
   - Evidence: `fig.text` cites the immgenT published figures, which carry no subtitle. Their rule has no cited figure.
   - Reconcile: their content requirements (build, transcript, transform, sampling) are right. Put them in the `report.md` caption, not in the plot.
2. **UMAP: shapes and axis names.** Split verdict.
   - Theirs: "UMAP with cells colored by cell type ... and shapes for patient ID" (`hadley-viz/SKILL.md:72`); axes "UMAP1", "UMAP2" (`.../high-dimensional-embeddings.md:54`); "Do not label axes with original variable names" (`:16`).
   - Ours: colour-only, `s=1`, rasterised, ticks removed, then `ax.set_xlabel(x)` where `x` is the raw column `umap_1` (`L/templates/plots.py:98-111`, `fig.raster`, `house-style.md:123`). Same at `plots.R:75`.
   - Shapes: ours is better. A shape per donor on 10^4 to 10^5 points at `s=1` is unreadable, and `fig.wilke` limits "colour plus shape" to categorical panels.
   - Labels: theirs is better. `umap_1` as an axis title breaks our own `fig.labels` ("never a raw column or variable name", `:129`). `plot_box` has the same issue (`plots.py:93-94`).
   - Fix: write "UMAP 1" and "UMAP 2" in `plot_embedding`.
3. **Error bars.** Mostly compatible; the real gap is size.
   - Theirs: "Plot shows points and error bars labeled 'SD'" (`.../uncertainty-estimates-and-statistical-inference.md:49`); "summarizing with mean±error over time" (`.../longitudinal-paired-and-repeated-measures.md:18`).
   - Ours: "each donor a point over any summary, never bar plus error bar" (`house-style.md:126`).
   - Verdict: theirs allows points plus an interval, which ours also allows; only "bars" are banned. For many subjects, ours has no rule for a mean±error line. Theirs would permit it; flag it as an open question.
   - Their "SD for variability, SE/CI for inference" (`:17`) matches our "name SD, SE or CI in the caption".
4. **Significance stars.** Compatible, ours stricter.
   - Theirs: "do not use '*' ... unless a statistical test was done and described" (`.../uncertainty-estimates-and-statistical-inference.md:20`).
   - Ours: thresholds on adjusted p, omit "ns" (`fig.stars`, `house-style.md:128`).
   - Reconcile: add "caption names the test and the adjustment".
5. **Log scales and zeros.** Ours is better supported.
   - Theirs: "Y-axis label 'log10(Cell count + 1)' with caption 'Zero values shown as 0 on plot (log scale)'" (`.../transformations-scales-and-zero.md:50`); suggest `log(x+1)` (`:30`); bars "start at zero" (`:18`).
   - Ours: "bars start at 0, log-scale bars at 1, otherwise dots" (`house-style.md:126`).
   - Evidence: a zero cannot be drawn on a log axis, so their example contradicts itself. Their baseline rule gives no case for log bars.
   - Gap in ours: `sc-qc` plots UMI counts on a log scale (`sc_qc.py:502-505`) with no stated zero rule. Not verified whether a zero count can reach that plot.
   - Adopt only "state the log base and the zero treatment in the caption" (`:19`, `:26`).
6. **Colour: grayscale test and viridis.** Ours is better supported on categories, theirs adds a check.
   - Theirs: "Use palettes where every category can be differentiated in grayscale" (`.../visual-encoding-colour-and-figure-integrity.md:16`); "WCAG ... >=4.5:1" for text and lines (`:17`).
   - Ours: Tol-style label-keyed maps (`L/templates/palettes.py:13-14`; sc-embed uses a 10-colour list, `sc-embed/SKILL.md:10`); `fig.colourblind` asks for a simulator but has "0 corpus witnesses" (`house-style.md:130`).
   - Verdict: a 10-colour list cannot pass their grayscale rule, so the rule is not workable for clusters. Keep ours. The gap is that we require a simulator but ship none.
7. **Diverging centre.** No conflict; recorded because the brief asked.
   - Theirs: "two contrasting hues with a white midpoint" (`.../visual-encoding-colour-and-figure-integrity.md:16`); "midpoint indicated" in the legend (`.../heatmaps-matrices-and-clustering.md:18`).
   - Ours: `#2166AC` / white / `#B2182B`, `midpoint = 0`, symmetric limits, squish (`plots.R:98-104`; `plots.py:25-26`, `:124`).
   - Their extra: label the midpoint on the colour key. Ours does not label it.
8. **Gene IDs: symbols vs stable IDs.** Theirs is the better general rule; ours is forced by the resources.
   - Theirs: "A human-readable gene symbol ... may change over time. Stable identifiers and mapping provenance are safer for machine joins" (`bioinformatics/references/entity-mapping-and-join-cardinality.md:101`); "Keep stable machine identifiers separate from human-readable display labels" (`bioinformatics/SKILL.md:32`); "Preserve the original source identifier and the mapping resource or release" (`...join-cardinality.md:26`).
   - Ours: symbols matched by exact string (`signature-scoring/SKILL.md:101`; `tf-activity/SKILL.md:106`; `sc-annotate/SKILL.md:148`; `sc-qc/SKILL.md:136`). `enrichment/SKILL.md:110` allows either type but requires the same in both files.
   - Evidence: CollecTRI and PROGENy are symbol-keyed (`tf-activity/SKILL.md:106`), so switching to Ensembl would need a mapping step with a cardinality policy (`...join-cardinality.md:20`).
   - Reconcile: keep symbols; adopt Adopt items 1 and 2 so the symbol join leaves a trace.
   - Also: `ccc/SKILL.md:234` bans uppercasing `var_names` to force a match, which agrees with their no-silent-mapping rule. `sc-qc.py:207` uppercases `var_names` for the ribosomal prefix test only.
9. **Heatmap clustering.** No conflict; theirs is a caution, ours a constraint.
   - Theirs: "Applying hierarchical clustering heatmaps without biological rationale can suggest patterns that are methodological" (`hadley-viz/SKILL.md:51`); cluster "a subset of genes or ... heatmap of principal components" for 500 x 20,000 (`.../aggregation-sampling-and-large-data.md:52`).
   - Ours: fixed caller order (`plots.R:92`).
   - Where they differ: their heatmap-of-PCs advice would break our named-gene row rule. Neither side has evidence; keep ours.

## 3. Already covered

- Experimental unit is the donor, not the cell: theirs `hadley-viz/references/experimental-units-and-visual-independence.md:16`; ours `fig.wilke` ("each donor a point", `house-style.md:126`) and `plot-first/SKILL.md` ("the replication unit").
- Paired data connected: theirs `.../experimental-units-and-visual-independence.md:18`; ours `fig.wilke` slopegraph / x = y and `plot_box(pair=)` (`plots.py:69`).
- No violin on few points: theirs `.../distributions-summaries-and-sample-size.md:38`; ours `fig.wilke` (`house-style.md:126`).
- Name the interval type: theirs `.../uncertainty-estimates-and-statistical-inference.md:16`; ours `fig.wilke` "name SD, SE or CI in the caption" and `plot-first/SKILL.md` ("interval types").
- Raw points over bars, n shown: theirs `.../distributions-summaries-and-sample-size.md:16`, `:21`; ours `fig.wilke` and `fig.labels` `(n = 12)` (`house-style.md:129`).
- Overplotting: theirs `.../aggregation-sampling-and-large-data.md:16`; ours `fig.wilke` "smaller points, then alpha, then hex" and `fig.raster` (`house-style.md:123`).
- Diverging heatmap centred, legend: theirs `.../heatmaps-matrices-and-clustering.md:18`; ours `fig.heatmap` (`house-style.md:127`).
- Heatmap scaling stated: theirs `.../heatmaps-matrices-and-clustering.md:16`; ours `fig.wilke` "row-scaled heatmap, stated".
- Colour not the only cue: theirs `.../visual-encoding-colour-and-figure-integrity.md:18`; ours `fig.wilke` "colour plus shape, no red-green".
- Direct labels past a category count: theirs `.../visual-encoding-colour-and-figure-integrity.md:20`; ours `fig.wilke` `legend_loc="on data"`.
- Batch visible in an embedding: theirs `.../high-dimensional-embeddings.md:37`; ours `umap_batch` (`sc-embed/SKILL.md:216`, `:239`).
- Do not read cell types from the UMAP: theirs `.../high-dimensional-embeddings.md:18`; ours `sc-embed/SKILL.md:262` ("name clusters as cell types from the UMAP. Do not.").
- Fixed seed for the embedding: theirs `.../high-dimensional-embeddings.md:20`; ours `sc-embed/SKILL.md:166` (`--seed`) and `L/docs/pipeline-skills.md` ("Determinism").
- Common denominator and axis label for composition: theirs `.../categorical-compositional-and-denominator-data.md:17`; ours `plot_composition` y label `% of {count}` (`plots.py:170-173`).
- Do not assume an ID match is correct; record missing genes: theirs `...join-cardinality.md:24`; ours `signature-scoring/SKILL.md:148` (record missing genes per set, stop below `--min-genes`).
- Do not silently remap species: theirs `...join-cardinality.md:23`; ours `ccc/SKILL.md:234`.
- Keep planned and observed values apart; keep source and aggregate apart: theirs `bioinformatics/SKILL.md:30`; ours `L/docs/pipeline-skills.md` (the `data/` h5ad is the handoff, not a checksummed result) and per-skill output trees.
- Absence of detection is not absence: theirs `genomics/references/gene-expression.md:22`; ours `plot-first/SKILL.md` ("Show null, mixed, and adverse results").
- Save, inspect and cite a figure: theirs `hadley-viz/SKILL.md:41` (document figure metadata); ours `plot-first/SKILL.md` (inspect before interpreting; cite path, groups, n, axes, interval).
