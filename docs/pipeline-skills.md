# Pipeline skill contract

A pipeline skill turns one standard analysis into a tested script that the agent runs instead
of writing its own code. It ports ClawBio's determinism framework
(`docs/clawbio-determinism/parts/F_clawbio_patterns.md`, F-01 to F-33).
`skills/pseudobulk-de/` is the worked example: copy its shape.

## Layout

```
skills/<name>/
  SKILL.md                       # contract, not a prompt; only what the agent needs
  references/maintainer.md       # corpus sources, ClawBio IDs, tests, citations, plot spec
  scripts/<name_underscored>.py  # the only entry point
  scripts/*.R                    # R engines, called by the .py via Rscript (optional)
  scripts/skillkit.py            # vendored, byte-identical to templates/skillkit.py
  scripts/plots.py theme.py palettes.py   # vendored, byte-identical to templates/
  examples/make_demo.py          # seeded generator of the demo data
  examples/demo.*                # committed demo data, under 1 MB
  tests/test_<name_underscored>.py
```

The folder is self-contained: acumen copies only this folder into the agent sandbox.

## Script rules

- **Shared helpers.** Use `skillkit.py` for `parse_required`, `fail` and `ERRORS`, `make_outdir`, `seed_everything`, `write_report`, `write_result`, `write_reproducibility` and `run`. Re-implement none of them.
- **No defaults** for any analysis parameter: thresholds, resolutions, gene sets, designs and seeds. Each is a required argument and a grill question. `--demo` alone pins demo-only values. Paths to bundled files are not analysis parameters.
- **Fixed failures.** Bad input prints `ERROR: <message>` and exits 2, with no traceback. Error strings live in one `ERRORS` dict, and each one has a test.
- **Fixed step order**, listed in SKILL.md. One skill does one task (ClawBio F-27); hand off to the next skill through its output tables.
- **Fixed output tree:**
  - `report.md` and `result.json` (summary fields plus `chat_summary_lines` and `preferred_artifacts`);
  - `tables/*.tsv`;
  - `data/<name>.h5ad` for an object the next skill reads. It is not checksummed, because h5ad bytes are not stable. Create `data/` in the script.
  - `figures/<name>.{pdf,png,csv}` through `plots.py`/`save_panel`. Write a form plots.py lacks with matplotlib in the same theme.
  - `reproducibility/{commands.sh, environment.txt, checksums.sha256}`.
- **Figures** carry no title, subtitle or caption. Each caption goes in `report.md`, opens with what is shown, and stays under 60 words (house-style `fig.text`). Set a layout engine on every figure, because `save_panel` does not crop.
- **Determinism.**
  - `--seed` is required and passed to every stochastic call.
  - Sort every output table by a stated key.
  - Pin joblib/BLAS threads to 1 where a library would fan out.
  - `SOURCE_DATE_EPOCH=0` keeps PDF checksums stable; `skillkit.run` sets it.
- **Offline demo.** No network at demo or test time. Bundle any reference the demo needs: a tiny model, gene sets or a network.
- **R engines.** R engines run through `Rscript`, found from `$RSCRIPT` or else PATH. Missing Rscript or a missing R package fails with a fixed message. R 4.5.2 here: `/share/software/statistical/R/R-4.5.2/gcc-8.5.0/R-4.5.2/bin/Rscript`, with `R_LIBS_USER=/exports/para-lipg-hpc/mdmanurung/R/4.5` and `LD_LIBRARY_PATH=/home/mdmanurung/.local/lib/conda-gcc-runtime`. Never install R packages: ask the parent session to.

## SKILL.md sections, in order

1. **Frontmatter:** `name`, and a `description` with trigger phrases and "Do NOT fire when".
2. **Agent boundary:** run the script; never re-implement it, change a threshold, or pick a value the user has not given; ask the user for each unknown argument (grill question); use `--demo` only on request.
3. **Input contract table:** file, format, required fields.
4. **Parameter table:** flag, meaning, grill question, demo-only value.
5. **Steps,** in order.
6. **Output tree.**
7. **"The model will want to … Do not",** with at least 4 gotchas.
8. **Maintainer notes:** one line pointing to `references/maintainer.md`.

`references/maintainer.md` holds what serves the maintainer, not the agent running the script:

1. **Corpus sources,** as `path:line` rows from `CATALOGUE.md` / `parts/*.md`.
2. **ClawBio mechanisms used,** with their F-IDs.
3. **Running and testing,** including the exact `uv run --with ...` line.
4. **Citations, fixed plot spec, reference skills** and any demo notes.

They moved out of SKILL.md on 2026-10-27: the acumen bench showed that a slimmer SKILL.md cut
the cost per run by about 20%, with the same accuracy and speed
(`docs/clawbio-determinism/BENCH.md`).

## Tests (pytest, CLI through subprocess)

- **Demo golden values:** exact counts and labels from `--demo`.
- **Output contract:** every path in the SKILL.md tree exists and is non-empty, and `sha256sum -c` passes.
- **Failure paths:** one test per `ERRORS` entry, checking exit 2 and no `Traceback`.
- **Missing argument:** a missing required argument without `--demo` fails.
- **Vendored files:** they equal `../../templates/` when that directory exists.

Run tests with
`UV_CACHE_DIR=/tmp/claude-149488/-exports-para-lipg-hpc-mdmanurung-lucid-coder/66eadf91-8195-442c-9410-99e83561a8a9/scratchpad/uvcache PYTHONDONTWRITEBYTECODE=1 uv run --python 3.12 --with <deps>,pytest python -m pytest -q -p no:cacheprovider tests`.
Format with `ruff format --config ../../templates/ruff.toml`, lint with `ruff check`, and format R with `air format`.

## Required checks (adopted 2026-10-27)

Adopted from the biology-skills v0.6.0 review
(`docs/clawbio-determinism/biology-skills-review-{stats,viz}.md`) after
the user approved each one. The `skillkit` helpers carry the fixed messages.

- **Unique identifiers.** Every skill that reads an h5ad or a table runs
  `skillkit.check_unique` on gene names and on sample or cell ids before any join.
- **Design rank.** Every skill that fits a design runs `skillkit.check_design_rank` on its
  fixed-effect terms. A rank-deficient design, such as batch fully confounded with
  condition, stops the run.
- **Donor versus sample.**
  - Where a donor can give several samples, a `--donor-col` is required and
    `skillkit.check_replication` runs.
  - Several samples per donor with no donor term in the design stops the run.
  - `result.json` reports `n_donors` and `n_samples`.
- **Cells as the unit.** Cells are the unit only when the condition was assigned to each
  cell (a per-cell perturbation). In every other case the donor is the unit; say so in the
  sc-de and sc-markers gotchas.
- **Reference record.** A gene-level skill takes a required `--reference` (genome build
  and annotation release, for example `GRCh38, GENCODE 44`) and passes it to
  `write_result(reference=...)`. It is a grill question and never a default.
- **Resource overlap.** A skill that joins a resource (gene sets, a network, a
  ligand-receptor list or a model) reports `n_matched` out of `n_in_resource` in
  `result.json`.
- **Drops per group.** QC and filter removals are reported per contrast group, not only as
  a total. The `result.json` key names what was dropped and the grouping, as
  `<what>_dropped_by_<group>` (`dropped_cells_by_level`, `genes_dropped_by_align`), since
  samples, cells, genes and rows are not comparable counts.
- **Captions say how.** A caption names the reduction (PCs, neighbours, `min_dist`, seed),
  any downsampling, and any transform. A heatmap caption names its order basis: clustered
  (metric and linkage) or fixed by label.
- **Plot arguments.** `plot_box(..., value_label=)` takes the y title with its unit, and
  `plot_embedding(..., axis_name=)` takes the embedding name, giving "UMAP 1" and "UMAP 2".
  Never pass a raw column name as a label. Under 10 points per group, `plot_box` draws
  points and a median, with no box.
- **Wording gotchas.** Add each one where it applies:
  - Non-detection is not absence: a filtered gene is "not tested".
  - A score is not a probability, and an enriched set is not a change in every member.
  - A bulk change is not a cell-intrinsic change: run composition next.
- **Longitudinal figures.** Never connect a line across a missing visit, and note
  dropped subjects. A mean line with a named interval (SD, SE or CI) may sit beside the
  per-subject lines.

## Biology-core checks (adopted 2026-10-28)

From the installed `biology-core` skill, approved by the user ("1/2/3/4/5 please").
The `skillkit` helpers carry the fixed behaviour.

1. **Context.** Every skill takes a required `--context` (species, tissue or cell
   type, condition, dataset; for example `human PBMC, IFN-beta vs control, Kang 2018`).
   It is a grill question. `write_report` prints it under the title, and
   `write_result(context=)` puts it first in `chat_summary_lines`. The demo value
   names the synthetic data.
2. **Status column.** Every per-feature result table (genes, gene sets, TFs, pairs)
   keeps every input feature, built with `skillkit.with_status`. Untested rows keep NA
   statistics and a `status` of `not_tested` (filtered), `not_detected`,
   `below_detection` or `not_measured`. `result.json` counts each status.
3. **Claim wording.** Report lines and captions say "associated with", "scored" or
   "predicted", never "drives" or "causes" (house-style `claims.strength`). Add a
   gotcha where a skill's output invites the stronger claim.
4. **Gene ids.** A gene-level skill writes `gene_id` beside `gene_symbol`:
   - h5ad input: a required `--gene-id-col` names the `var` column of Ensembl ids,
     and `var_names` are the symbols;
   - table input: the table must hold `gene_id` and `gene_symbol` columns;
   - `--gene-id-col none` is the explicit answer when the source has no ids (Kang
     2018 has symbols only); `gene_id` is then written as NA and `result.json`
     has `gene_id_col: none`, never a guessed mapping;
   - `skillkit.check_gene_ids` runs first and stops on an ambiguous map.
5. **Compatible inputs.** A skill that combines several inputs runs
   `skillkit.check_compatible` on the ids it joins them by: gene ids across
   mashr-sharing `--de` tables (mixed id types also stop), sample ids across mofa
   `--view` files (views share samples, not features). No shared ids stops the run;
   `result.json` reports `n_shared_ids`. group-tests reads one input, so it has none.

## Reference skills to check against

These installed skills document each tool's API. A pipeline skill's calls and gotchas
must agree with them; cite the one used in `references/maintainer.md`. `S` is
`~/.claude/skills/synced/f6b88190-5f7e-45de-8798-8d1f5b4ccea4_25160026-14a0-4c5c-af37-b6bb232e66d4`.

| Pipeline skill | Reference skill(s) |
|---|---|
| sc-qc | `S/single-cell-rna-qc`, `S/single-cell-best-practices` |
| sc-normalize, sc-embed | `S/single-cell-best-practices`, `single-cell-recipes`, `S/scvi-tools` (the later scVI engine) |
| sc-annotate | `S/celltype-annotation-workflow`, `S/celltypist-immune-markers`, `S/azimuth-pbmc-markers` |
| sc-markers | `scanpy-plotting`, `S/single-cell-best-practices` |
| pseudobulk-de (dreamlet), variance-partition, mashr-sharing | `dreamlet` |
| sc-de | `muscat` |
| composition | `crumblr-compositional` |
| milo-da | `milo-r`, `pertpy-tools`, `S/pertpy-perturbation-analysis` |
| group-tests, longitudinal | `S/glmm-practical-guide`, `S/emmeans-statistical-testing`, `S/marginaleffects-inference` |
| enrichment, signature-scoring, tf-activity | `decoupler`, `S/decoupler-enrichment` |
| gene-programs | `cnmf`, `hdwgcna` |
| ccc | `liana-ccc`, `S/liana-cell-communication` |
| trajectory, velocity | `cellrank` |
| mofa | `mofaflex` |

Skills without a pipeline yet: `partipy` (archetypes) and `complex-heatmap` (R heatmaps).
