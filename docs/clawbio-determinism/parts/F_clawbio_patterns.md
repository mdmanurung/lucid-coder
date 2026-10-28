# F. ClawBio determinism patterns vs Mycelium / mycelium-extra

Paths: C = clawbio plugin 0.3.0 root. M = mycelium 0.7.2 root. X = mycelium-extra 0.16.1 root.
Line numbers are from the cached copies read on 2026-10-27.

Core finding: ClawBio moves every analysis and plot choice into one Python script per skill. The agent only dispatches and explains (`C/templates/SKILL-TEMPLATE.md:243-247`). Mycelium and mycelium-extra instead govern process, approval and provenance, and leave the analysis code and plot specs to the agent.

## 1. ClawBio determinism mechanisms

| ID | Mechanism | How ClawBio does it | Source |
|---|---|---|---|
| F-01 | Agent boundary | "The agent (LLM) dispatches and explains. The skill (Python) executes. The agent must NOT override thresholds or invent associations." | `C/templates/SKILL-TEMPLATE.md:243-247` |
| F-02 | Mandatory run path | Read SKILL.md, run the script with CLI args, show output. Skills with a script are never re-implemented by the agent. | `C/CLAUDE.md:96-102` |
| F-03 | Demo fallback | No input file means run `--demo` at once. Never refuse or ask for a file. | `C/CLAUDE.md:102` |
| F-04 | Bundled demo data | Each script resolves demo inputs from its own `examples/` or skill dir. Demo also fixes formula and contrast. | `C/skills/rnaseq-de/rnaseq_de.py:626-631`, `C/skills/sample-qc-triage/sample_qc_triage.py:177-181` |
| F-05 | Input contract table | SKILL.md lists format, extension, required columns, example per input. Detection by header. | `C/skills/diff-visualizer/SKILL.md:63-75`, `C/skills/rnaseq-de/SKILL.md:61-66`, `C/skills/bio-orchestrator/orchestrator.py:261` |
| F-06 | Validation with fixed fail strings | `ValueError` with exact messages: missing `sample_id`, non-numeric or negative counts, contrast level absent, groups under 2 samples, too few genes after filter. | `C/skills/rnaseq-de/rnaseq_de.py:53-138` |
| F-07 | Clean CLI failure | Bad input prints `ERROR: ...` to stderr, exit 2, no traceback. A test asserts it. | `C/skills/sample-qc-triage/sample_qc_triage.py:187-190`, `C/skills/sample-qc-triage/tests/test_sample_qc_triage.py:51-62` |
| F-08 | Fixed output tree | `report.md`, `figures/`, `tables/`, `reproducibility/`, plus `result.json`. File names are fixed in SKILL.md. | `C/skills/rnaseq-de/SKILL.md:68-85`, `C/skills/diff-visualizer/SKILL.md:129-156` |
| F-09 | Output contract test | The template says a `TestOutputContract` runs `--demo` and asserts every file in the SKILL.md tree exists. | `C/templates/SKILL-TEMPLATE.md:209-214` |
| F-10 | Fixed plot geometry | Hard-coded `figsize=(7, 5)` and `dpi=200` for PCA, volcano and MA. Other skills fix dpi 150, 160 or 180. | `C/skills/rnaseq-de/rnaseq_de.py:338,346,355,365,375,383`, `C/skills/scrna-orchestrator/scrna_orchestrator.py:963-991`, `C/skills/pathway-enricher/pathway_enricher.py:265` |
| F-11 | Fixed significance rule drawn on plot | Volcano marks `padj < 0.05` and `abs(log2FC) >= 1.0`. The MA plot marks `padj < 0.05`. The thresholds are baked in, not agent-chosen. | `C/skills/rnaseq-de/rnaseq_de.py:352-375` |
| F-12 | Threshold parameters with defaults | diffviz volcano takes `padj_threshold` and `lfc_threshold` as arguments. rnaseq-de CLI defaults: `--min-count 10`, `--min-samples 2`, `--formula "~ condition"`, `--contrast "condition,treated,control"`, `--backend auto`. | `C/skills/diff-visualizer/diff_visualizer.py:537-550`, `C/skills/rnaseq-de/rnaseq_de.py:609-613` |
| F-13 | Fixed style and palette | One `plt.rcParams.update` with named color dict (ink, muted, panel, spine). Heatmap `RdBu_r` with `vmin=-2.5, vmax=2.5`. Pathway heatmap `YlOrRd` 0 to 1. `viridis` elsewhere. `_PALETTE` list. | `C/skills/diff-visualizer/diff_visualizer.py:54-77,842`, `C/skills/pathway-enricher/pathway_enricher.py:192,233,365` |
| F-14 | Headless backend | `matplotlib.use("Agg")` in each plotting script. | `C/skills/diff-visualizer/diff_visualizer.py:23`, `C/skills/pathway-enricher/pathway_enricher.py:40`, `C/skills/scrna-orchestrator/scrna_orchestrator.py:946` |
| F-15 | Adaptive but formula-set sizes | Figure height is `max(min, k*n + c)` of the number of rows. Layout is a function of data, not of agent taste. | `C/skills/pathway-enricher/pathway_enricher.py:227`, `C/skills/diff-visualizer/diff_visualizer.py:616,899` |
| F-16 | Seed as CLI arg with default | `--random-state` default 0 is passed to PCA, UMAP and Leiden. Demo data built from `np.random.default_rng(random_state)`. Proteomics sets `np.random.seed(42)`. | `C/skills/scrna-orchestrator/scrna_orchestrator.py:94-98,518-521,1646`, `C/skills/proteomics-de/proteomics_de.py:175` |
| F-17 | Seed echoed in provenance | The seed is written back to `commands.sh` as a flag and into `result.json`. | `C/skills/scrna-orchestrator/scrna_orchestrator.py:1289,1504` |
| F-18 | commands.sh written by the script | The exact CLI that ran is written to `reproducibility/commands.sh`. Portable variant uses `ReproPath` anchors (`$CLAWBIO_ROOT`, `$OUTPUT_DIR`). | `C/skills/rnaseq-de/rnaseq_de.py:405-415`, `C/clawbio/common/reproducibility.py:17-60` |
| F-19 | environment.yml | Written per run. rnaseq-de writes unpinned names; the shared helper also offers `write_conda_lock` and RO-Crate. | `C/skills/rnaseq-de/rnaseq_de.py:417-428`, `C/clawbio/common/reproducibility.py:1-7` |
| F-20 | checksums.sha256 | SHA-256 of tables and figures, `sha256_file` shared helper. | `C/skills/rnaseq-de/rnaseq_de.py:430-437`, `C/clawbio/common/checksums.py:9` |
| F-21 | Shared report helpers | `generate_report_header`, `generate_report_footer`, `write_result_json`, `write_audit_log` give every report the same frame. | `C/clawbio/common/report.py:22-95` |
| F-22 | Disclaimer constant | One `DISCLAIMER` string embedded in report and result.json. A test asserts it is in `report.md`. | `C/skills/rnaseq-de/rnaseq_de.py:26-30,499-501,590`, `C/skills/sample-qc-triage/tests/test_sample_qc_triage.py:88` |
| F-23 | Thresholds documented, none hidden | Sample-QC rules are a numbered list in SKILL.md (e.g. `contamination_pct > 5`, `mapped_pct < 80`). Gotcha: "Do not use hidden thresholds". | `C/skills/sample-qc-triage/SKILL.md:125-133,168` |
| F-24 | Fixed-weight scoring | CRISPR triage uses fixed weights for depletion, druggability, essentiality. No external calls so demo and tests stay deterministic. | `C/skills/crispr-screen-triage/SKILL.md:85,105,168` |
| F-25 | Gotchas: "The model will want to ... Do not" | Names the specific agent failure and the reason (over-interpret p-values, submit thousands of genes, query unlisted databases). | `C/skills/pathway-enricher/SKILL.md:149-153` |
| F-26 | Closed option sets | Pathway enricher supports 6 hard-coded databases. `--backend` is a `choices` list. | `C/skills/pathway-enricher/SKILL.md:153`, `C/skills/rnaseq-de/rnaseq_de.py:611` |
| F-27 | Single-task scope | "One skill, one task": ranks hits, does not design guides or call statistics. Diffviz "recomputes no DE statistics". | `C/skills/crispr-screen-triage/SKILL.md:90`, `C/skills/diff-visualizer/SKILL.md:162` |
| F-28 | Demo/golden tests | Per-skill tests: parser unit tests, demo analysis asserts exact counts and labels (5 samples, 3 flagged, named dominant issue), CLI subprocess run, output files exist. | `C/skills/sample-qc-triage/tests/test_sample_qc_triage.py:65-90`, `C/skills/rnaseq-de/tests/test_rnaseq_de.py:30-135` |
| F-29 | Fixture beyond demo | Separate pseudo-bulk fixture dir for a second test path. | `C/skills/rnaseq-de/tests/fixtures/`, `C/skills/rnaseq-de/tests/test_rnaseq_de.py:100` |
| F-30 | trigger_keywords frontmatter | Lowercase phrases per skill. Template adds "Do NOT fire when" and says the trigger must be loud. | `C/skills/rnaseq-de/SKILL.md:36-45`, `C/templates/SKILL-TEMPLATE.md:76-86`, `C/skills/pathway-enricher/SKILL.md:47` |
| F-31 | Orchestrator routing | `KEYWORD_MAP` with longest-keyword-wins, file-type and tabular-header detection, chain-aware scRNA rule, optional LLM fallback only for ambiguous queries. | `C/skills/bio-orchestrator/orchestrator.py:60,242,261,314-365,386` |
| F-32 | Routing table in CLAUDE.md | Intent phrase to skill dir to action. Backed by catalog.json and optional `INTENTS.json` with schema, confined paths, demo only on explicit request. | `C/CLAUDE.md:26-80`, `C/docs/skill-intents.md:1-40`, `C/clawbio/skill_intents.py:25` |
| F-33 | Structured result handoff | Skill writes `chat_summary_lines`, `preferred_artifacts`, `suggested_actions`, `contract_alerts` to result.json. Runner promotes them. The agent does not invent next steps. | `C/docs/skill-action-contract.md:48-59` |
| F-34 | Negative case: repro-enforcer | Reproducibility export skill is only "Planned". No implementation. Reproducibility lives in each skill's script plus `clawbio/common/reproducibility.py`. | `C/skills/repro-enforcer/SKILL.md:49-51` |
| F-35 | Pre-fetched data for region plots | locuscompare uses a bundled config per demo, `INPUT_SCHEMA.md`, fixed 4-bin LocusZoom r2 palette, dpi 160, committed `commands.sh` and `environment.yml`. | `C/skills/locuscompare-region-render/cli.py:56-76,120-139`, `C/skills/locuscompare-region-render/regional_plot.py:38,752,864` |
| F-36 | SKILL-only skills are the weak case | When no script exists the agent "applies the methodology using your own capabilities". Determinism drops to prose guidance there. | `C/CLAUDE.md:104-108` |

## 2. Mycelium and mycelium-extra current state

| Mechanism | State | Source |
|---|---|---|
| Who writes analysis code | Agent-judgment. `analyze` routes to conventions and says "plain Python scripts for reproducible pipelines". It supplies no analysis script. | `M/skills/analyze/SKILL.md:31-39` |
| Reproduce entry point | Fixed requirement, free content. "Every analysis must have a `run.sh` or `run.py`." `new-analysis` scaffolds `run.sh` plus Snakefile. | `M/skills/analyze/SKILL.md:40`, `X/skills/new-analysis/templates/run.sh:1-9` |
| Analysis folder layout | Fixed. `analysis/<name>/` with doc, `scripts/`, `outputs/`, `reports/`. `new-analysis` adds numbered stubs, PLAN.md, TRACKER.md, flat `outputs/` with step-number prefix. | `M/skills/analyze/SKILL.md:29`, `X/skills/new-analysis/SKILL.md:10-16,80` |
| Output names inside outputs/ | Mostly agent-judgment. Numbered prefix is fixed. The rest is "descriptive names". `figures/` subfolders are optional and need a repo-local convention. | `X/skills/new-analysis/SKILL.md:80-82`, `M/skills/core/references/writing-conventions.md:31` |
| Input contract (required columns, fail messages) | Absent as a fixed contract. `data-contract-check` verifies a plan's assumptions against the sample table, read-only. Stub only says "Fail loudly: check shapes and counts". | `X/skills/data-contract-check/SKILL.md`, `X/skills/new-analysis/templates/stub.py:16` |
| Seeds | Agent-judgment with a nudge. Stub has `SEED = 0`. Core conventions say "set and document". `scilintr` flags unseeded RNGs. | `X/skills/new-analysis/templates/stub.py:13`, `M/skills/core/references/statistical-conventions.md:45-55`, `M/skills/analyze/SKILL.md:43` |
| Package versions / environment | Recorded after the fact. `verify` records conda and pip packages of runs. Core says "record versions". No per-run `environment.yml` writer. | `X/skills/verify/SKILL.md:3`, `M/skills/core/references/statistical-conventions.md:46` |
| Commands record | Fixed by the gate, not by the script. Pre/Post hooks write run receipts with exit status. Launch matching, Slurm `sacct` check. | `X/skills/verify/SKILL.md:41-46,69` |
| Input checksums | Fixed when the plan pins inputs. `Inputs:` line is hashed at approval. A changed pin blocks. No checksum of outputs unless the plan lists them. | `X/skills/grill/SKILL.md:68`, `X/skills/verify/SKILL.md:47-53` |
| Exact args frozen | Optional. `<!-- exact-execution -->` freezes chosen steps' arguments. `validation-bindings` pins tests. | `X/skills/grill/SKILL.md:71` |
| Plan before run | Fixed. grill produces a brief and plan table. Gate needs approval. verify checks what ran against the plan. | `X/skills/grill/SKILL.md:57-78`, `X/skills/verify/SKILL.md:6-8` |
| Numbers in prose | Fixed mechanism. `register_value` writes `outputs/numbers.json`. `<!-- claims -->` blocks are checked against output cells. | `M/skills/analyze/SKILL.md:41`, `X/skills/verify/SKILL.md:37` |
| Lint of analysis code | Fixed. `scilintr` after every edit, drive findings to zero or waive with a reasoned `ANALYSIS_OK`. Verify re-runs it. | `M/skills/analyze/SKILL.md:43-48`, `X/skills/verify/SKILL.md:33` |
| Which plot to make | Agent-judgment. plot-first: "name the plot, its output path, source data". The agent picks the plot type per question. | `X/skills/plot-first/SKILL.md:16-18` |
| Plot spec (size, dpi, palette, thresholds) | Absent from plot-first. Core says PDF preferred, min 300 dpi raster, colorblind palette, descriptive names. Bioinformatics pack says volcano uses padj 0.05, optional LFC cut, dashed lines at thresholds. Nothing sets figsize, style file, or exact palette. | `M/skills/core/references/writing-conventions.md:25-31`, `M/network/conventions/bioinformatics/statistical-conventions.md:14-15,77-83` |
| Plot inspection | Instruction only. "Open the rendered figure before interpreting it." The text says the gate does not enforce visual support. | `X/skills/plot-first/SKILL.md:36-39,76-78` |
| Conclusion to figure link | Fixed as a rule for the agent; unenforced. Verify counts `[agent-derived: file]` tags and claims blocks, not figure links. | `X/skills/plot-first/SKILL.md:44-62`, `X/skills/new-analysis/SKILL.md:71` |
| Table vs plot | Fixed preference: plot first, interpretive table needs a stated reason. | `X/skills/plot-first/SKILL.md:64-70` |
| Significance thresholds | Defaults documented (padj < 0.05, optional |log2FC| > 1) but chosen by the agent. A threshold in code passes lint only if the plan's Choice cell holds the number. | `M/network/conventions/bioinformatics/statistical-conventions.md:14-15`, `X/skills/verify/SKILL.md:33` |
| QC thresholds | Agent-judgment with documentation duty. "Document all QC thresholds in the analysis README. Common defaults (adjust per dataset)." | `M/network/conventions/bioinformatics/analysis-conventions.md:55-59` |
| Demo mode / bundled data | Absent. No analysis skill ships demo data. mycelium-extra has a labelled synthetic `verify exercise`, which tests the reviewer, not analysis. | `X/skills/verify/SKILL.md:10-12` |
| Golden / smoke tests for analysis output | Absent for analysis. Plugin tests cover Mycelium's own scripts. Validation bindings let a plan pin a test for its own run. | `X/skills/verify/SKILL.md:60` |
| Report | Router to a convention pack with planning brief, manifest, three blind reviewers and recompile gate. Wording is agent-written; numbers are bound by register_value. | `M/skills/report/SKILL.md:19-25` |
| Routing / triggers | Description-text triggers with "Do not use for" clauses, no keyword list or router. | `M/skills/analyze/SKILL.md:3-9`, `X/skills/plot-first/SKILL.md:3` |
| Disclaimer / safety text | Absent as a constant. verify emits a `disclose` paragraph on request. | `X/skills/verify/SKILL.md:3` |
| Sensitivity / null / adversarial checks | Described in robust-analysis packs (threshold sweeps, seed stability, method comparison), executed by the agent. | `M/network/conventions/robust-analysis/analysis-conventions.md:25-37` |

## 3. Gaps: mechanisms most worth porting

1. **F-01 + F-28 pattern: a tested, parameterised script per standard analysis, with the agent barred from re-implementing it.** This is the root of ClawBio's determinism and Mycelium has nothing equivalent for DE, PCA, QC or enrichment.
2. **F-10, F-11, F-13, F-14: one shared plot-spec module (rcParams, palette, figsize and dpi per plot type, threshold lines, Agg backend) used by plot-first.** plot-first says which evidence to show but not how it must look, so every run draws differently.
3. **F-05, F-06, F-07: input contract per analysis type with fixed fail strings and a clean exit code.** Today `data-contract-check` checks a plan, but nothing validates the data inside the analysis script itself.
4. **F-08, F-09: fixed output tree named in the skill, plus an output-contract test that runs on demo data.** `new-analysis` fixes the folder but not the file names, and nothing asserts the promised files exist.
5. **F-16, F-17, F-18, F-19, F-20: script-written provenance (seed flag, commands.sh, pinned environment, checksums of outputs) per run.** The gate records what ran, but the seed, the environment and the output hashes are not written beside the results.
6. **F-23, F-12, F-24: thresholds and defaults as named CLI parameters, listed in one table in the skill, with no hidden constants.** Mycelium lets the agent pick cutoffs and relies on lint plus the Choice cell to catch them afterwards.
7. **F-03, F-04, F-28: bundled `--demo` data and a golden-value smoke test.** It lets an agent or a CI job prove the pipeline still gives the same counts before touching real data.
8. **F-25, F-30, F-31: "The model will want to ... Do not" gotchas, "Do NOT fire when" triggers, and a keyword router.** Cheap prose changes that name the exact agent failure per analysis and cut mis-routing between analyze, plot-first and report.

Not worth porting: F-33 (chat-adapter action contract), F-21/F-22 medical disclaimer text, F-34 (repro-enforcer is unbuilt), F-36 (a SKILL-only skill is the weakness to avoid).
