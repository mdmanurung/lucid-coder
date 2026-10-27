# The user's blinded picks, three rounds, 2026-10-27

The user is a computational biologist (immunology, single-cell, statistics). Each panel showed one analysis step in 2 or 4 versions; the user picked one and wrote why. Arms: A = baseline Opus, B = ponytail (lazy senior dev), C = lucid-coder rules, D = lucid-coder rules + 3 learned rules. Column letters have been unblinded to arms below. Human panels name the repo.

## Round 1: task 01sk (pseudobulk DE skill), four arms

| Step | Pick | Reason |
|---|---|---|
| setup | D | "I like Y because there is no unncesseray custom functions. but the name format parameters of W [A, UPPER_SNAKE] I like more, it makes it clear." |
| load | D | "to the point. but i do think the stopifnot call in Y [A] is important sanity check that I will also write anyway for my own analysis to ensure no silent error: stopifnot('every count column has a metadata row' = identical(cells[['cell_id']], colnames(counts_mat)))" |
| filter/aggregate | A | "to the point!" |
| fit | D | "lengthy but necessary to keep the whole gene lists in the dataset. code are sectioned well by linebreak smartly and appropriately." |
| write | B | "i like how you are saving the key parameters/arguments, but can be toned down a bit" |
| report notebook | B | "concise" |
| SKILL.md workflow | D | "key words first, detail second. make it easier for me to understand the major steps." |

## Round 1: human corpus panels (repo named after the pick)

| Step | Pick | Reason |
|---|---|---|
| R setup | immgenT-GP-analysis `script/Figure3.R:42-56` (library calls, then `source()` lines each with a comment naming the functions used) | "W is concise, but maybe add the #### section headers from Y [ALTRA]" |
| R load/join | dreamlet_analysis `Kfoury_CancerCell_2021_get_data.Rmd:105-124` (base R, each block preceded by a one-line comment) | "concise explanation of what code goes ahead by the comments" |
| R wrangle | sound-life `15-R_assemble_type_frequencies.ipynb` cell 23 (one long pipe, a comment above each verb) | "explained the code well" |
| R model many | ALTRA `RA_06-flow_statistical_modeling...Rmd:200-214` (group_by + group_modify, comments) | "adequate code comments, but used old pipe. prefer native R pipe |>" |
| R plot | none | "W/X/Z actually okay. though please do use native R pipe. Y is too dense, no whitespace" |
| R write | ALTRA `...Rmd:236-247` (#### export results #### banner, dir.create, write_csv per table) | "straightforward, clear section header, but a bit too crowded the header, simplify header" |
| R report | GTestimate `pub_figures.Rmd:28-49` | (no reason) |
| Py setup | cns-haematopoiesis `13_standard_reintegation_script.py:1-17` (imports, then `# ---- Defaults ----` block of UPPER_SNAKE parameters) | "goood organization. clear parameters and key variabels" |
| Py load/join | cns-haematopoiesis `05_basic_qc_filtering_add_metadata.py:30-48` (numbered step comments 1..6, map metadata, report missing) | (no reason) |
| Py model many | aifi `07-Python_B_cells.ipynb` cell 56 (plain for loop over tasks, sc.tl + sc.get) | (no reason) |
| Py plot | skin_fibroblast_atlas `extended_fig1.ipynb` cell 12 | (no reason) |
| Py notebook report | none | (no reason) |

## Round 2: task 01sk again after new rules (UPPER_SNAKE, step comments), four arms

| Step | Pick | Reason |
|---|---|---|
| setup | D | "header is informative and the code straightforward, and parameters are clear upfront." |
| load | D | "key checkpoints done, concise" |
| filter/aggregate | C | "use existing function instead of inventing new ones" |
| fit | D | "white space good, sectioning the code lines. comments adequate." |
| write | C | (no reason) |
| report notebook | D | "parameters up front, also settings and libraries" |
| SKILL.md workflow | A | (no reason) |

## Round 3: tasks 02-04, A vs D, like-for-like panels only

| Step | Pick | Reason |
|---|---|---|
| 02 whole script | D | (no reason) |
| 03 setup | D | (no reason) |
| 03 load | D | "comments informative and well positioned" |
| 03 figure | none | (no reason) |
| 04 whole script | none | "i like headers of X [A] better and the section dividers also X. but the snake case of parameters i like W [D, UPPER_SNAKE]. and also the functions I think are more concise in W for some of them. for writing table i prefer X because despite W is more 'efficient', but repeated write_output to me is much more intuitive and directly understandable. however i like how W use comments to explain the subsequent code. so need a mix" |

## Earlier explicit statements by the user (same project)
- "within a script, can you please not make any custom function if it is only used once or twice? if the function is used at least three times then do make the custom function... or if the function is used once per script but used a lot across many scripts then you can do that"
- "use R/tidyverse as much as possible instead of e.g. data.table or base-r"; "%>% should be superseded in favour of native R pipe |>"
- "i really want you to learn ... not in terms of the code formatting/linter, but in the brevity/conciseness of the code"; "not only the brevity, but more on the larger patterns of code writing"
- "Rscript for data prep/model fit --> ipynb for visualizations/reporting"; the writeLines report block at the end of a script is "very distracting"
- "i am done with the checks" (after the check budget was set: few checks, one line, able to fail)
- "procedurally correct, but i would not do this analysis, i would use specific DE analysis package such as limma/edgeR"
