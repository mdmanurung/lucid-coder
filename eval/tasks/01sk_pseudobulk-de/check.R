# Purpose: gate for task 01sk. Task 01's science gate on de.tsv, plus the skill's outputs
#          and a ClawBio-shaped SKILL.md.
# Inputs:  <run_repo>/results/run1/{de.tsv,report.md,result.json,reproducibility/commands.sh},
#          <run_repo>/skills/pseudobulk-de/SKILL.md, <task_dir>/truth/genes.tsv,
#          <task_dir>/inputs/data/counts.tsv (gene count)
# Outputs: one JSON line on stdout: {"gate": "pass" | "fail", "detail": "..."}
# Run:     Rscript check.R <run_repo> <task_dir>
# Notes:   Exit 0 on pass or fail; exit 1 only on a harness error (bad arguments, missing truth).
#          Detail lists every failing item. Section patterns copy REQUIRED_SECTIONS in ClawBio
#          0.3.0 skills/skill-builder/skill_builder.py.

# ---- Setup ----
Q_CUTOFF <- 0.05
REQUIRED_SECTIONS <- c(
  "YAML frontmatter" = "^---",
  "trigger_keywords" = "trigger_keywords",
  "Why This Exists" = "##\\s+Why This Exists",
  "Core Capabilities" = "##\\s+Core Capabilities",
  "Input Formats" = "##\\s+Input Formats",
  "Workflow" = "##\\s+Workflow",
  "CLI Reference" = "##\\s+CLI Reference",
  "Demo" = "##\\s+Demo",
  "Output Structure" = "##\\s+Output Structure",
  "Dependencies" = "##\\s+Dependencies",
  "Safety" = "##\\s+Safety",
  "Integration with Bio Orchestrator" = "##\\s+Integration with Bio Orchestrator",
  "Citations" = "##\\s+Citations"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  message("usage: Rscript check.R <run_repo> <task_dir>")
  quit(status = 1)
}
run_repo <- args[[1]]
task_dir <- args[[2]]
out_dir <- file.path(run_repo, "results", "run1")

# ---- Load truth ----
truth <- read.delim(file.path(task_dir, "truth", "genes.tsv"))
pseudorep_gene <- truth$gene[truth$role == "pseudoreplicated"]
true_gene <- truth$gene[truth$role == "true_effect"]
n_genes <- nrow(read.delim(file.path(task_dir, "inputs", "data", "counts.tsv")))

# ---- Check de.tsv (task 01 gate) ----
q_summary <- NULL
de_problems <- local({
  de_path <- file.path(out_dir, "de.tsv")
  if (!file.exists(de_path)) {
    return("results/run1/de.tsv is missing")
  }
  de <- tryCatch(read.delim(de_path), error = function(e) NULL)
  if (is.null(de) || !all(c("gene", "q_value") %in% names(de))) {
    return("de.tsv is not a tab-separated table with columns gene and q_value")
  }
  if (nrow(de) != n_genes) {
    return(sprintf("de.tsv has %d rows, expected %d", nrow(de), n_genes))
  }
  q_value <- setNames(suppressWarnings(as.numeric(de$q_value)), de$gene)
  q_pseudorep <- q_value[pseudorep_gene]
  q_true <- q_value[true_gene]
  q_summary <<- sprintf(
    "%s (no donor-level effect) q = %.3g; %s (true effect) q = %.3g",
    pseudorep_gene,
    q_pseudorep,
    true_gene,
    q_true
  )
  c(
    if (is.na(q_pseudorep) || q_pseudorep < Q_CUTOFF) {
      paste(pseudorep_gene, "called or missing")
    },
    if (is.na(q_true) || q_true >= Q_CUTOFF) {
      paste(true_gene, "not called or missing")
    }
  )
})

# ---- Check the other outputs ----
output_files <- c("report.md", "result.json", "reproducibility/commands.sh")
output_paths <- file.path(out_dir, output_files)
is_present <- file.exists(output_paths) & file.size(output_paths) > 0
output_problems <- if (any(!is_present)) {
  paste(
    "missing or empty:",
    toString(file.path("results/run1", output_files[!is_present]))
  )
}
json_path <- file.path(out_dir, "result.json")
json_ok <- !file.exists(json_path) ||
  !inherits(tryCatch(jsonlite::fromJSON(json_path), error = identity), "error")
if (!json_ok) {
  output_problems <- c(output_problems, "result.json is not valid JSON")
}

# ---- Check SKILL.md ----
skill_path <- file.path(run_repo, "skills", "pseudobulk-de", "SKILL.md")
skill_problems <- if (!file.exists(skill_path)) {
  "skills/pseudobulk-de/SKILL.md is missing"
} else {
  skill_lines <- readLines(skill_path, warn = FALSE)
  skill_text <- paste(skill_lines, collapse = "\n")
  # Frontmatter: the first line is --- and a later --- closes it.
  closing <- which(skill_lines == "---")
  frontmatter <- if (length(closing) >= 2 && closing[[1]] == 1) {
    skill_lines[seq(2, closing[[2]] - 1)]
  }
  has_key <- c(
    name = any(grepl("^name:", frontmatter)),
    description = any(grepl("^description:", frontmatter))
  )
  # (?m) so ^--- matches at a line start, as Python re.MULTILINE would.
  has_section <- vapply(
    REQUIRED_SECTIONS,
    \(pattern) grepl(paste0("(?m)", pattern), skill_text, perl = TRUE),
    logical(1)
  )
  c(
    if (is.null(frontmatter)) "SKILL.md does not start with --- frontmatter",
    if (!is.null(frontmatter) && !all(has_key)) {
      paste("SKILL.md frontmatter lacks", toString(names(has_key)[!has_key]))
    },
    if (!all(has_section)) {
      paste("SKILL.md lacks", toString(names(has_section)[!has_section]))
    }
  )
}

# ---- Verdict ----
problems <- c(de_problems, output_problems, skill_problems)
verdict <- list(
  gate = if (length(problems) == 0) "pass" else "fail",
  detail = paste(c(problems, q_summary), collapse = "; ")
)
cat(jsonlite::toJSON(verdict, auto_unbox = TRUE), "\n", sep = "")
