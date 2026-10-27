# Purpose: lmerTest engine for longitudinal.py; never run by hand.
# Inputs:  long TSV (feature, subject, time, value, covariates) cleaned by longitudinal.py
# Outputs: <out_tsv> with feature, n_obs, n_subjects, estimate, se, df, t, pvalue, singular,
#          converged, fit_message; R, lme4 and lmerTest versions on stdout.
# Run:     Rscript lmer.R long.tsv "value ~ time + age + (1|subject)" time3 factor 1 0 out.tsv
# Notes:   exit 4 = an R package is missing, exit 3 = --term is not a fixed-effect coefficient;
#          longitudinal.py maps both to fixed error strings.

for (pkg in c(
  "lme4",
  "lmerTest",
  "reformulas",
  "readr",
  "dplyr",
  "purrr",
  "forcats"
)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("missing R package:", pkg, "\n", file = stderr())
    quit(status = 4)
  }
}
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

args <- commandArgs(trailingOnly = TRUE)
names(args) <- c(
  "input",
  "formula",
  "term",
  "time_type",
  "reference",
  "seed",
  "out"
)
# lmer draws no random numbers; the seed is set for parity with the Python record.
set.seed(as.integer(args[["seed"]]))

long <- readr::read_tsv(
  args[["input"]],
  col_types = readr::cols(
    feature = "c",
    subject = "c",
    time = "c",
    value = "d",
    .default = readr::col_guess()
  )
)
reference <- args[["reference"]]
# Factor time: the reference visit is the first level, so each time<level> is level minus
# reference. Numeric time is shifted so the intercept sits at the reference time.
long <- if (args[["time_type"]] == "factor") {
  mutate(
    long,
    time = forcats::fct(
      time,
      levels = c(reference, setdiff(sort(unique(time)), reference))
    )
  )
} else {
  mutate(long, time = as.numeric(time) - as.numeric(reference))
}
long <- mutate(
  long,
  across(where(is.character) & !c(feature, subject), \(x) {
    forcats::fct(x, levels = sort(unique(x)))
  })
)
model_formula <- as.formula(args[["formula"]])
term <- args[["term"]]

fixed_terms <- colnames(model.matrix(reformulas::nobars(model_formula), long))
if (!term %in% fixed_terms) {
  cat(
    "coefficients:",
    paste(fixed_terms, collapse = ", "),
    "\n",
    file = stderr()
  )
  quit(status = 3)
}

fit_feature <- \(d) {
  fit_messages <- character()
  # Singular-fit messages and convergence warnings are kept in fit_message and the
  # singular / converged columns, not printed.
  fit <- tryCatch(
    withCallingHandlers(
      lmerTest::lmer(model_formula, data = d),
      warning = \(w) {
        fit_messages <<- c(fit_messages, conditionMessage(w))
        invokeRestart("muffleWarning")
      },
      message = \(m) {
        fit_messages <<- c(fit_messages, trimws(conditionMessage(m)))
        invokeRestart("muffleMessage")
      }
    ),
    error = \(e) {
      fit_messages <<- c(fit_messages, conditionMessage(e))
      NULL
    }
  )
  row <- tibble::tibble(
    n_obs = nrow(d),
    n_subjects = n_distinct(d$subject),
    estimate = NA_real_,
    se = NA_real_,
    df = NA_real_,
    t = NA_real_,
    pvalue = NA_real_,
    singular = NA,
    converged = FALSE,
    fit_message = paste(unique(fit_messages), collapse = "; ")
  )
  if (is.null(fit)) {
    return(row)
  }
  coefs <- coef(summary(fit))
  row <- mutate(
    row,
    singular = lme4::isSingular(fit),
    # lme4 lists the singular-fit note among its convergence messages; only the
    # optimizer code and "failed to converge" count against convergence.
    converged = fit@optinfo$conv$opt == 0 &&
      !any(grepl(
        "failed to converge",
        fit@optinfo$conv$lme4$messages,
        fixed = TRUE
      ))
  )
  # A level missing from this feature drops the coefficient: the row stays NA.
  if (!term %in% rownames(coefs)) {
    return(row)
  }
  mutate(
    row,
    estimate = coefs[term, "Estimate"],
    se = coefs[term, "Std. Error"],
    df = coefs[term, "df"],
    t = coefs[term, "t value"],
    pvalue = coefs[term, "Pr(>|t|)"]
  )
}

results <- long |>
  split(long$feature) |>
  map(fit_feature) |>
  list_rbind(names_to = "feature")
readr::write_tsv(results, args[["out"]], na = "")
cat(R.version.string, "\n", sep = "")
cat("lme4 ", as.character(packageVersion("lme4")), "\n", sep = "")
cat("lmerTest ", as.character(packageVersion("lmerTest")), "\n", sep = "")
