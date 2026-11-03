# Purpose: lmerTest engine for longitudinal.py; never run by hand.
# Inputs:  long TSV (feature, subject, time, value, covariates) cleaned by longitudinal.py
# Outputs: <out_tsv> with feature, n_obs, n_subjects, estimate, se, df, t, num_df, f, pvalue,
#          singular, converged, fit_message; R, lme4, lmerTest (and mice) versions on stdout.
# Run:     Rscript lmer.R long.tsv "value ~ time + age + (1|subject)" time3 factor 1 0 out.tsv \
#            "" "" complete-case
# Notes:   exit 4 = an R package is missing, exit 3 = --term is not a fixed-effect coefficient;
#          longitudinal.py maps both to fixed error strings. Spline time: `time` becomes the
#          natural-spline basis (columns time1..timeK), and --term time is their joint F test.
#          Impute: mice per feature on a subject-level table, the fit on each imputation, the
#          --term coefficient pooled with Rubin's rules.

# Load packages ----------------------------------------------------------------
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
  "out",
  "knots",
  "boundary",
  "missing"
)
seed <- as.integer(args[["seed"]])
# lmer draws no random numbers; mice does and takes the seed per feature.
set.seed(seed)
impute <- args[["missing"]] == "impute"
if (impute && !requireNamespace("mice", quietly = TRUE)) {
  cat("missing R package: mice \n", file = stderr())
  quit(status = 4)
}

# Load data --------------------------------------------------------------------
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
# reference. Numeric time is shifted so the intercept sits at the reference time. Spline
# time is the ns() basis on knots fixed across features (longitudinal.py chose them), so
# every feature is fitted on the same curve family.
long <- switch(
  args[["time_type"]],
  factor = mutate(
    long,
    time = forcats::fct(
      time,
      levels = c(reference, setdiff(sort(unique(time)), reference))
    )
  ),
  numeric = mutate(long, time = as.numeric(time) - as.numeric(reference)),
  spline = mutate(
    long,
    time = matrix(
      splines::ns(
        as.numeric(time),
        knots = as.numeric(strsplit(args[["knots"]], ",")[[1]]),
        Boundary.knots = as.numeric(strsplit(args[["boundary"]], ",")[[1]])
      ),
      nrow = n()
    )
  )
)
long <- mutate(
  long,
  across(where(is.character) & !c(feature, subject), \(x) {
    forcats::fct(x, levels = sort(unique(x)))
  })
)
model_formula <- as.formula(args[["formula"]])
term <- args[["term"]]
# Covariates with a missing value: only present under --missing impute.
imputed_cols <- names(long)[map_lgl(long, anyNA)]

# Check the term ---------------------------------------------------------------
fixed_terms <- colnames(model.matrix(reformulas::nobars(model_formula), long))
spline_cols <- grep("^time[0-9]+$", fixed_terms, value = TRUE)
joint <- args[["time_type"]] == "spline" && term == "time"
if (!term %in% fixed_terms && !joint) {
  cat(
    "coefficients:",
    paste(
      c(fixed_terms, if (args[["time_type"]] == "spline") "time"),
      collapse = ", "
    ),
    "\n",
    file = stderr()
  )
  quit(status = 3)
}

# Fit helpers ------------------------------------------------------------------
# One lmer fit; singular-fit messages and convergence warnings are kept, not printed.
fit_lmer <- \(d) {
  fit_messages <- character()
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
  list(fit = fit, messages = fit_messages)
}

# The --term test of one fit: a Satterthwaite t-test of one coefficient, or the joint
# Satterthwaite F test of the spline columns. NA when the coefficient was dropped.
na_test <- tibble::tibble(
  estimate = NA_real_,
  se = NA_real_,
  df = NA_real_,
  t = NA_real_,
  num_df = NA_real_,
  f = NA_real_,
  pvalue = NA_real_
)
test_term <- \(fit) {
  row <- na_test
  beta <- lme4::fixef(fit)
  if (joint) {
    if (!all(spline_cols %in% names(beta))) {
      return(row)
    }
    contrast <- diag(length(beta))[
      match(spline_cols, names(beta)),
      ,
      drop = FALSE
    ]
    ftest <- lmerTest::contest(fit, contrast, joint = TRUE)
    return(mutate(
      row,
      df = ftest[["DenDF"]],
      num_df = ftest[["NumDF"]],
      f = ftest[["F value"]],
      pvalue = ftest[["Pr(>F)"]]
    ))
  }
  if (!term %in% names(beta)) {
    return(row)
  }
  coefs <- coef(summary(fit))
  mutate(
    row,
    estimate = coefs[term, "Estimate"],
    se = coefs[term, "Std. Error"],
    df = coefs[term, "df"],
    t = coefs[term, "t value"],
    pvalue = coefs[term, "Pr(>|t|)"]
  )
}

# lme4 lists the singular-fit note among its convergence messages; only the optimizer
# code and "failed to converge" count against convergence.
is_converged <- \(fit) {
  fit@optinfo$conv$opt == 0 &&
    !any(grepl(
      "failed to converge",
      fit@optinfo$conv$lme4$messages,
      fixed = TRUE
    ))
}

# Fit each feature -------------------------------------------------------------
# Complete-case: one fit per feature.
fit_feature <- \(d) {
  one <- fit_lmer(d)
  row <- tibble::tibble(
    n_obs = nrow(d),
    n_subjects = n_distinct(d$subject),
    singular = NA,
    converged = FALSE,
    fit_message = paste(unique(one$messages), collapse = "; ")
  )
  if (is.null(one$fit)) {
    return(bind_cols(row, na_test))
  }
  row |>
    mutate(
      singular = lme4::isSingular(one$fit),
      converged = is_converged(one$fit)
    ) |>
    bind_cols(test_term(one$fit))
}

# Impute: mice on one row per subject (the imputed covariates and the subject's mean
# value, so the outcome informs the imputation), the model on each completed table,
# then Rubin's rules with the mean Satterthwaite df as the complete-data df.
# Ceiling: other covariates are not imputation predictors; mice 2lonly.pmm on the long
# table (or 2l.pan for a varying covariate) is the upgrade.
impute_feature <- \(d) {
  subjects <- summarise(
    d,
    mean_value = mean(value),
    .by = c(subject, all_of(imputed_cols))
  )
  imputation <- mice::mice(
    select(subjects, -subject),
    m = 5, # the mice default; longitudinal.py reports it
    seed = seed,
    printFlag = FALSE
  )
  fits <- map(seq_len(imputation$m), \(i) {
    completed <- bind_cols(
      select(subjects, subject),
      select(mice::complete(imputation, i), all_of(imputed_cols))
    )
    d |>
      select(-all_of(imputed_cols)) |>
      left_join(
        completed,
        by = join_by(subject),
        relationship = "many-to-one"
      ) |>
      fit_lmer()
  })
  messages <- unlist(map(fits, "messages"))
  row <- tibble::tibble(
    n_obs = nrow(d),
    n_subjects = n_distinct(d$subject),
    singular = NA,
    converged = FALSE,
    fit_message = paste(unique(messages), collapse = "; ")
  )
  fitted <- map(fits, "fit")
  if (any(map_lgl(fitted, is.null))) {
    return(bind_cols(row, na_test))
  }
  tests <- map(fitted, test_term) |>
    list_rbind()
  row <- mutate(
    row,
    singular = any(map_lgl(fitted, lme4::isSingular)),
    converged = all(map_lgl(fitted, is_converged))
  )
  if (anyNA(tests$estimate)) {
    return(bind_cols(row, na_test))
  }
  pooled <- mice::pool.scalar(
    tests$estimate,
    tests$se^2,
    n = mean(tests$df) + 1,
    k = 1
  )
  bind_cols(
    row,
    tibble::tibble(
      estimate = pooled$qbar,
      se = sqrt(pooled$t),
      df = pooled$df,
      t = pooled$qbar / sqrt(pooled$t),
      num_df = NA_real_,
      f = NA_real_,
      pvalue = 2 * pt(-abs(pooled$qbar / sqrt(pooled$t)), pooled$df)
    )
  )
}

results <- long |>
  split(long$feature) |>
  map(if (impute) impute_feature else fit_feature) |>
  list_rbind(names_to = "feature")

# Write outputs ----------------------------------------------------------------
readr::write_tsv(results, args[["out"]], na = "")
cat(R.version.string, "\n", sep = "")
cat("lme4 ", as.character(packageVersion("lme4")), "\n", sep = "")
cat("lmerTest ", as.character(packageVersion("lmerTest")), "\n", sep = "")
if (impute) {
  cat("mice ", as.character(packageVersion("mice")), "\n", sep = "")
}
