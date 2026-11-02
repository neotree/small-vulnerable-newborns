# ==============================================================================
# Purpose: Table 10 -- case fatality by resuscitation requirement at birth,
# on the full analysis cohort.
#
# Definitions: resuscitation categories come from resus_category() in
# 00a_setup_shared.R. "Unknown" (the form's Unknown option, or blank) is shown
# as its own row so every admission is accounted for, rather than folding it
# into "None". Reports n, deaths and CFR per 1,000 admissions with 95% Wilson CI.
#
# Output: table_10_cfr_by_resuscitation.csv
# ==============================================================================
# --- Locate and source the shared setup (self-contained; runs from anywhere) --
SCRIPT_DIR <- local({
    a <- commandArgs(FALSE); f <- a[grepl("^--file=", a)]
    src <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
    if (length(src)) dirname(normalizePath(tail(src, 1)[[1]]))
    else if (length(f)) dirname(normalizePath(sub("^--file=", "", f[1])))
    else if (file.exists("00a_setup_shared.R")) getwd()
    else normalizePath("scripts")
})
setwd(SCRIPT_DIR)
source("00a_setup_shared.R")
OUT2 <- OUTPUT_DIR
suppressPackageStartupMessages({library(dplyr); library(stringr)})
cat("\n=== Table_10_cfr_by_resuscitation.R ===\n")

# Wilson 95% CI for a proportion, returned on the per-1,000 scale.
wilson1000 <- function(x, n){
    if (is.na(n) || n == 0) return(c(NA_real_, NA_real_))
    p <- x / n; z <- 1.96; d <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / d
    hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
    c((ctr - hw) * 1000, (ctr + hw) * 1000)
}

PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")

# ------------------------------------------------------------------------------
# Analysis cohort with resuscitation categorised (Unknown kept)
# ------------------------------------------------------------------------------
nnu <- process_nnu_data(load_nnu_data())
if (!"resus" %in% names(nnu$svn_filtered))
    stop("Column 'resus' not found in svn_filtered.")

cohort <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(
        resus_cat = recode(as.character(resus_category(resus)),
            "Stimulation/Suction/O2" = "Stimulation / suction / oxygen",
            "BVM"                    = "Bag-valve-mask ventilation",
            "CPR"                    = "Cardiopulmonary resuscitation"),
        resus_cat = factor(resus_cat,
            levels = c("None", "Stimulation / suction / oxygen",
                       "Bag-valve-mask ventilation", "Cardiopulmonary resuscitation",
                       "Unknown")))

# ------------------------------------------------------------------------------
# n, deaths, CFR per 1,000 (95% Wilson CI) by category + Total
# ------------------------------------------------------------------------------
by_cat <- cohort %>%
    group_by(resus_cat, .drop = FALSE) %>%
    summarise(n = n(), deaths = sum(outcome_bin == 1), .groups = "drop")

tot <- data.frame(resus_cat = "Total", n = nrow(cohort),
                  deaths = sum(cohort$outcome_bin == 1), stringsAsFactors = FALSE)

tab <- bind_rows(as.data.frame(by_cat), tot) %>%
    rowwise() %>%
    mutate(cfr_per1000 = ifelse(n > 0, round(1000 * deaths / n), NA_real_),
           ci_lo = round(wilson1000(deaths, n)[1]),
           ci_hi = round(wilson1000(deaths, n)[2]),
           cfr_ci = ifelse(is.na(cfr_per1000), "-",
                           sprintf("%d (%d-%d)", cfr_per1000, ci_lo, ci_hi))) %>%
    ungroup() %>%
    transmute(`Resuscitation at birth` = as.character(resus_cat),
              `Admissions, n` = n,
              `Deaths, n` = deaths,
              `CFR per 1,000 (95% CI)` = cfr_ci)

write.csv(tab, file.path(OUT2, "table_10_cfr_by_resuscitation.csv"), row.names = FALSE)
cat("Saved: table_10_cfr_by_resuscitation.csv\n")
print(tab, row.names = FALSE)

# Reconciliation check
stopifnot(sum(by_cat$n) == nrow(cohort),
          sum(by_cat$deaths) == sum(cohort$outcome_bin == 1))
cat(sprintf("Reconciles: admissions %d, deaths %d\n",
            nrow(cohort), sum(cohort$outcome_bin == 1)))
cat("Table_10_cfr_by_resuscitation.R complete.\n")
