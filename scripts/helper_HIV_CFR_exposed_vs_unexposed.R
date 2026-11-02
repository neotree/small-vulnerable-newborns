# ==============================================================================
# Purpose: source for the Results sentence on maternal HIV testing rate,
# HIV positivity by site and case fatality in exposed versus unexposed
# neonates.
#
# Definitions: cohort = nnu$svn_filtered with a recorded outcome, the same base
# as helper_HIV_SVN_distribution.R. HIV status is classified with the same rule
# as that script:
#   HIV Positive = tested and result "R"; HIV Negative = tested and "NR";
#   Tested (unknown result); Not Tested; Missing.
# Exposed = HIV Positive; unexposed = HIV Negative. Neonates whose mothers
# were not tested or whose result is unknown are excluded from the CFR
# comparison. Case fatality = deaths per 1,000 admissions with Wilson 95% CI.
# p-values: Pearson chi-square with continuity correction (R default for a
# 2x2 table) and Fisher's exact test; crude OR (exposed vs unexposed) with
# Wald 95% CI.
#
# Output: helper_hiv_cfr_exposed_vs_unexposed.csv,
#         helper_hiv_testing_and_positivity.csv
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
suppressPackageStartupMessages(library(dplyr))
cat("\n=== helper_HIV_CFR_exposed_vs_unexposed.R ===\n")

# ------------------------------------------------------------------------------
# 1. Cohort (same base as helper_HIV_SVN_distribution.R)
# ------------------------------------------------------------------------------
nnu  <- process_nnu_data(load_nnu_data())
base <- nnu$svn_filtered %>% filter(!is.na(outcome_bin))
cat(sprintf("Cohort: %d (KCH %d, SMCH %d)\n", nrow(base),
            sum(base$facility == "KCH"), sum(base$facility == "SMCH")))

# ------------------------------------------------------------------------------
# 2. HIV classification (same rule as helper_HIV_SVN_distribution.R)
# ------------------------------------------------------------------------------
classify_hiv <- function(mathivtest_col, hivtestresult_col) {
  tested  <- toupper(trimws(as.character(mathivtest_col))) %in% c("TRUE","T","YES","Y","1")
  not_tested <- toupper(trimws(as.character(mathivtest_col))) %in% c("FALSE","F","NO","N","0")
  result  <- toupper(trimws(as.character(hivtestresult_col)))

  dplyr::case_when(
    tested  & result == "R"                          ~ "HIV Positive",
    tested  & result == "NR"                         ~ "HIV Negative",
    tested  & !result %in% c("R","NR")               ~ "Tested (unknown result)",
    not_tested                                        ~ "Not Tested",
    TRUE                                              ~ "Missing"
  )
}

df <- base %>% mutate(hiv_status = classify_hiv(mathivtest, hivtestresult))
cat("\nHIV status by site:\n")
print(table(df$facility, df$hiv_status, useNA = "ifany"))

# ------------------------------------------------------------------------------
# 3. Testing rate and positivity (first half of the sentence)
#    testing rate = tested (any result) / all neonates in the cohort
#    positivity   = HIV Positive / (HIV Positive + HIV Negative)
# ------------------------------------------------------------------------------
testing <- function(d, site) {
    tested <- d$hiv_status %in% c("HIV Positive", "HIV Negative", "Tested (unknown result)")
    pos <- sum(d$hiv_status == "HIV Positive"); neg <- sum(d$hiv_status == "HIV Negative")
    data.frame(site = site, n_cohort = nrow(d), n_tested = sum(tested),
               testing_rate_pct = round(100 * mean(tested), 1),
               n_positive = pos, n_negative = neg,
               positivity_pct = round(100 * pos / (pos + neg), 1))
}
tab_test <- rbind(testing(df[df$facility == "KCH", ], "KCH"),
                  testing(df[df$facility == "SMCH", ], "SMCH"),
                  testing(df, "Both sites"))

# ------------------------------------------------------------------------------
# 4. Case fatality, exposed vs unexposed (second half of the sentence)
# ------------------------------------------------------------------------------
wilson <- function(x, n, z = 1.959964) {
    p <- x / n; den <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / den
    hw  <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
    c(ctr - hw, ctr + hw)
}
cfr_compare <- function(d, site) {
    e <- d[d$hiv_status == "HIV Positive", ]; u <- d[d$hiv_status == "HIV Negative", ]
    de <- sum(e$outcome_bin == 1); ne <- nrow(e)
    du <- sum(u$outcome_bin == 1); nu <- nrow(u)
    tab <- matrix(c(de, ne - de, du, nu - du), nrow = 2, byrow = TRUE)
    ci_e <- 1000 * wilson(de, ne); ci_u <- 1000 * wilson(du, nu)
    or <- (de * (nu - du)) / ((ne - de) * du)
    se <- sqrt(1 / de + 1 / (ne - de) + 1 / du + 1 / (nu - du))
    data.frame(
        site = site,
        exposed_deaths = de, exposed_n = ne,
        exposed_cfr_per1000 = round(1000 * de / ne, 1),
        exposed_ci = sprintf("%.0f-%.0f", ci_e[1], ci_e[2]),
        unexposed_deaths = du, unexposed_n = nu,
        unexposed_cfr_per1000 = round(1000 * du / nu, 1),
        unexposed_ci = sprintf("%.0f-%.0f", ci_u[1], ci_u[2]),
        OR = round(or, 2),
        OR_ci = sprintf("%.2f-%.2f", exp(log(or) - 1.959964 * se), exp(log(or) + 1.959964 * se)),
        p_chisq_yates = signif(suppressWarnings(chisq.test(tab))$p.value, 3),
        p_chisq_no_correction = signif(suppressWarnings(chisq.test(tab, correct = FALSE))$p.value, 3),
        p_fisher = signif(fisher.test(tab)$p.value, 3))
}
tab_cfr <- rbind(cfr_compare(df[df$facility == "KCH", ], "KCH"),
                 cfr_compare(df[df$facility == "SMCH", ], "SMCH"),
                 cfr_compare(df, "Both sites"))

# ------------------------------------------------------------------------------
# 5. Save and print
# ------------------------------------------------------------------------------
write.csv(tab_test, file.path(OUTPUT_DIR, "helper_hiv_testing_and_positivity.csv"), row.names = FALSE)
write.csv(tab_cfr,  file.path(OUTPUT_DIR, "helper_hiv_cfr_exposed_vs_unexposed.csv"), row.names = FALSE)

cat("\nTesting rate and positivity:\n"); print(tab_test, row.names = FALSE)
cat("\nCase fatality, HIV-exposed vs unexposed:\n"); print(tab_cfr, row.names = FALSE)
b <- tab_cfr[tab_cfr$site == "Both sites", ]
cat(sprintf("\nManuscript check: %.0f per 1000 (exposed) vs %.0f (unexposed), p=%.2f (chi-square)\n",
            b$exposed_cfr_per1000, b$unexposed_cfr_per1000, b$p_chisq_yates))
cat(sprintf("Saved: %s\n", file.path(OUTPUT_DIR, "helper_hiv_cfr_exposed_vs_unexposed.csv")))
cat("helper_HIV_CFR_exposed_vs_unexposed.R complete.\n")
