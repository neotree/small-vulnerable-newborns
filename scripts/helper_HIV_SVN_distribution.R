# ==============================================================================
# Purpose: source for the Results sentences on the SVN phenotype and
# birthweight distribution by maternal HIV status, and HIV positivity within
# each phenotype.
#
# Definitions: cohort = nnu$svn_filtered with a recorded outcome, the same base
# as helper_HIV_CFR_exposed_vs_unexposed.R, using the same HIV classification.
# Only neonates of tested mothers with a known result are included (HIV
# Positive = exposed, HIV Negative = unexposed). Distributions are column
# percentages within each HIV group; positivity = HIV Positive / tested with a
# known result within each phenotype or birthweight category. Association
# between HIV status and phenotype (or birthweight category) is tested with
# Pearson's chi-square test, by site and with both sites combined.
#
# Output: helper_hiv_svn_distribution.csv, helper_hiv_bw_distribution.csv,
#         helper_hiv_positivity_by_svn.csv, helper_hiv_positivity_by_bw.csv,
#         helper_hiv_distribution_chisq.csv
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
cat("\n=== helper_HIV_SVN_distribution.R ===\n")

# ------------------------------------------------------------------------------
# 1. Cohort and HIV classification (as in helper_HIV_CFR_exposed_vs_unexposed.R)
# ------------------------------------------------------------------------------
nnu  <- process_nnu_data(load_nnu_data())
base <- nnu$svn_filtered %>% filter(!is.na(outcome_bin))

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

tested <- base %>%
  mutate(hiv_status = classify_hiv(mathivtest, hivtestresult)) %>%
  filter(hiv_status %in% c("HIV Positive", "HIV Negative")) %>%
  mutate(hiv_group = ifelse(hiv_status == "HIV Positive", "HIV-exposed", "HIV-unexposed"),
         svn_category = as.character(svn_category),
         bw_cat = as.character(bw_cat))
cat(sprintf("Tested mothers with a known result: %d (KCH %d, SMCH %d)\n", nrow(tested),
            sum(tested$facility == "KCH"), sum(tested$facility == "SMCH")))

with_combined <- function(d) bind_rows(d, mutate(d, facility = "Combined"))

# ------------------------------------------------------------------------------
# 2. Distribution within HIV groups and positivity within categories
# ------------------------------------------------------------------------------
distribution <- function(d, var) {
  with_combined(d) %>%
    count(facility, hiv_group, category = .data[[var]], name = "n") %>%
    group_by(facility, hiv_group) %>%
    mutate(total_in_hiv_group = sum(n), pct_within_hiv_group = round(100 * n / total_in_hiv_group, 1)) %>%
    ungroup() %>%
    arrange(facility, hiv_group, category)
}
positivity <- function(d, var) {
  with_combined(d) %>%
    group_by(facility, category = .data[[var]]) %>%
    summarise(n_tested = n(), n_positive = sum(hiv_group == "HIV-exposed"),
              pct_positive = round(100 * n_positive / n_tested, 1), .groups = "drop") %>%
    arrange(facility, category)
}
svn_dist <- distribution(tested, "svn_category")
bw_dist  <- distribution(tested, "bw_cat")
svn_pos  <- positivity(tested, "svn_category")
bw_pos   <- positivity(tested, "bw_cat")

# ------------------------------------------------------------------------------
# 3. Chi-square tests of association
# ------------------------------------------------------------------------------
chisq_row <- function(d, var, site) {
  ct <- table(d$hiv_group, d[[var]])
  ct <- ct[rowSums(ct) > 0, colSums(ct) > 0, drop = FALSE]
  cs <- suppressWarnings(chisq.test(ct))
  data.frame(comparison = paste("HIV status ~", var), facility = site, n = sum(ct),
             chi_sq = round(unname(cs$statistic), 3), df = unname(cs$parameter),
             p_value = signif(cs$p.value, 3),
             note = ifelse(any(cs$expected < 5), "Some expected counts < 5", "All expected counts >= 5"))
}
chisq <- do.call(rbind, lapply(c("svn_category", "bw_cat"), function(v) rbind(
  chisq_row(tested, v, "Combined"),
  chisq_row(tested[tested$facility == "KCH", ], v, "KCH"),
  chisq_row(tested[tested$facility == "SMCH", ], v, "SMCH"))))

# ------------------------------------------------------------------------------
# 4. Save and print
# ------------------------------------------------------------------------------
write.csv(svn_dist, file.path(OUTPUT_DIR, "helper_hiv_svn_distribution.csv"), row.names = FALSE)
write.csv(bw_dist,  file.path(OUTPUT_DIR, "helper_hiv_bw_distribution.csv"), row.names = FALSE)
write.csv(svn_pos,  file.path(OUTPUT_DIR, "helper_hiv_positivity_by_svn.csv"), row.names = FALSE)
write.csv(bw_pos,   file.path(OUTPUT_DIR, "helper_hiv_positivity_by_bw.csv"), row.names = FALSE)
write.csv(chisq,    file.path(OUTPUT_DIR, "helper_hiv_distribution_chisq.csv"), row.names = FALSE)

cat("\nSVN distribution within HIV groups (both sites):\n")
print(as.data.frame(filter(svn_dist, facility == "Combined")), row.names = FALSE)
cat("\nBirthweight distribution within HIV groups (both sites):\n")
print(as.data.frame(filter(bw_dist, facility == "Combined")), row.names = FALSE)
cat("\nHIV positivity by phenotype (both sites):\n")
print(as.data.frame(filter(svn_pos, facility == "Combined")), row.names = FALSE)
cat("\nChi-square tests:\n"); print(chisq, row.names = FALSE)
cat("helper_HIV_SVN_distribution.R complete.\n")
