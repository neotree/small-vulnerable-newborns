# ==============================================================================
# Purpose: source for numbers quoted in the manuscript text (Abstract, Methods,
# Results, Discussion and supplementary table notes) that are derived from,
# rather than printed in, the other output tables: annual averages, exclusion
# counts, site shares, the normal-birthweight SVN count, timing-of-death
# denominators, ANC exact-count missingness, the S7 steroid-status split, the
# HIV distribution sentence and the adjusted HIV odds ratio.
#
# Definitions: every value is read from the output tables written earlier in
# config/scripts.csv, except the maternal delivery totals before exclusions,
# which are counted from the locked maternal inputs within the study period.
# Annual averages divide by the number of calendar years in the study period.
# Percentages are rounded to one decimal place unless stated in the item.
#
# Output: manuscript_text_values.csv (section, item, value, reported, source)
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
cat("\n=== Manuscript_text_values.R ===\n")

tab <- function(name) read.csv(file.path(OUTPUT_DIR, name), stringsAsFactors = FALSE, check.names = FALSE)
lead_n <- function(x) as.numeric(gsub(",", "", sub("^\\s*([0-9,]+).*$", "\\1", x)))
fmt_n <- function(x) formatC(x, format = "d", big.mark = ",")
fmt_p <- function(x, d = 1) formatC(x, format = "f", digits = d)
rows <- list()
add <- function(section, item, value, reported, source)
    rows[[length(rows) + 1]] <<- data.frame(section = section, item = item, value = value,
                                            reported = reported, source = source)
n_years <- as.numeric(format(STUDY_END, "%Y")) - as.numeric(format(STUDY_START, "%Y")) + 1

# ------------------------------------------------------------------------------
# 1. Methods: fetal estimated-weight proxies at 22+0 and 23+0 weeks
# ------------------------------------------------------------------------------
add("Methods", "Fetal EFW 10th centile proxy at 22+0 weeks (g)", PCT10_UNISEX[["22"]],
    fmt_n(PCT10_UNISEX[["22"]]), "00a_setup_shared.R")
add("Methods", "Fetal EFW 10th centile proxy at 23+0 weeks (g)", PCT10_UNISEX[["23"]],
    fmt_n(PCT10_UNISEX[["23"]]), "00a_setup_shared.R")

# ------------------------------------------------------------------------------
# 2. Maternal deliveries: totals, exclusions and annual averages
# ------------------------------------------------------------------------------
mat_raw <- load_maternal_data()
t1a <- tab("table_01a_maternal_characteristics.csv")
incl <- t1a[t1a$Variable == "Deliveries (valid BW + GA), n", ]
for (site in c("KCH", "SMCH", "Combined")) {
    total <- if (site == "Combined") nrow(mat_raw) else sum(mat_raw$facility == site)
    kept <- as.numeric(incl[[site]]); excl <- total - kept
    add("Results", paste0("Deliveries before exclusions, ", site), total, fmt_n(total),
        "maternal inputs (study period)")
    add("Results", paste0("Deliveries excluded for missing/invalid birthweight or gestation, ", site),
        excl, sprintf("%s (%s%%)", fmt_n(excl), fmt_p(100 * excl / total)),
        "maternal inputs; table_01a_maternal_characteristics.csv")
    if (site != "Combined")
        add("Results", paste0("Average deliveries per year, ", site), total / n_years,
            fmt_n(round(total / n_years)), "maternal inputs (study period)")
}

allb <- tab("helper_table01_svn_by_site_allbirths.csv")
for (site in c("KCH", "SMCH", "Combined")) {
    n <- lead_n(allb[[site]]); names(n) <- allb$svn_category
    svn <- sum(n[c("Term-SGA", "Preterm-AGA", "Preterm-SGA")])
    add("Abstract/Results", paste0("Births meeting the SVN definition, ", site), 100 * svn / sum(n),
        sprintf("%s%% (%s/%s)", fmt_p(100 * svn / sum(n)), fmt_n(svn), fmt_n(sum(n))),
        "helper_table01_svn_by_site_allbirths.csv")
}

# ------------------------------------------------------------------------------
# 3. Neonatal-unit admissions: annual averages and SVN share by site
# ------------------------------------------------------------------------------
t2 <- tab("table_02_nnu_characteristics.csv")
adm <- t2[t2$Variable == "N admissions", ]
svn_rows <- t2[t2$Variable == "SVN phenotype", ]
for (site in c("KCH", "SMCH", "Combined")) {
    n_adm <- as.numeric(adm[[site]])
    if (site != "Combined")
        add("Results", paste0("Average admissions per year, ", site), n_adm / n_years,
            fmt_n(round(n_adm / n_years)), "table_02_nnu_characteristics.csv")
    svn <- sum(lead_n(svn_rows[[site]][svn_rows$Category != "Term-AGA"]))
    add("Abstract/Results", paste0("SVN admissions, ", site), 100 * svn / n_adm,
        sprintf("%s%% (%s/%s)", fmt_p(100 * svn / n_adm), fmt_n(svn), fmt_n(n_adm)),
        "table_02_nnu_characteristics.csv")
}

# ------------------------------------------------------------------------------
# 4. Normal-birthweight admissions classified as SVN (Supplementary Table S1)
# ------------------------------------------------------------------------------
s1 <- tab("manuscript_supplementary_s1_counts_percentages.csv")
nbw <- s1[s1$bw_category == "NBW", ]; nb <- setNames(nbw$n, nbw$phenotype)
nbw_svn <- sum(nb[c("Term-SGA", "Preterm-AGA", "Preterm-SGA")])
add("Abstract/Results/Discussion", "Normal-birthweight admissions classified as SVN", nbw_svn,
    sprintf("%s%% (%s/%s)", fmt_p(100 * nbw_svn / nb[["Row total"]]), fmt_n(nbw_svn), fmt_n(nb[["Row total"]])),
    "manuscript_supplementary_s1_counts_percentages.csv")
for (p in c("Term-SGA", "Preterm-AGA"))
    add("Results", paste0("Normal-birthweight SVN admissions that were ", p), nb[[p]],
        sprintf("%s%% (%s/%s)", fmt_p(100 * nb[[p]] / nbw_svn), fmt_n(nb[[p]]), fmt_n(nbw_svn)),
        "manuscript_supplementary_s1_counts_percentages.csv")

# ------------------------------------------------------------------------------
# 5. Timing of death: recorded-age denominators (Results and S2a note)
# ------------------------------------------------------------------------------
t5 <- tab("table_05_timing_by_site.csv")
deaths <- setNames(lead_n(t2[t2$Variable == "Deaths, n (per 1,000)", c("KCH", "SMCH", "Combined")]),
                   c("KCH", "SMCH", "Both sites"))
admissions <- setNames(as.numeric(adm[c("KCH", "SMCH", "Combined")]), c("KCH", "SMCH", "Both sites"))
for (site in c("KCH", "SMCH", "Both sites")) {
    r <- t5[t5$group == site, ]
    unknown <- lead_n(r[["Unknown age"]]); known <- deaths[[site]] - unknown
    add("Results/S2a", paste0("Deaths with recorded age, ", site), known,
        sprintf("%s of %s (%s%%)", fmt_n(known), fmt_n(deaths[[site]]), fmt_p(100 * known / deaths[[site]])),
        "table_05_timing_by_site.csv; table_02_nnu_characteristics.csv")
    add("Results/S2a", paste0("Deaths with unknown age, ", site), unknown,
        sprintf("%s (%s%%)", fmt_n(unknown), fmt_p(100 * unknown / deaths[[site]])),
        "table_05_timing_by_site.csv")
    add("S2a", paste0("Recorded-age CFR per 1,000 admissions, ", site), 1000 * known / admissions[[site]],
        fmt_p(1000 * known / admissions[[site]]), "table_05_timing_by_site.csv; table_02_nnu_characteristics.csv")
    if (site == "Both sites") {
        d01 <- lead_n(r[["0-1 days"]]); d13 <- lead_n(r[["1-3 days"]]); d37 <- lead_n(r[["3-7 days"]])
        add("Results", "Deaths on day 0-1, % of recorded-age deaths", d01,
            sprintf("%s%% (n=%s)", fmt_p(100 * d01 / known), fmt_n(d01)), "table_05_timing_by_site.csv")
        add("Results", "Deaths on day 1-3, % of recorded-age deaths", d13,
            sprintf("%s%% (n=%s)", fmt_p(100 * d13 / known), fmt_n(d13)), "table_05_timing_by_site.csv")
        add("Results/Discussion", "Deaths within three days, % of recorded-age deaths", 100 * (d01 + d13) / known,
            sprintf("%s%%", fmt_p(100 * (d01 + d13) / known)), "table_05_timing_by_site.csv")
        add("Results/Discussion", "Deaths within the first week, % of recorded-age deaths",
            100 * (d01 + d13 + d37) / known, sprintf("%s%%", fmt_p(100 * (d01 + d13 + d37) / known)),
            "table_05_timing_by_site.csv")
    }
}

# ------------------------------------------------------------------------------
# 6. ANC: exact visit count missing
# ------------------------------------------------------------------------------
anc <- tab("manuscript_anc_missingness_check.csv")
add("Results", "ANC exact visit count missing", sum(anc$exact_missing),
    sprintf("%s (%s%%)", fmt_n(sum(anc$exact_missing)), fmt_p(100 * sum(anc$exact_missing) / sum(anc$n))),
    "manuscript_anc_missingness_check.csv")
add("Results", "ANC category recovered without an exact count", sum(anc$category_only),
    fmt_n(sum(anc$category_only)), "manuscript_anc_missingness_check.csv")

# ------------------------------------------------------------------------------
# 7. Antenatal steroids: status split in the eligible window (S7 note)
# ------------------------------------------------------------------------------
comp <- tab("table_11_ga24to34wk_06_steroid_completeness.csv")
both <- comp[comp$facility == "Both", ]; window_n <- sum(both$n)
known <- sum(both$n[both$steroid_record %in% c("Yes (Y)", "No (N)")])
s7 <- tab("supplementary_table_s7_steroids_with_pooled.csv")
s7b <- s7[grepl("^Both sites", s7$Group), ]
den <- function(x) as.numeric(gsub(",", "", sub(".*/([0-9,]+)\\).*", "\\1", x)))
svn_known <- sum(den(s7b[[2]]) + den(s7b[[3]]))
add("S7", "Eligible-window admissions (24+0 to 34+6 weeks)", window_n, fmt_n(window_n),
    "table_11_ga24to34wk_06_steroid_completeness.csv")
add("S7", "SVN-classifiable neonates with known steroid status", svn_known,
    sprintf("%s (%s%%)", fmt_n(svn_known), fmt_p(100 * svn_known / window_n)),
    "supplementary_table_s7_steroids_with_pooled.csv")
add("S7", "Known steroid status but not SVN-classifiable", known - svn_known, fmt_n(known - svn_known),
    "table_11_ga24to34wk_06_steroid_completeness.csv; supplementary_table_s7_steroids_with_pooled.csv")
add("S7", "Steroid status unknown or missing", window_n - known,
    sprintf("%s (%s%%)", fmt_n(window_n - known), fmt_p(100 * (window_n - known) / window_n)),
    "table_11_ga24to34wk_06_steroid_completeness.csv")

# ------------------------------------------------------------------------------
# 8. Maternal HIV: distribution sentence and adjusted odds ratio
# ------------------------------------------------------------------------------
hc <- tab("helper_hiv_distribution_chisq.csv")
hc <- hc[hc$facility == "Combined" & hc$comparison == "HIV status ~ svn_category", ]
add("Results", "HIV status by SVN phenotype, chi-square (both sites)", hc$chi_sq,
    sprintf("chi-square=%s, p=%s", fmt_p(hc$chi_sq), fmt_p(hc$p_value, 3)), "helper_hiv_distribution_chisq.csv")
hd <- tab("helper_hiv_svn_distribution.csv"); hd <- hd[hd$facility == "Combined", ]
hb <- tab("helper_hiv_bw_distribution.csv");  hb <- hb[hb$facility == "Combined", ]
pair <- function(d, cat) sprintf("%s%% vs %s%%",
    fmt_p(d$pct_within_hiv_group[d$category == cat & d$hiv_group == "HIV-exposed"]),
    fmt_p(d$pct_within_hiv_group[d$category == cat & d$hiv_group == "HIV-unexposed"]))
for (p in c("Term-AGA", "Term-SGA"))
    add("Results", paste0(p, ", % of HIV-exposed vs unexposed"), NA, pair(hd, p), "helper_hiv_svn_distribution.csv")
add("Results", "LBW (1,500-2,499 g), % of HIV-exposed vs unexposed", NA, pair(hb, "LBW"),
    "helper_hiv_bw_distribution.csv")
hp <- tab("helper_hiv_positivity_by_svn.csv"); hp <- hp[hp$facility == "Combined", ]
add("Results", "HIV positivity range across phenotypes", NA,
    sprintf("%s %s%% to %s %s%%", hp$category[which.min(hp$pct_positive)], fmt_p(min(hp$pct_positive)),
            hp$category[which.max(hp$pct_positive)], fmt_p(max(hp$pct_positive))),
    "helper_hiv_positivity_by_svn.csv")
mv <- tab("41b_wholecohort_multivariable_tempcontinuous.csv")
hiv <- mv[mv$term == "hiv_statusHIV Positive", ]
add("Results", "Adjusted OR, HIV-positive vs negative (whole cohort)", hiv$OR,
    sprintf("aOR %.2f, 95%% CI %.2f-%.2f", hiv$OR, hiv$ci_lower, hiv$ci_upper),
    "41b_wholecohort_multivariable_tempcontinuous.csv")

# ------------------------------------------------------------------------------
# 9. Save and print
# ------------------------------------------------------------------------------
values <- do.call(rbind, rows)
write.csv(values, file.path(OUTPUT_DIR, "manuscript_text_values.csv"), row.names = FALSE)
print(values[c("section", "item", "reported")], row.names = FALSE, right = FALSE)
cat("Manuscript_text_values.R complete.\n")
