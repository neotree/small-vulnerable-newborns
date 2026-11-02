# ==============================================================================
# Purpose: significance tests for the distribution comparisons summarised in
# Table 3 (case fatality by SVN phenotype and by birth-weight category) and
# Table 6 (length of stay by SVN phenotype and by birth-weight category),
# both sites combined.
#
# Definitions: Kruskal-Wallis tests the rank distribution of the outcome
# across the k categories of a grouping variable; for case fatality this is
# the outcome_bin (0/1) distribution, for length of stay it is the LOS-days
# distribution among survivors (population and exclusions identical to
# Table_06_length_of_stay.R: outcome_bin == 0, stays > 180 days excluded).
# A Pearson chi-square test on the death/survive counts is also reported for
# case fatality, alongside Kruskal-Wallis, since a rank test on a binary
# outcome is a less conventional choice than a proportions test for that
# specific comparison. Population = svn_filtered, both sites combined.
# Respects POPULATION ("all" default, "USS" for the ultrasound subset).
#
# Output: helper_table03_06_kruskal_wallis_cfr_los.csv
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

suppressPackageStartupMessages({library(dplyr)})
cat("\n=== 58_kruskal_wallis_cfr_los.R ===\n")
MAX_LOS_DAYS <- 180
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
BW <- c("ELBW","VLBW","LBW","NBW","HBW")

nnu <- process_nnu_data(load_nnu_data())

# ---- CFR population (same construction as Table_03_04_distribution_and_CFR.R) ----
cfr_base <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(svn_category = factor(as.character(svn_category), levels = PH),
           bw_kg = bw_kg_canonical(.),
           bw_category = factor(ifelse(is.finite(bw_kg) & bw_kg >= BW_MIN_KG & bw_kg <= BW_MAX_KG,
                                       bw_category_from_kg(bw_kg), NA_character_), levels = BW))

# ---- LOS population (same construction as Table_06_length_of_stay.R) ----
los_base <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), outcome_bin == 0, svn_category %in% PH) %>%
    mutate(svn_category = factor(as.character(svn_category), levels = PH),
           bw_kg = bw_kg_canonical(.),
           bw_category = factor(ifelse(is.finite(bw_kg) & bw_kg >= BW_MIN_KG & bw_kg <= BW_MAX_KG,
                                       bw_category_from_kg(bw_kg), NA_character_), levels = BW),
           adm_ts = suppressWarnings(lubridate::parse_date_time(datetimeadmission, orders = c("Ymd HMS","Ymd HM","Ymd"))),
           dis_ts = suppressWarnings(lubridate::parse_date_time(datetimedischarge, orders = c("Ymd HMS","Ymd HM","Ymd"))),
           los = as.numeric(difftime(dis_ts, adm_ts, units = "days")),
           los = ifelse(is.finite(los) & los >= 0 & los <= MAX_LOS_DAYS, los, NA_real_))

run_kw <- function(df, outcome_col, group_col, outcome_label, group_label) {
    d <- df %>% filter(!is.na(.data[[outcome_col]]), !is.na(.data[[group_col]]))
    d[[group_col]] <- droplevels(d[[group_col]])
    kw <- kruskal.test(d[[outcome_col]] ~ d[[group_col]])
    tibble(
        outcome        = outcome_label,
        grouping       = group_label,
        test           = "Kruskal-Wallis",
        n              = nrow(d),
        k_groups       = nlevels(d[[group_col]]),
        statistic      = round(unname(kw$statistic), 3),
        df             = unname(kw$parameter),
        p_value        = signif(kw$p.value, 4),
        p_fmt          = ifelse(kw$p.value < 0.001, "<0.001", sprintf("%.3f", kw$p.value))
    )
}

run_chisq <- function(df, outcome_col, group_col, outcome_label, group_label) {
    d <- df %>% filter(!is.na(.data[[outcome_col]]), !is.na(.data[[group_col]]))
    d[[group_col]] <- droplevels(d[[group_col]])
    tab <- table(d[[group_col]], d[[outcome_col]])
    cs  <- suppressWarnings(chisq.test(tab, correct = FALSE))
    tibble(
        outcome        = outcome_label,
        grouping       = group_label,
        test           = "Pearson chi-square",
        n              = nrow(d),
        k_groups       = nrow(tab),
        statistic      = round(unname(cs$statistic), 3),
        df             = unname(cs$parameter),
        p_value        = signif(cs$p.value, 4),
        p_fmt          = ifelse(cs$p.value < 0.001, "<0.001", sprintf("%.3f", cs$p.value))
    )
}

results <- bind_rows(
    run_kw(cfr_base, "outcome_bin", "svn_category", "Case fatality", "SVN phenotype"),
    run_kw(filter(cfr_base, !is.na(bw_category)), "outcome_bin", "bw_category", "Case fatality", "Birth-weight category"),
    run_chisq(cfr_base, "outcome_bin", "svn_category", "Case fatality", "SVN phenotype"),
    run_chisq(filter(cfr_base, !is.na(bw_category)), "outcome_bin", "bw_category", "Case fatality", "Birth-weight category"),
    run_kw(los_base, "los", "svn_category", "Length of stay (survivors)", "SVN phenotype"),
    run_kw(filter(los_base, !is.na(bw_category)), "los", "bw_category", "Length of stay (survivors)", "Birth-weight category")
)

write.csv(results, file.path(OUT2, "helper_table03_06_kruskal_wallis_cfr_los.csv"), row.names = FALSE)

cat("\nCase fatality and length-of-stay distribution tests, both sites combined:\n")
print(as.data.frame(results), row.names = FALSE)
cat("\nSaved helper_table03_06_kruskal_wallis_cfr_los.csv\n")
cat("58_kruskal_wallis_cfr_los.R complete.\n")
