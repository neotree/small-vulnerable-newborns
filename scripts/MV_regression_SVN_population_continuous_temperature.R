# ==============================================================================
# Purpose: variant of MV_regression_SVN_population.R -- same adjustment set
# and SVN population, but admission temperature entered as a continuous
# per-1degC term (temp_val) instead of the categorical temp_cat, giving the
# adjusted per-degree odds ratio to sit alongside the univariable per-degree
# estimate already reported in Supplementary_Table_S1_validation_CFR_full_vs_USS.R.
#
# Definitions: population restricted to the three SVN phenotypes (Term-SGA,
# Preterm-AGA, Preterm-SGA); Term-AGA is excluded. Reference phenotype =
# Term-SGA. Covariates and coding are otherwise identical to
# MV_regression_SVN_population.R (same exclusions: birth weight, gestation
# band, pre-eclampsia). See MV_regression_whole_cohort_continuous_temperature.R
# for the same model with Term-AGA as the phenotype reference.
#
# Output: 41_svn_population_counts.csv, 41_svn_univariable_tempcontinuous.csv,
#         41_svn_multivariable_tempcontinuous.csv (+ .xlsx),
#         41_svn_multivariable_tempcontinuous_missingness.csv
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

cat("\n=== MV_regression_SVN_population_continuous_temperature.R ===\n")

OUT2 <- OUTPUT_DIR
TEMP_VALID_MIN <- 30.0
TEMP_VALID_MAX <- 45.0

# ==============================================================================
# LOAD DATA AND RESTRICT TO THE SVN POPULATION
# ==============================================================================
cat("Loading NNU data...\n")
raw  <- load_nnu_data()
nnu  <- process_nnu_data(raw)

SVN3 <- c("Term-SGA", "Preterm-AGA", "Preterm-SGA")   # Term-AGA EXCLUDED

data <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin)) %>%
    filter(svn_category %in% SVN3) %>%
    mutate(
        svn_cat3 = factor(as.character(svn_category), levels = SVN3),
        facility = factor(as.character(facility), levels = c("SMCH", "KCH"))
    )

# Verify the denominator (expect total 17,931; KCH 6,758; SMCH 11,173)
counts <- data %>%
    group_by(facility, svn_cat3) %>%
    summarise(n = n(), deaths = sum(outcome_bin == 1), .groups = "drop")
totals <- data %>%
    summarise(facility = "TOTAL", svn_cat3 = "ALL SVN",
              n = n(), deaths = sum(outcome_bin == 1))
site_tot <- data %>%
    group_by(facility) %>%
    summarise(svn_cat3 = "ALL SVN", n = n(),
              deaths = sum(outcome_bin == 1), .groups = "drop")
count_tbl <- bind_rows(counts, site_tot, totals)
write.csv(count_tbl, file.path(OUT2, "41_svn_population_counts.csv"),
          row.names = FALSE)

cat(sprintf("SVN population: %d  (KCH %d, SMCH %d); deaths %d\n",
            nrow(data),
            sum(data$facility == "KCH"), sum(data$facility == "SMCH"),
            sum(data$outcome_bin == 1)))
cat("  (Expect total 17,931; KCH 6,758; SMCH 11,173 - verify against 41_svn_population_counts.csv)\n")

# ==============================================================================
# DERIVE RISK-FACTOR VARIABLES (same coding as MV_regression_SVN_population.R;
# temp_val kept continuous here rather than collapsed into temp_cat). Birth-
# weight and gestation-band predictors are intentionally NOT derived here.
# ==============================================================================
data <- data %>%
    mutate(
        # HIV
        hiv_status = case_when(
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y")  &
                toupper(trimws(as.character(hivtestresult))) == "R"  ~ "HIV Positive",
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y")  &
                toupper(trimws(as.character(hivtestresult))) == "NR" ~ "HIV Negative",
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y")  ~ "HIV Tested (unknown result)",
            toupper(trimws(as.character(mathivtest))) %in% c("FALSE","N") ~ "Not Tested",
            TRUE ~ NA_character_
        ),
        hiv_status = factor(hiv_status,
                            levels = c("HIV Negative","HIV Positive",
                                       "HIV Tested (unknown result)","Not Tested")),

        # Parity
        par_num = suppressWarnings(as.numeric(par)),
        parity = factor(case_when(
            is.na(par_num) ~ NA_character_,
            par_num == 1   ~ "Primiparous (1st delivery)",
            par_num <= 4   ~ "Multiparous (2nd-4th)",
            TRUE           ~ "Grand Multipara (5th+)"
        ), levels = c("Primiparous (1st delivery)",
                      "Multiparous (2nd-4th)","Grand Multipara (5th+)")),

        # Maternal age (reference 20-24 years)
        matage_num = suppressWarnings(as.numeric(matageyrs)),
        matage_cat = case_when(
            is.na(matage_num)  ~ NA_character_,
            matage_num < 16    ~ "<16 years",
            matage_num < 20    ~ "16-19 years",
            matage_num < 25    ~ "20-24 years",
            matage_num < 35    ~ "25-34 years",
            TRUE               ~ ">=35 years"
        ),
        matage_cat = factor(matage_cat,
                            levels = c("20-24 years","<16 years","16-19 years",
                                       "25-34 years",">=35 years")),

        # ANC (WHO threshold 4 visits)
        anc_num = suppressWarnings(as.numeric(antenatalcare)),
        anc = anc_factor(antenatalcare, antenatalcare_category),

        # Type of pregnancy
        typebirth_cat = case_when(
            is.na(typebirth) | trimws(typebirth) == "" ~ NA_character_,
            toupper(trimws(typebirth)) %in% c("S","SINGLETON") ~ "Singleton",
            toupper(trimws(typebirth)) %in% c("TW1","TW2","TWIN1","TWIN2","TWIN") ~ "Twin",
            toupper(trimws(typebirth)) %in% c("TR1","TR2","TR3","TRIPLET") ~ "Triplet",
            TRUE ~ NA_character_
        ),
        typebirth_cat = factor(typebirth_cat, levels = c("Singleton","Twin","Triplet")),

        # Mode of delivery
        modedelivery_num = suppressWarnings(as.numeric(modedelivery)),
        mode_del = factor(case_when(
            is.na(modedelivery_num) ~ NA_character_,
            modedelivery_num == 1   ~ "SVD",
            modedelivery_num %in% c(2,3,6,7) ~ "Assisted Vaginal",
            modedelivery_num == 4   ~ "ELCS",
            modedelivery_num %in% c(5,8) ~ "EMCS",
            TRUE ~ NA_character_
        ), levels = c("SVD","Assisted Vaginal","ELCS","EMCS")),

        # Resuscitation
        resus_cat = resus_category(resus),

        # Inborn / Outborn
        inout_cat = case_when(
            !is.na(inorout) ~ case_when(
                toupper(trimws(as.character(inorout))) %in% c("TRUE","INBORN","1","IN","YES") ~ "Inborn",
                toupper(trimws(as.character(inorout))) %in% c("FALSE","OUTBORN","0","OUT","NO") ~ "Outborn",
                TRUE ~ NA_character_),
            !is.na(admittedfrom) ~ case_when(
                toupper(trimws(as.character(admittedfrom))) %in% c("LW","LABOUR WARD","LW ") ~ "Inborn",
                TRUE ~ "Outborn"),
            TRUE ~ NA_character_
        ),
        inout_cat = factor(inout_cat, levels = c("Inborn","Outborn")),

        # Apgar 5 min
        apgar5_num = suppressWarnings(as.numeric(apgar5)),
        apgar5_num = ifelse(apgar5_num >= 0 & apgar5_num <= 10, apgar5_num, NA_real_),
        apgar5_cat = factor(case_when(
            is.na(apgar5_num)  ~ NA_character_,
            apgar5_num <= 3    ~ "Severely depressed (0-3)",
            apgar5_num <= 6    ~ "Moderately depressed (4-6)",
            TRUE               ~ "Normal (7-10)"
        ), levels = c("Normal (7-10)","Moderately depressed (4-6)","Severely depressed (0-3)")),

        # Sex
        sex_cat = case_when(
            toupper(trimws(as.character(gender_std))) %in% c("M","MALE","BOY")    ~ "Male",
            toupper(trimws(as.character(gender_std))) %in% c("F","FEMALE","GIRL") ~ "Female",
            TRUE ~ NA_character_
        ),
        sex_cat = factor(sex_cat, levels = c("Male","Female")),

        # Temperature -- CONTINUOUS (per 1 degC), not categorical
        temp_val = suppressWarnings(as.numeric(temperature)),
        temp_val = ifelse(temp_val >= TEMP_VALID_MIN & temp_val <= TEMP_VALID_MAX,
                          temp_val, NA_real_)
    )

# ==============================================================================
# COVARIATE SET -- identical to MV_regression_SVN_population.R except temp_cat
# is replaced by the continuous temp_val.
# ==============================================================================
selected_vars <- c("svn_cat3", "facility", "hiv_status", "parity", "matage_cat",
                    "anc", "typebirth_cat", "mode_del",
                    "resus_cat", "inout_cat", "apgar5_cat", "sex_cat", "temp_val")

NOTREC_VARS <- c("parity", "matage_cat", "anc", "apgar5_cat")

cat("\nCovariates in the multivariable model (adjustment set):\n")
cat(paste0("  - ", selected_vars, collapse = "\n"), "\n")
writeLines(selected_vars, file.path(OUT2, "41_covariate_list.txt"))

# ==============================================================================
# UNIVARIABLE ESTIMATES -- generalised to handle temp_val (continuous), which
# has no reference level to bind; all other predictors are handled exactly as
# in MV_regression_SVN_population.R.
# ==============================================================================
run_univar <- function(df, predictor_var, outcome_var = "outcome_bin") {
    df2 <- df %>%
        select(all_of(c(outcome_var, predictor_var))) %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[predictor_var]]))
    is_num <- is.numeric(df2[[predictor_var]])
    insufficient <- nrow(df2) < 10 || length(unique(df2[[outcome_var]])) < 2 ||
        (is_num  && length(unique(df2[[predictor_var]])) < 2) ||
        (!is_num && (!is.factor(df2[[predictor_var]]) ||
                     length(levels(droplevels(df2[[predictor_var]]))) < 2))
    if (insufficient) {
        return(tibble(predictor = predictor_var, term = NA_character_,
                      n_total = nrow(df2), n_events = sum(df2[[outcome_var]] == 1),
                      OR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                      p_value = NA_real_, note = "Insufficient data"))
    }
    fmla <- as.formula(paste(outcome_var, "~", predictor_var))
    fit  <- glm(fmla, data = df2, family = binomial(link = "logit"))
    res  <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)")
    if (is_num) {
        res <- res %>% mutate(note = "OR per 1 degC increase")
    } else {
        ref_level <- levels(df2[[predictor_var]])[1]
        ref_row <- tibble(term = paste0(predictor_var, ref_level),
                          estimate = 1.0, conf.low = NA_real_,
                          conf.high = NA_real_, p.value = NA_real_)
        res <- res %>% mutate(note = "") %>% bind_rows(ref_row) %>%
            mutate(note = ifelse(is.na(note) & term == paste0(predictor_var, ref_level),
                                 "Reference", note))
    }
    res %>%
        mutate(predictor = predictor_var, n_total = nrow(df2),
               n_events = sum(df2[[outcome_var]] == 1),
               OR = round(estimate, 3), ci_lower = round(conf.low, 3),
               ci_upper = round(conf.high, 3), p_value = round(p.value, 4)) %>%
        select(predictor, term, n_total, n_events, OR, ci_lower, ci_upper, p_value, note)
}

cat("\nRunning univariable regressions...\n")
uni <- bind_rows(lapply(selected_vars, function(v) run_univar(data, v)))
write.csv(uni, file.path(OUT2, "41_svn_univariable_tempcontinuous.csv"), row.names = FALSE)
cat("  Saved: 41_svn_univariable_tempcontinuous.csv\n")

# ==============================================================================
# MULTIVARIABLE MODEL (all predictors together; "Not recorded" kept as a level
# for NOTREC_VARS; temp_val is complete-case like any other low-missingness
# covariate)
# ==============================================================================
miss_rows <- lapply(selected_vars, function(v) {
    n_miss <- sum(is.na(data[[v]]))
    data.frame(variable = v, n_total = nrow(data), n_missing = n_miss,
               pct_missing = round(n_miss / nrow(data) * 100, 1),
               stringsAsFactors = FALSE)
})
miss_summary <- do.call(rbind, miss_rows)
write.csv(miss_summary, file.path(OUT2, "41_svn_multivariable_tempcontinuous_missingness.csv"),
          row.names = FALSE)

add_not_recorded <- function(x) {
    if (is.factor(x)) { x <- addNA(x); levels(x)[is.na(levels(x))] <- "Not recorded" }
    x
}
data <- data %>% mutate(across(all_of(NOTREC_VARS), add_not_recorded))

n_cc <- sum(complete.cases(data[, selected_vars, drop = FALSE]))
cat(sprintf("Analysis sample after missing-data handling: %d / %d (%.1f%%)\n",
            n_cc, nrow(data), n_cc / nrow(data) * 100))
cat("  ('Not recorded' kept for: ", paste(NOTREC_VARS, collapse = ", "),
    "; other variables, including temp_val, complete-case.)\n", sep = "")

run_multivar <- function(df, vars, outcome_var = "outcome_bin", label = "Overall") {
    df2 <- df %>%
        select(all_of(c(outcome_var, vars))) %>%
        filter(!is.na(.data[[outcome_var]])) %>%
        tidyr::drop_na()
    cat(sprintf("  %s: %d complete cases, %d events\n",
                label, nrow(df2), sum(df2[[outcome_var]] == 1)))
    if (nrow(df2) < 20 || length(unique(df2[[outcome_var]])) < 2) return(NULL)

    active_vars <- vars[vapply(vars, function(v) {
        col <- df2[[v]]
        if (is.factor(col)) length(levels(droplevels(col))) >= 2
        else length(unique(col[!is.na(col)])) >= 2
    }, logical(1))]
    dropped <- setdiff(vars, active_vars)
    if (length(dropped) > 0)
        cat(sprintf("  Dropped (single level): %s\n", paste(dropped, collapse = ", ")))
    if (length(active_vars) == 0) return(NULL)

    fmla <- as.formula(paste(outcome_var, "~", paste(active_vars, collapse = " + ")))
    fit  <- glm(fmla, data = df2, family = binomial(link = "logit"))
    tidied <- tryCatch(
        broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE),
        error = function(e) {
            cat(sprintf("  Profiled CI failed (%s); using Wald CIs.\n", e$message))
            co <- summary(fit)$coefficients
            ci <- suppressMessages(confint.default(fit))
            tibble::tibble(term = rownames(co), estimate = exp(co[, "Estimate"]),
                           std.error = co[, "Std. Error"], statistic = co[, "z value"],
                           p.value = co[, "Pr(>|z|)"],
                           conf.low = exp(ci[, 1]), conf.high = exp(ci[, 2]))
        })
    tidied %>%
        filter(term != "(Intercept)") %>%
        mutate(model_group = label, n_total = nrow(df2),
               n_events = sum(df2[[outcome_var]] == 1),
               OR = round(estimate, 3), ci_lower = round(conf.low, 3),
               ci_upper = round(conf.high, 3), p_value = round(p.value, 4),
               OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, ci_lower, ci_upper),
               # For the continuous temperature term, also express the OR per
               # 1 degC DECREASE (the clinically-relevant direction: falling
               # temperature), by inverting the point estimate and swapping +
               # inverting the CI bounds.
               OR_per_1degC_decrease = ifelse(term == "temp_val", round(1 / estimate, 3), NA_real_),
               ci_decrease_lower = ifelse(term == "temp_val", round(1 / conf.high, 3), NA_real_),
               ci_decrease_upper = ifelse(term == "temp_val", round(1 / conf.low, 3), NA_real_)) %>%
        select(model_group, term, n_total, n_events,
               OR, ci_lower, ci_upper, p_value, OR_CI,
               OR_per_1degC_decrease, ci_decrease_lower, ci_decrease_upper)
}

cat("\nFitting multivariable model (SVN population, overall, continuous temperature)...\n")
m1 <- run_multivar(data, selected_vars, label = "SVN population (overall, temp continuous)")

if (!is.null(m1)) {
    write.csv(m1, file.path(OUT2, "41_svn_multivariable_tempcontinuous.csv"), row.names = FALSE)
    cat("  Saved: 41_svn_multivariable_tempcontinuous.csv\n")

    wb <- createWorkbook()
    addWorksheet(wb, "Multivariable")
    writeData(wb, "Multivariable", m1)
    addWorksheet(wb, "Univariable")
    writeData(wb, "Univariable", uni)
    addWorksheet(wb, "Covariates")
    writeData(wb, "Covariates", data.frame(covariate = selected_vars))
    addWorksheet(wb, "Missingness")
    writeData(wb, "Missingness", miss_summary)
    addWorksheet(wb, "Population_counts")
    writeData(wb, "Population_counts", count_tbl)
    saveWorkbook(wb, file.path(OUT2, "41_svn_multivariable_tempcontinuous.xlsx"), overwrite = TRUE)
    cat("  Saved: 41_svn_multivariable_tempcontinuous.xlsx\n")

    cat("\nKey multivariable ORs (SVN population, continuous temperature):\n")
    print(m1 %>% select(term, OR_CI, p_value))
    temp_row <- m1 %>% filter(term == "temp_val")
    if (nrow(temp_row) == 1) {
        cat(sprintf("\nAdjusted temperature effect: OR %.2f per 1 degC INCREASE (%.2f-%.2f), p%s\n",
                    temp_row$OR, temp_row$ci_lower, temp_row$ci_upper,
                    ifelse(temp_row$p_value < 0.001, "<0.001", sprintf("=%.4f", temp_row$p_value))))
        cat(sprintf("Equivalently: OR %.2f per 1 degC DECREASE (%.2f-%.2f)\n",
                    temp_row$OR_per_1degC_decrease, temp_row$ci_decrease_lower, temp_row$ci_decrease_upper))
    }
}

cat("\n41_svn_multivariable_tempcontinuous_regression.R complete.\n")
