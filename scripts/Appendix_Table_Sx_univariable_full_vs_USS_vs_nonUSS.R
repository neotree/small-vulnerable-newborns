# ==============================================================================
# Purpose: compares the univariable SVN risk-factor results between the full
# cohort and the ultrasound-dated (USS) subset, to test whether ultrasound
# dating changes the associations.
#
# Definitions: the USS subset is nested within the full cohort, so a direct
# "full vs USS" contrast is not a comparison of independent samples. Instead
# this partitions the cohort into USS (methodestgest in USS/EUSS) and non-USS
# (all other GA methods) -- two disjoint groups -- and reports the univariable
# OR for every factor in all three views (ALL, USS, non-USS), plus a
# likelihood-ratio interaction p-value per factor testing whether the
# association differs between USS-dated and non-USS-dated newborns. A
# non-significant p means the estimates are robust to GA dating method.
# All estimates are within the pooled All-SVN group (Term-SGA, Preterm-AGA,
# Preterm-SGA); Term-AGA is the comparator only for the headline estimate.
#
# Output: appendix_tablesx_group_counts_by_dating.csv,
#         appendix_tablesx_headline_full_vs_uss.csv,
#         appendix_tablesx_univariable_full_vs_uss.csv, 44_compare.xlsx
# ==============================================================================

POPULATION <- "all"
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

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(broom) })

cat("\n=== 44_compare_full_vs_uss.R ===\n")

# OUT2 defined in the setup block above (this folder's 'outputs')
dir.create(OUT2, recursive = TRUE, showWarnings = FALSE)

USS_METHODS <- c("USS", "EUSS")
PHENOTYPES  <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")

cat("Loading NNU data (full cohort)...\n")
raw <- load_nnu_data()
nnu <- process_nnu_data(raw)

data <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin)) %>%
    filter(svn_category %in% PHENOTYPES) %>%
    mutate(
        svn_category = factor(as.character(svn_category), levels = PHENOTYPES),
        svn_group = factor(
            ifelse(as.character(svn_category) == "Term-AGA", "Term-AGA", "All-SVN"),
            levels = c("Term-AGA", "All-SVN")),
        facility = factor(as.character(facility), levels = c("SMCH", "KCH")),
        uss_flag = factor(
            ifelse(!is.na(methodestgest) &
                   toupper(trimws(methodestgest)) %in% USS_METHODS,
                   "USS", "Non-USS"),
            levels = c("Non-USS", "USS"))
    )

# ==============================================================================
# DERIVE RISK-FACTOR VARIABLES (identical coding to Table_08_univariable_predictors_AllSVN.R)
# ==============================================================================
data <- data %>%
    mutate(
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
        par_num = suppressWarnings(as.numeric(par)),
        parity = factor(case_when(
            is.na(par_num) ~ NA_character_,
            par_num == 1   ~ "Primiparous (1st delivery)",
            par_num <= 4   ~ "Multiparous (2nd-4th)",
            TRUE           ~ "Grand Multipara (5th+)"
        ), levels = c("Primiparous (1st delivery)",
                      "Multiparous (2nd-4th)","Grand Multipara (5th+)")),
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
        anc_num = suppressWarnings(as.numeric(antenatalcare)),
        anc = anc_factor(antenatalcare, antenatalcare_category),
        typebirth_cat = case_when(
            is.na(typebirth) | trimws(typebirth) == "" ~ NA_character_,
            toupper(trimws(typebirth)) %in% c("S","SINGLETON") ~ "Singleton",
            toupper(trimws(typebirth)) %in% c("TW1","TW2","TWIN1","TWIN2","TWIN") ~ "Twin",
            toupper(trimws(typebirth)) %in% c("TR1","TR2","TR3","TRIPLET") ~ "Triplet",
            TRUE ~ NA_character_
        ),
        typebirth_cat = factor(typebirth_cat, levels = c("Singleton","Twin","Triplet")),
        modedelivery_num = suppressWarnings(as.numeric(modedelivery)),
        mode_del = factor(case_when(
            is.na(modedelivery_num) ~ NA_character_,
            modedelivery_num == 1   ~ "SVD",
            modedelivery_num %in% c(2,3,6,7) ~ "Assisted Vaginal",
            modedelivery_num == 4   ~ "ELCS",
            modedelivery_num %in% c(5,8) ~ "EMCS",
            TRUE ~ NA_character_
        ), levels = c("SVD","Assisted Vaginal","ELCS","EMCS")),
        resus_cat = resus_category(resus),
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
        apgar5_num = suppressWarnings(as.numeric(apgar5)),
        apgar5_num = ifelse(apgar5_num >= 0 & apgar5_num <= 10, apgar5_num, NA_real_),
        apgar5_cat = factor(case_when(
            is.na(apgar5_num)  ~ NA_character_,
            apgar5_num <= 3    ~ "Severely depressed (0-3)",
            apgar5_num <= 6    ~ "Moderately depressed (4-6)",
            TRUE               ~ "Normal (7-10)"
        ), levels = c("Normal (7-10)","Moderately depressed (4-6)","Severely depressed (0-3)")),
        sex_cat = case_when(
            toupper(trimws(as.character(gender_std))) %in% c("M","MALE","BOY")    ~ "Male",
            toupper(trimws(as.character(gender_std))) %in% c("F","FEMALE","GIRL") ~ "Female",
            TRUE ~ NA_character_
        ),
        sex_cat = factor(sex_cat, levels = c("Male","Female")),
        temp_val = suppressWarnings(as.numeric(temperature)),
        temp_val = ifelse(temp_val >= 30 & temp_val <= 45, temp_val, NA_real_),
        temp_cat = case_when(
            is.na(temp_val)    ~ NA_character_,
            temp_val <  32.0   ~ "Severe Hypothermia (<32 degC)",
            temp_val <  36.0   ~ "Moderate Hypothermia (32-35.9 degC)",
            temp_val <  36.5   ~ "Mild Hypothermia (36.0-36.4 degC)",
            temp_val <= 37.5   ~ "Normal (36.5-37.5 degC)",
            TRUE               ~ "Hyperthermia (>37.5 degC)"
        ),
        temp_cat = factor(temp_cat,
                          levels = c("Normal (36.5-37.5 degC)",
                                     "Severe Hypothermia (<32 degC)",
                                     "Moderate Hypothermia (32-35.9 degC)",
                                     "Mild Hypothermia (36.0-36.4 degC)",
                                     "Hyperthermia (>37.5 degC)"))
    )

RISK_FACTORS <- c("facility", "hiv_status", "parity", "matage_cat", "anc",
                  "typebirth_cat", "mode_del", "resus_cat", "inout_cat",
                  "apgar5_cat", "sex_cat", "temp_cat")

# ==============================================================================
# COUNTS BY GROUP x DATING METHOD
# ==============================================================================
counts <- data %>%
    group_by(svn_group, uss_flag) %>%
    summarise(n = n(), deaths = sum(outcome_bin == 1),
              cfr_per1000 = round(1000 * mean(outcome_bin == 1), 1),
              .groups = "drop")
write.csv(counts, file.path(OUT2, "appendix_tablesx_group_counts_by_dating.csv"), row.names = FALSE)
cat("\nCounts by group x dating method:\n"); print(counts)

# ==============================================================================
# UNIVARIABLE OR PER FACTOR WITHIN A DATA SUBSET (returns per-level rows)
# ==============================================================================
univ_or <- function(df, predictor_var, label, outcome_var = "outcome_bin") {
    df2 <- df %>%
        select(all_of(c(outcome_var, predictor_var))) %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[predictor_var]])) %>%
        mutate(across(all_of(predictor_var), droplevels))
    ok <- nrow(df2) >= 10 && length(unique(df2[[outcome_var]])) == 2 &&
          is.factor(df2[[predictor_var]]) && nlevels(df2[[predictor_var]]) >= 2
    ref_level <- if (is.factor(df[[predictor_var]])) levels(df[[predictor_var]])[1] else NA
    if (!ok) {
        return(tibble(predictor = predictor_var, level = NA_character_,
                      grp = label, OR = NA_real_, ci_lower = NA_real_,
                      ci_upper = NA_real_, p_value = NA_real_,
                      n_total = nrow(df2),
                      n_events = sum(df2[[outcome_var]] == 1, na.rm = TRUE)))
    }
    fit <- glm(as.formula(paste(outcome_var, "~", predictor_var)),
               data = df2, family = binomial(link = "logit"))
    res <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)")
    ref <- tibble(term = paste0(predictor_var, levels(df2[[predictor_var]])[1]),
                  estimate = 1, conf.low = NA_real_, conf.high = NA_real_,
                  p.value = NA_real_)
    bind_rows(res, ref) %>%
        transmute(predictor = predictor_var,
                  level = sub(paste0("^", predictor_var), "", term),
                  grp = label,
                  OR = round(estimate, 3), ci_lower = round(conf.low, 3),
                  ci_upper = round(conf.high, 3), p_value = round(p.value, 4),
                  n_total = nrow(df2),
                  n_events = sum(df2[[outcome_var]] == 1))
}

or_ci_str <- function(OR, lo, hi, p) {
    ifelse(is.na(OR), "-",
    ifelse(is.na(lo) & OR == 1, "1 (ref)",
           sprintf("%.2f (%.2f-%.2f)", OR, lo, hi)))
}

# All-SVN subsets by dating method
svn_all    <- data %>% filter(svn_group == "All-SVN")
svn_uss    <- svn_all %>% filter(uss_flag == "USS")
svn_nonuss <- svn_all %>% filter(uss_flag == "Non-USS")

build_block <- function(df, lab)
    bind_rows(lapply(RISK_FACTORS, function(v) univ_or(df, v, lab)))

or_all    <- build_block(svn_all,    "ALL")
or_uss    <- build_block(svn_uss,    "USS")
or_nonuss <- build_block(svn_nonuss, "Non-USS")

mk <- function(df) df %>%
    mutate(OR_CI = or_ci_str(OR, ci_lower, ci_upper, p_value)) %>%
    select(predictor, level, OR_CI, p_value)

# ==============================================================================
# HETEROGENEITY: death ~ factor * uss_flag within All-SVN (LRT interaction)
# Tests whether each factor's association differs between USS and Non-USS.
# ==============================================================================
het_test <- function(df, predictor_var, outcome_var = "outcome_bin") {
    df2 <- df %>%
        select(all_of(c(outcome_var, predictor_var, "uss_flag"))) %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[predictor_var]]),
               !is.na(uss_flag)) %>%
        mutate(across(all_of(predictor_var), droplevels),
               uss_flag = droplevels(uss_flag))
    ok <- nrow(df2) >= 20 && is.factor(df2[[predictor_var]]) &&
          nlevels(df2[[predictor_var]]) >= 2 && nlevels(df2$uss_flag) == 2
    if (!ok) return(tibble(predictor = predictor_var,
                           p_interaction = NA_real_, note = "Not estimable"))
    out <- tryCatch({
        add <- glm(as.formula(paste(outcome_var, "~", predictor_var, "+ uss_flag")),
                   data = df2, family = binomial(link = "logit"))
        inx <- glm(as.formula(paste(outcome_var, "~", predictor_var, "* uss_flag")),
                   data = df2, family = binomial(link = "logit"))
        p <- anova(add, inx, test = "LRT")$`Pr(>Chi)`[2]
        tibble(predictor = predictor_var, p_interaction = round(p, 4),
               note = ifelse(is.na(p), "Not estimable", ""))
    }, error = function(e) tibble(predictor = predictor_var,
                                  p_interaction = NA_real_, note = "Model failed"))
    out
}

het <- bind_rows(lapply(RISK_FACTORS, function(v) het_test(svn_all, v)))

# ==============================================================================
# ASSEMBLE THE SIDE-BY-SIDE COMPARISON TABLE
# ==============================================================================
comp <- mk(or_all)    %>% rename(OR_ALL = OR_CI,     p_ALL = p_value) %>%
    left_join(mk(or_uss)    %>% rename(OR_USS = OR_CI,     p_USS = p_value),
              by = c("predictor", "level")) %>%
    left_join(mk(or_nonuss) %>% rename(OR_NonUSS = OR_CI,  p_NonUSS = p_value),
              by = c("predictor", "level")) %>%
    left_join(het %>% select(predictor, p_interaction), by = "predictor") %>%
    select(predictor, level, OR_ALL, OR_USS, OR_NonUSS,
           p_ALL, p_USS, p_NonUSS, p_interaction)
write.csv(comp, file.path(OUT2, "appendix_tablesx_univariable_full_vs_uss.csv"), row.names = FALSE)

# ==============================================================================
# HEADLINE: svn_group (All-SVN vs Term-AGA) in ALL / USS / Non-USS + interaction
# ==============================================================================
headline_or <- function(df, lab) {
    fit <- glm(outcome_bin ~ svn_group, data = df, family = binomial(link = "logit"))
    broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE) %>%
        filter(term != "(Intercept)") %>%
        transmute(grp = lab, OR = round(estimate, 3),
                  ci_lower = round(conf.low, 3), ci_upper = round(conf.high, 3),
                  p_value = signif(p.value, 3),
                  OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, ci_lower, ci_upper))
}
hl <- bind_rows(headline_or(data, "ALL"),
                headline_or(filter(data, uss_flag == "USS"),     "USS"),
                headline_or(filter(data, uss_flag == "Non-USS"), "Non-USS"))
# interaction: does the SVN-vs-Term-AGA OR differ by dating method?
hl_add <- glm(outcome_bin ~ svn_group + uss_flag, data = data,
              family = binomial(link = "logit"))
hl_inx <- glm(outcome_bin ~ svn_group * uss_flag, data = data,
              family = binomial(link = "logit"))
hl_p <- anova(hl_add, hl_inx, test = "LRT")$`Pr(>Chi)`[2]
hl$p_interaction <- c(round(hl_p, 4), NA, NA)
write.csv(hl, file.path(OUT2, "appendix_tablesx_headline_full_vs_uss.csv"), row.names = FALSE)

cat(sprintf("\nHeadline All-SVN vs Term-AGA: ALL %s | USS %s | Non-USS %s | interaction p=%s\n",
            hl$OR_CI[1], hl$OR_CI[2], hl$OR_CI[3], format(round(hl_p, 4))))

# ==============================================================================
# EXCEL WORKBOOK
# ==============================================================================
wb <- createWorkbook()
addWorksheet(wb, "Counts");        writeData(wb, "Counts", counts)
addWorksheet(wb, "Headline");      writeData(wb, "Headline", hl)
addWorksheet(wb, "Univariable");   writeData(wb, "Univariable", comp)
addWorksheet(wb, "Heterogeneity"); writeData(wb, "Heterogeneity", het)
saveWorkbook(wb, file.path(OUT2, "44_compare.xlsx"), overwrite = TRUE)

cat("\nSaved to ", OUT2, " (44_ prefix):\n", sep = "")
cat("  appendix_tablesx_group_counts_by_dating.csv\n  appendix_tablesx_headline_full_vs_uss.csv\n")
cat("  appendix_tablesx_univariable_full_vs_uss.csv\n  44_compare.xlsx\n")

n_sig <- sum(het$p_interaction < 0.05, na.rm = TRUE)
cat(sprintf("\nHeterogeneity summary: %d of %d factors show a significant USS vs Non-USS\n",
            n_sig, sum(!is.na(het$p_interaction))))
cat("difference in effect (p_interaction < 0.05). Factors with no significant\n")
cat("difference are consistent between the full cohort and the USS subset.\n")
print(het)

cat("\n44_compare_full_vs_uss.R complete.\n")
