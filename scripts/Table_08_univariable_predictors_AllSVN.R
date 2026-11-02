# ==============================================================================
# Purpose: Table 8 -- univariable risk factors for neonatal death, comparing
# the pooled small-vulnerable-newborn (All-SVN: Term-SGA, Preterm-AGA,
# Preterm-SGA) population against the Term-AGA reference group.
#
# Definitions: univariable only -- each factor is fitted in its own simple
# model (death ~ factor), since several predictors are strongly
# inter-related (Apgar, resuscitation, temperature) and a multivariable
# model would over-adjust and dilute the estimates for this descriptive
# question. The univariable set is run separately in the All-SVN and
# Term-AGA groups, with a factor-by-group interaction test for whether a
# risk factor behaves differently in vulnerable vs non-vulnerable newborns.
# Birth-weight category and gestation band are excluded as predictors
# (definitional components of the phenotype). Headline estimate:
# death ~ svn_group (All-SVN vs Term-AGA), univariable.
#
# Output: table_08_group_population_counts.csv, table_08_svn_group_headline.csv,
#         table_08_univariable_svn.csv, table_08_univariable_termaga.csv,
#         table_08_interaction_pvalues.csv, table_08_missingness.csv,
#         table_08_univariable_comparison.csv, 42_univariable.xlsx
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

cat("\n=== 42_svn_grouped_univariable.R ===\n")

# OUT2 defined in the setup block above (this folder's 'outputs')
dir.create(OUT2, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# LOAD DATA AND BUILD THE TWO GROUPS
# ==============================================================================
cat("Loading NNU data...\n")
raw <- load_nnu_data()
nnu <- process_nnu_data(raw)

PHENOTYPES <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")
SVN3       <- c("Term-SGA", "Preterm-AGA", "Preterm-SGA")   # the vulnerable group

data <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin)) %>%
    filter(svn_category %in% PHENOTYPES) %>%
    mutate(
        svn_category = factor(as.character(svn_category), levels = PHENOTYPES),
        # Binary grouping: Term-AGA (reference) vs the pooled SVN phenotypes
        svn_group = factor(
            ifelse(as.character(svn_category) == "Term-AGA", "Term-AGA", "All-SVN"),
            levels = c("Term-AGA", "All-SVN")),
        facility = factor(as.character(facility), levels = c("SMCH", "KCH"))
    )

# ==============================================================================
# DERIVE RISK-FACTOR VARIABLES
# Coding is identical to script 40_ (Version A categorisations) so the estimates
# are directly comparable with the rest of the pipeline. Birth-weight and
# gestation-band predictors are intentionally NOT derived.
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

        # Parity (Version A)
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

        # ANC (Version A: WHO threshold 4 visits)
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

        # Mode of delivery (Version A)
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

        # Apgar 5 min (Version A: 3 groups)
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

        # Temperature (categorical; reference Normal 36.5-37.5 degC)
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

# Risk factors to run one at a time (site included so it can also be reported)
RISK_FACTORS <- c("facility", "hiv_status", "parity", "matage_cat", "anc",
                  "typebirth_cat", "mode_del", "resus_cat", "inout_cat",
                  "apgar5_cat", "sex_cat", "temp_cat")

# ==============================================================================
# POPULATION COUNTS (denominators and events by group and site)
# ==============================================================================
cnt <- data %>%
    group_by(svn_group, facility) %>%
    summarise(n = n(), deaths = sum(outcome_bin == 1),
              cfr_per1000 = round(1000 * mean(outcome_bin == 1), 1),
              .groups = "drop")
cnt_grp <- data %>%
    group_by(svn_group) %>%
    summarise(facility = "BOTH", n = n(), deaths = sum(outcome_bin == 1),
              cfr_per1000 = round(1000 * mean(outcome_bin == 1), 1),
              .groups = "drop")
cnt_all <- data %>%
    summarise(svn_group = "TOTAL", facility = "BOTH", n = n(),
              deaths = sum(outcome_bin == 1),
              cfr_per1000 = round(1000 * mean(outcome_bin == 1), 1))
count_tbl <- bind_rows(cnt, cnt_grp, cnt_all)
write.csv(count_tbl, file.path(OUT2, "table_08_group_population_counts.csv"),
          row.names = FALSE)

cat(sprintf("All-SVN : n=%d, deaths=%d\n",
            sum(data$svn_group == "All-SVN"),
            sum(data$outcome_bin[data$svn_group == "All-SVN"] == 1)))
cat(sprintf("Term-AGA: n=%d, deaths=%d\n",
            sum(data$svn_group == "Term-AGA"),
            sum(data$outcome_bin[data$svn_group == "Term-AGA"] == 1)))
cat("  (See table_08_group_population_counts.csv for the full breakdown.)\n")

# ==============================================================================
# HEADLINE: death ~ svn_group  (All-SVN vs Term-AGA, univariable)
# ==============================================================================
fit_grp <- glm(outcome_bin ~ svn_group, data = data,
               family = binomial(link = "logit"))
grp_res <- broom::tidy(fit_grp, conf.int = TRUE, exponentiate = TRUE) %>%
    filter(term != "(Intercept)") %>%
    transmute(term = "All-SVN vs Term-AGA",
              OR = round(estimate, 3),
              ci_lower = round(conf.low, 3),
              ci_upper = round(conf.high, 3),
              p_value = signif(p.value, 3),
              OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, ci_lower, ci_upper))
write.csv(grp_res, file.path(OUT2, "table_08_svn_group_headline.csv"), row.names = FALSE)
cat("\nHeadline (All-SVN vs Term-AGA): ", grp_res$OR_CI,
    " p=", grp_res$p_value, "\n", sep = "")

# ==============================================================================
# UNIVARIABLE OR PER FACTOR, WITHIN A GIVEN GROUP
# ==============================================================================
run_univar <- function(df, predictor_var, outcome_var = "outcome_bin") {
    df2 <- df %>%
        select(all_of(c(outcome_var, predictor_var))) %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[predictor_var]])) %>%
        mutate(across(all_of(predictor_var), droplevels))
    ok <- nrow(df2) >= 10 &&
          length(unique(df2[[outcome_var]])) == 2 &&
          is.factor(df2[[predictor_var]]) &&
          nlevels(df2[[predictor_var]]) >= 2
    if (!ok) {
        return(tibble(predictor = predictor_var, term = NA_character_,
                      level = NA_character_, n_total = nrow(df2),
                      n_events = sum(df2[[outcome_var]] == 1),
                      OR = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                      p_value = NA_real_, note = "Insufficient data"))
    }
    fmla <- as.formula(paste(outcome_var, "~", predictor_var))
    fit  <- glm(fmla, data = df2, family = binomial(link = "logit"))
    res  <- broom::tidy(fit, conf.int = TRUE, exponentiate = TRUE)
    ref_level <- levels(df2[[predictor_var]])[1]
    ref_row <- tibble(term = paste0(predictor_var, ref_level),
                      estimate = 1.0, conf.low = NA_real_,
                      conf.high = NA_real_, p.value = NA_real_)
    res %>%
        filter(term != "(Intercept)") %>%
        bind_rows(ref_row) %>%
        mutate(predictor = predictor_var,
               level = sub(paste0("^", predictor_var), "", term),
               n_total = nrow(df2),
               n_events = sum(df2[[outcome_var]] == 1),
               OR = round(estimate, 3), ci_lower = round(conf.low, 3),
               ci_upper = round(conf.high, 3), p_value = round(p.value, 4),
               note = ifelse(term == paste0(predictor_var, ref_level),
                             "Reference", "")) %>%
        select(predictor, term, level, n_total, n_events,
               OR, ci_lower, ci_upper, p_value, note)
}

cat("\nRunning univariable regressions (All-SVN group)...\n")
svn_df  <- data %>% filter(svn_group == "All-SVN")
uni_svn <- bind_rows(lapply(RISK_FACTORS, function(v) run_univar(svn_df, v)))
uni_svn$group <- "All-SVN"
write.csv(uni_svn, file.path(OUT2, "table_08_univariable_svn.csv"), row.names = FALSE)

cat("Running univariable regressions (Term-AGA group)...\n")
ref_df  <- data %>% filter(svn_group == "Term-AGA")
uni_ref <- bind_rows(lapply(RISK_FACTORS, function(v) run_univar(ref_df, v)))
uni_ref$group <- "Term-AGA"
write.csv(uni_ref, file.path(OUT2, "table_08_univariable_termaga.csv"), row.names = FALSE)

# ==============================================================================
# INTERACTION TEST: does each factor's effect differ between the two groups?
# For every non-reference level we fit  death ~ factor * svn_group  and take the
# likelihood-ratio p-value for the factor:group interaction block (compared with
# the additive model death ~ factor + svn_group). One p-value per factor.
# ==============================================================================
run_interaction <- function(df, predictor_var, outcome_var = "outcome_bin") {
    df2 <- df %>%
        select(all_of(c(outcome_var, predictor_var, "svn_group"))) %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[predictor_var]]),
               !is.na(svn_group)) %>%
        mutate(across(all_of(predictor_var), droplevels),
               svn_group = droplevels(svn_group))
    ok <- nrow(df2) >= 20 &&
          is.factor(df2[[predictor_var]]) &&
          nlevels(df2[[predictor_var]]) >= 2 &&
          nlevels(df2$svn_group) == 2
    if (!ok) return(tibble(predictor = predictor_var,
                           p_interaction = NA_real_, note = "Not estimable"))
    add <- glm(as.formula(paste(outcome_var, "~", predictor_var, "+ svn_group")),
               data = df2, family = binomial(link = "logit"))
    inx <- glm(as.formula(paste(outcome_var, "~", predictor_var, "* svn_group")),
               data = df2, family = binomial(link = "logit"))
    lrt <- anova(add, inx, test = "LRT")
    p   <- lrt$`Pr(>Chi)`[2]
    tibble(predictor = predictor_var, p_interaction = round(p, 4),
           note = ifelse(is.na(p), "Not estimable", ""))
}

cat("Testing SVN vs Term-AGA interaction per factor...\n")
inter <- bind_rows(lapply(RISK_FACTORS, function(v) run_interaction(data, v)))
write.csv(inter, file.path(OUT2, "table_08_interaction_pvalues.csv"), row.names = FALSE)

# ==============================================================================
# SIDE-BY-SIDE COMPARISON TABLE (SVN | Term-AGA | p_interaction)
# ==============================================================================
fmt <- function(df) df %>%
    mutate(OR_CI = ifelse(note == "Reference", "1 (ref)",
                   ifelse(is.na(OR), "-",
                          sprintf("%.2f (%.2f-%.2f)", OR, ci_lower, ci_upper)))) %>%
    select(predictor, term, level, OR_CI, p_value)

comp <- fmt(uni_svn) %>%
    rename(OR_CI_SVN = OR_CI, p_SVN = p_value) %>%
    left_join(fmt(uni_ref) %>%
                  rename(OR_CI_TermAGA = OR_CI, p_TermAGA = p_value),
              by = c("predictor", "term", "level")) %>%
    left_join(inter, by = "predictor") %>%
    select(predictor, level, OR_CI_SVN, p_SVN,
           OR_CI_TermAGA, p_TermAGA, p_interaction)
write.csv(comp, file.path(OUT2, "table_08_univariable_comparison.csv"), row.names = FALSE)

# ==============================================================================
# MISSINGNESS PER FACTOR BY GROUP
# ==============================================================================
miss <- bind_rows(lapply(RISK_FACTORS, function(v) {
    tibble(
        variable = v,
        svn_missing     = sum(is.na(svn_df[[v]])),
        svn_pct         = round(100 * mean(is.na(svn_df[[v]])), 1),
        termaga_missing = sum(is.na(ref_df[[v]])),
        termaga_pct     = round(100 * mean(is.na(ref_df[[v]])), 1))
}))
write.csv(miss, file.path(OUT2, "table_08_missingness.csv"), row.names = FALSE)

# ==============================================================================
# EXCEL WORKBOOK (all sheets in one file for the manuscript tables)
# ==============================================================================
wb <- createWorkbook()
addWorksheet(wb, "Population_counts"); writeData(wb, "Population_counts", count_tbl)
addWorksheet(wb, "Headline_group");    writeData(wb, "Headline_group", grp_res)
addWorksheet(wb, "Univariable_SVN");   writeData(wb, "Univariable_SVN", uni_svn)
addWorksheet(wb, "Univariable_TermAGA"); writeData(wb, "Univariable_TermAGA", uni_ref)
addWorksheet(wb, "Comparison");        writeData(wb, "Comparison", comp)
addWorksheet(wb, "Interaction_p");     writeData(wb, "Interaction_p", inter)
addWorksheet(wb, "Missingness");       writeData(wb, "Missingness", miss)
addWorksheet(wb, "Risk_factors");      writeData(wb, "Risk_factors",
                                                  data.frame(risk_factor = RISK_FACTORS))
saveWorkbook(wb, file.path(OUT2, "42_univariable.xlsx"), overwrite = TRUE)

cat("\nSaved to ", OUT2, ":\n", sep = "")
cat("  table_08_group_population_counts.csv\n")
cat("  table_08_svn_group_headline.csv\n")
cat("  table_08_univariable_svn.csv\n")
cat("  table_08_univariable_termaga.csv\n")
cat("  table_08_univariable_comparison.csv\n")
cat("  table_08_interaction_pvalues.csv\n")
cat("  table_08_missingness.csv\n")
cat("  42_univariable.xlsx\n")

cat("\nUnivariable ORs, All-SVN group (headline factors):\n")
print(uni_svn %>%
      filter(note != "Reference", !is.na(OR)) %>%
      mutate(OR_CI = sprintf("%.2f (%.2f-%.2f)", OR, ci_lower, ci_upper)) %>%
      select(predictor, level, n_total, n_events, OR_CI, p_value),
      n = 100)

cat("\n42_svn_grouped_univariable.R complete.\n")
