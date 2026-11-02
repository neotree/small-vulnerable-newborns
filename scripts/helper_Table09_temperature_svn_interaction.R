# ==============================================================================
# Purpose: tests whether the per-1degC admission-temperature effect on
# neonatal death differs across SVN phenotypes -- the interaction underlying
# the "steeper in some phenotypes than others" claim that accompanies the
# per-degree odds ratios in Supplementary_Table_S1_validation_CFR_full_vs_USS.R.
#
# Definitions: one overall likelihood-ratio test (not per-phenotype), minimally
# adjusted for facility:
#   additive:     outcome ~ temp_val + svn_category + facility
#   interaction:  outcome ~ temp_val * svn_category + facility
# temp_val construction (valid range, coercion) is identical to
# dep_Table09_temperature_categorical_OR.R and
# Supplementary_Table_S1_validation_CFR_full_vs_USS.R. Population =
# svn_filtered, all four SVN phenotypes, both sites combined. Respects
# POPULATION ("all" default, "USS" for the ultrasound subset).
#
# Output: helper_table09_temp_svn_interaction.csv
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
cat("\n=== 57_temp_svn_interaction.R ===\n")
TEMP_VALID_MIN <- 30.0
TEMP_VALID_MAX <- 45.0

nnu <- process_nnu_data(load_nnu_data())
base <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin)) %>%
    mutate(
        facility     = factor(as.character(facility), levels = c("KCH","SMCH")),
        svn_category = factor(as.character(svn_category), levels = SVN_LEVELS),
        temp_val = suppressWarnings(as.numeric(temperature)),
        temp_val = ifelse(temp_val >= TEMP_VALID_MIN & temp_val <= TEMP_VALID_MAX,
                          temp_val, NA_real_)
    ) %>%
    filter(!is.na(temp_val), !is.na(svn_category))

cat(sprintf("Interaction-test population: %d (deaths %d)\n", nrow(base), sum(base$outcome_bin)))

intx <- tryCatch({
    fit_add <- glm(outcome_bin ~ temp_val + svn_category + facility,
                   data = base, family = binomial())
    fit_int <- glm(outcome_bin ~ temp_val * svn_category + facility,
                   data = base, family = binomial())
    an    <- anova(fit_add, fit_int, test = "Chisq")
    p_int <- an$`Pr(>Chi)`[2]
    tibble(
        test            = "Temperature (per 1 degC) x SVN-phenotype interaction",
        n               = nrow(base),
        n_events        = sum(base$outcome_bin),
        df              = an$Df[2],
        lr_chisq        = round(an$Deviance[2], 3),
        p_value         = signif(p_int, 4),
        p_fmt           = ifelse(p_int < 0.001, "<0.001", sprintf("%.3f", p_int)),
        additive_aic    = round(AIC(fit_add), 1),
        interaction_aic = round(AIC(fit_int), 1),
        note            = ifelse(p_int < 0.05,
                                 "Temperature effect differs significantly across SVN phenotypes (p < 0.05)",
                                 "No significant difference in the temperature effect across SVN phenotypes (p >= 0.05)")
    )
}, error = function(e) tibble(
    test = "Temperature (per 1 degC) x SVN-phenotype interaction", n = NA_integer_,
    n_events = NA_integer_, df = NA_integer_, lr_chisq = NA_real_, p_value = NA_real_,
    p_fmt = NA_character_, additive_aic = NA_real_, interaction_aic = NA_real_,
    note = paste("Error:", e$message)))

write.csv(intx, file.path(OUT2, "helper_table09_temp_svn_interaction.csv"), row.names = FALSE)

cat("\nTemperature x SVN-phenotype interaction (likelihood-ratio test):\n")
print(as.data.frame(intx), row.names = FALSE)
cat("\nSaved helper_table09_temp_svn_interaction.csv\n")
cat("57_temp_svn_interaction.R complete.\n")
