# ==============================================================================
# Purpose: Figure 7 -- categorical admission-temperature odds ratios by SVN
# group, reference 36.5-37.5 degC. Shows that hypothermia raises risk in
# every phenotype, while only Term-AGA also shows raised risk at
# hyperthermia (U-shape).
#
# Reads: dep_table09_temp_cat_or_by_svn.csv, dep_table09_temp_cat_sample_sizes.csv
#   (written by dep_Table09_temperature_categorical_OR.R -- run that first)
#
# The black reference line is "All-SVN" (pooled Term-SGA + Preterm-AGA +
# Preterm-SGA), the same column as Supplementary Table S6. All-SVN ORs are
# computed from the pooled counts as in helper_Table09_temperature_AllSVN_column.R
# (2x2 OR vs Normal = glm MLE).
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

cat("\n===== Temperature categorical OR by SVN group =====\n")

# --- All-SVN: pooled Term-SGA + Preterm-AGA + Preterm-SGA (Supplementary Table S6) --
ss <- read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_sample_sizes.csv"),
               show_col_types = FALSE) %>%
  filter(svn_category %in% c("Term-SGA", "Preterm-AGA", "Preterm-SGA")) %>%
  group_by(temp_cat) %>%
  summarise(n = sum(n), d = sum(n_deaths), .groups = "drop")
ref_d <- ss$d[ss$temp_cat == "4_Normal"]
ref_a <- ss$n[ss$temp_cat == "4_Normal"] - ref_d
allsvn <- ss %>%
  filter(temp_cat != "4_Normal") %>%
  mutate(a = n - d,
         OR = (d / a) / (ref_d / ref_a),
         se = sqrt(1 / d + 1 / a + 1 / ref_d + 1 / ref_a),
         ci_lower = exp(log(OR) - 1.959964 * se),
         ci_upper = exp(log(OR) + 1.959964 * se),
         p_value  = 2 * pnorm(-abs(log(OR) / se))) %>%
  transmute(group = "All-SVN", temp_cat, OR, ci_lower, ci_upper, p_value) %>%
  bind_rows(tibble(group = "All-SVN", temp_cat = "4_Normal (Reference)", OR = 1,
                   ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_))

ov <- allsvn %>% transmute(series = group, temp_cat, OR)

bg <- read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_by_svn.csv"),
               show_col_types = FALSE) %>%
  transmute(series = group, temp_cat, OR)

cat_levels <- c("1_SevereHypothermia", "2_Hypothermia", "3_LowNormal",
                "4_Normal (Reference)", "5_Hyperthermia")
cat_labels <- c("Severe\nhypothermia", "Hypothermia", "Mild\nhypothermia",
                "Normal\n(reference)", "Hyperthermia")

plot_df <- bind_rows(ov, bg) %>%
  filter(temp_cat %in% cat_levels, !is.na(OR)) %>%   # drop not-estimable cells
  mutate(
    temp_f = factor(temp_cat, levels = cat_levels, labels = cat_labels),
    series = factor(series, levels = c("All-SVN", SVN_LEVELS))
  )

series_cols <- c("All-SVN" = "#222222", SVN_COLOURS)

fig8 <- ggplot(plot_df, aes(x = temp_f, y = OR, colour = series, group = series)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = 1.1) +
  geom_point(size = 3) +
  scale_y_log10() +
  scale_colour_manual(values = series_cols, name = NULL) +
  labs(
    title    = "Temperature-mortality by SVN group",
    subtitle = "Odds of neonatal death vs normal range (36.5-37.5 \u00b0C); log scale. All-SVN = Term-SGA + Preterm-AGA + Preterm-SGA",
    x        = NULL,
    y        = "Odds ratio (log scale)",
    caption  = svn_caption("Reference: normal range. Values as in Supplementary Table S6.")
  ) +
  theme_svn(legend_position = "right")

save_figure(fig8, "figure_07_temperature_categorical_ors_by_svn", width = 9, height = 6)

# ------------------------------------------------------------------------------
# COMPANION VALUES TABLE: OR (95% CI) and Wald p vs the normal range, per group
# Written so the figure's underlying numbers and test are reproducible from R.
# ------------------------------------------------------------------------------
plain_lab <- c("1_SevereHypothermia" = "Severe hypothermia",
               "2_Hypothermia"        = "Hypothermia",
               "3_LowNormal"          = "Mild hypothermia",
               "4_Normal (Reference)" = "Normal (reference)",
               "5_Hyperthermia"       = "Hyperthermia")

tbl8 <- bind_rows(
  allsvn,
  read_csv(file.path(OUTPUT_DIR, "dep_table09_temp_cat_or_by_svn.csv"), show_col_types = FALSE)
) %>%
  filter(temp_cat %in% names(plain_lab)) %>%
  transmute(
    group                = factor(group, levels = c("All-SVN", SVN_LEVELS)),
    temperature_category = unname(plain_lab[temp_cat]),
    OR                   = round(OR, 3),
    ci_lower             = round(ci_lower, 3),
    ci_upper             = round(ci_upper, 3),
    p_value_vs_normal    = round(p_value, 4)
  ) %>%
  arrange(group)

write_csv(tbl8, file.path(OUTPUT_DIR, "figure_07_temperature_by_svn_table.csv"))
cat("  Saved: figure_07_temperature_by_svn_table.csv\n")
cat("Temperature categorical OR by SVN: complete.\n")
