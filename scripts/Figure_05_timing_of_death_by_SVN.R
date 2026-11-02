# ==============================================================================
# Purpose: Figure 5 -- timing of neonatal death by SVN phenotype, both sites
# combined, single-panel line chart.
#
# Definitions: uses the same population and definitions as
# Table_05_timing_of_death.R, so the numbers match that table exactly.
# y = % of that phenotype's recorded-age deaths in each time band
# (unknown-age deaths excluded from the denominator). Time bands: 0-1, 1-3,
# 3-7, 7-14, 14-21, >21 days.
#
# Output: figure_05_timing_of_death_by_svn.png/.pdf
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

suppressPackageStartupMessages(library(lubridate))

cat("\n===== Timing of neonatal death by SVN phenotype (lines, combined) =====\n")

# ------------------------------------------------------------------------------
# DATA BUILD (identical to Table_05_timing_of_death.R)
# ------------------------------------------------------------------------------
.AGE_MAX_HOURS <- 672
.PH        <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")
.INTERVALS <- c("0-1 days", "1-3 days", "3-7 days", "7-14 days", "14-21 days", ">21 days")
.classify_age <- function(days) {
  dplyr::case_when(
    is.na(days) | days < 0 ~ NA_character_,
    days <= 1.0  ~ "0-1 days",
    days <= 3.0  ~ "1-3 days",
    days <= 7.0  ~ "3-7 days",
    days <= 14.0 ~ "7-14 days",
    days <= 21.0 ~ "14-21 days",
    TRUE         ~ ">21 days"
  )
}
.nnu <- process_nnu_data(load_nnu_data())
.deaths <- .nnu$svn_filtered %>%
  filter(!is.na(outcome_bin), outcome_bin == 1, svn_category %in% .PH) %>%
  mutate(
    category     = as.character(svn_category),
    .age_adm_h   = suppressWarnings(as.numeric(age)),
    .age_adm_h   = ifelse(is.finite(.age_adm_h) & .age_adm_h <= .AGE_MAX_HOURS, .age_adm_h, NA_real_),
    .adm_ts      = ymd_hms(datetimeadmission, quiet = TRUE),
    .dth_ts      = ymd_hms(datetimedeath, quiet = TRUE),
    .age_death_d = .age_adm_h / 24 + as.numeric(difftime(.dth_ts, .adm_ts, units = "days")),
    interval     = factor(.classify_age(.age_death_d), levels = .INTERVALS)
  ) %>%
  filter(!is.na(interval))            # recorded-age deaths only (Unknown excluded)

# combined-sites: % of each phenotype's recorded-age deaths in each band
dod <- .deaths %>%
  group_by(category, interval, .drop = FALSE) %>%
  summarise(n_deaths = n(), .groups = "drop") %>%
  group_by(category) %>%
  mutate(n_total_deaths = sum(n_deaths),
         pct_deaths     = n_deaths / n_total_deaths * 100) %>%
  ungroup() %>%
  mutate(category = factor(category, levels = .PH))

# save the plotted data for transparency
write.csv(dod, file.path(OUTPUT_DIR_NNU, "figure_05_timing_of_death_by_svn_data.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# PLOT - single-panel lines
# ------------------------------------------------------------------------------
band_labels <- c("0-1 days" = "0-1 d", "1-3 days" = "1-3 d", "3-7 days" = "3-7 d",
                 "7-14 days" = "7-14 d", "14-21 days" = "14-21 d", ">21 days" = ">21 d")

fig5 <- ggplot(dod, aes(x = interval, y = pct_deaths,
                        colour = category, group = category)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 3) +
  scale_colour_manual(values = SVN_COLOURS, name = NULL) +
  scale_x_discrete(labels = band_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  labs(
    title = "Timing of neonatal death by SVN phenotype (both sites combined)",
    x     = "Age at death (days)",
    y     = "% of deaths in group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    legend.position  = c(0.87, 0.80),
    legend.background = element_rect(fill = scales::alpha("white", 0.6), colour = NA)
  )

save_figure(fig5, "figure_05_timing_of_death_by_svn", width = 9, height = 5.2)

# ------------------------------------------------------------------------------
# VERIFY rendered values
# ------------------------------------------------------------------------------
cat("\nRendered % of deaths in group (recorded-age denominator):\n")
print(as.data.frame(dod %>%
  transmute(category, interval, n_deaths, pct = round(pct_deaths, 1))))
cat("Timing of death by SVN (single-panel): complete.\n")
