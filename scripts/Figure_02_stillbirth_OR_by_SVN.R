# ==============================================================================
# Purpose: Figure 2 -- stillbirth odds ratio by SVN phenotype, both sites
# combined, reference = Term-AGA.
#
# Definitions: computed directly from the maternal data (all deliveries with
# valid birth weight and gestational age, 73,845), matching Table 1's
# population. Panel A = all stillbirths; Panel B = late stillbirths
# (>= 28 weeks).
#
# Output: figure_02_stillbirth_or_data.csv, figure_02_stillbirth_or_combined.png/.pdf
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
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(broom); library(ggplot2); library(patchwork)})
cat("\n=== Figure_02_stillbirth_OR_by_SVN.R ===\n")

STILLBIRTH_CODES <- c("SBF", "SBM", "STBF", "STBM")
GA_MIN_M <- if (exists("GA_MIN_MAT")) GA_MIN_MAT else 22
PH <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")

# ---- current maternal base: all deliveries with valid BW + GA (both sites) ----
base <- process_maternal_data(load_maternal_data()) %>%
    filter(is.finite(bw_kg), bw_kg >= BW_MIN_KG, bw_kg <= BW_MAX_KG,
           is.finite(gestation_weeks), gestation_weeks >= GA_MIN_M, gestation_weeks <= GA_MAX) %>%
    mutate(
        svn     = factor(make_svn_category(classify_sga(bw_kg, gestation_weeks, gender_std), gestation_weeks), levels = PH),
        sb      = as.integer(!is.na(neotreeoutcome) & toupper(trimws(neotreeoutcome)) %in% STILLBIRTH_CODES),
        late_sb = as.integer(sb == 1 & gestation_weeks >= 28)
    ) %>%
    filter(!is.na(svn))
cat(sprintf("Deliveries %d; stillbirths %d (late >=28w %d)\n",
            nrow(base), sum(base$sb), sum(base$late_sb)))

# ---- OR vs Term-AGA (univariable logistic; combined both sites) --------------
or_for <- function(df, outcome_col, scope_label) {
    df <- df %>% mutate(y = .data[[outcome_col]])
    df$svn <- relevel(df$svn, ref = "Term-AGA")
    m  <- glm(y ~ svn, data = df, family = binomial())
    est <- broom::tidy(m) %>% filter(term != "(Intercept)")
    ci  <- suppressMessages(confint.default(m))[-1, , drop = FALSE]
    tibble(sb_scope = scope_label,
           category = sub("^svn", "", est$term),
           OR = exp(est$estimate), CI_lower = exp(ci[, 1]), CI_upper = exp(ci[, 2]),
           p_value = est$p.value)
}
all_or  <- or_for(base,                                            "sb",      "All stillbirths")
# late: drop early stillbirths (sb & <28w); outcome = late stillbirth vs live birth
late_or <- or_for(base %>% filter(!(sb == 1 & gestation_weeks < 28)), "late_sb", "Late stillbirth (>=28w)")

ref <- crossing(sb_scope = c("All stillbirths", "Late stillbirth (>=28w)"), category = "Term-AGA") %>%
    mutate(OR = 1, CI_lower = 1, CI_upper = 1, p_value = NA_real_)

or_plot <- bind_rows(all_or, late_or, ref) %>%
    mutate(
        svn_cat = factor(category, levels = PH,
                         labels = c("Term\nAGA", "Term\nSGA", "Preterm\nAGA", "Preterm\nSGA")),
        sig = ifelse(!is.na(p_value) & p_value < 0.05, "*", ""),
        p_lbl = ifelse(is.na(p_value), "",
                       ifelse(p_value < 0.001, "p<0.001", sprintf("p=%.3f", p_value))),
        lbl = ifelse(category == "Term-AGA", "1", sprintf("%.1f%s\n%s", OR, sig, p_lbl))
    )
write.csv(or_plot %>% select(sb_scope, category, OR, CI_lower, CI_upper, p_value),
          file.path(OUTPUT_DIR, "figure_02_stillbirth_or_data.csv"), row.names = FALSE)

# ---- plot --------------------------------------------------------------------
COLS <- c("Term\nAGA" = "#9E9E9E", "Term\nSGA" = "#4A90D9",
          "Preterm\nAGA" = "#F2C744", "Preterm\nSGA" = "#D98CB3")
mk <- function(scope, ttl) {
    d <- or_plot %>% filter(sb_scope == scope)
    ggplot(d, aes(svn_cat, OR, colour = svn_cat)) +
        geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
        geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.15, linewidth = 1) +
        geom_point(size = 4) +
        geom_text(aes(y = CI_upper, label = lbl), vjust = -0.9, size = 3.9,
                  colour = "grey20", lineheight = 0.9) +
        scale_colour_manual(values = COLS, guide = "none") +
        scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.28))) +
        scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
        labs(title = ttl, x = "Newborn type", y = "Stillbirth odds ratio") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold"),
              axis.title = element_text(size = 13.5),
              axis.text  = element_text(size = 12),
              panel.grid.major.x = element_blank())
}
fig <- mk("All stillbirths", "A  All stillbirths") +
       mk("Late stillbirth (>=28w)", "B  Late stillbirths (>=28 weeks)") +
       plot_annotation(
           title = "Stillbirth odds ratio by SVN phenotype (combined sites, vs Term-AGA)",
           subtitle = "Reference: Term-AGA. Error bars = 95% CI; point = OR. Panel A values match Table 1.")
ggsave(file.path(OUTPUT_DIR_GRAPHS, "figure_02_stillbirth_or_combined.png"), fig, width = 13, height = 6.5, dpi = 300)
ggsave(file.path(OUTPUT_DIR_GRAPHS, "figure_02_stillbirth_or_combined.pdf"), fig, width = 13, height = 6.5)
cat("\nSaved figure_02_stillbirth_or_combined.png / .pdf and figure_02_stillbirth_or_data.csv\n")
print(as.data.frame(or_plot %>% filter(category != "Term-AGA") %>%
      select(sb_scope, category, OR, CI_lower, CI_upper)))
cat("Figure_02_stillbirth_OR_by_SVN.R complete.\n")
