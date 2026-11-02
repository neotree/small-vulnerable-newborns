# ==============================================================================
# Purpose: two maternal-population outputs, on all births with valid birth
# weight and gestational age:
#   (1) SVN phenotype distribution by site (KCH / SMCH / Combined), n and %
#       of deliveries -- feeds Table 1's "SVN phenotype (all births)" block.
#   (2) Pooled stillbirth odds ratio + 95% CI + p for each category vs its
#       reference (Term-AGA for phenotype; NBW for birth weight).
#
# Definitions: denominator = all births (live births + stillbirths).
#
# Output: helper_table01_svn_by_site_allbirths.csv,
#         helper_table01_stillbirth_or_svn.csv, helper_table01_stillbirth_or_bw.csv
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
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(broom)})
cat("\n=== 55_table1_svn_by_site_and_stillbirth_OR.R ===\n")

STILLBIRTH_CODES <- c("SBF", "SBM", "STBF", "STBM")
GA_MIN_M <- if (exists("GA_MIN_MAT")) GA_MIN_MAT else 22
PH  <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")
BWC <- c("ELBW", "VLBW", "LBW", "NBW", "HBW")

# ---- base maternal population: all births with valid BW + GA -----------------
base <- process_maternal_data(load_maternal_data()) %>%
    filter(is.finite(bw_kg), bw_kg >= BW_MIN_KG, bw_kg <= BW_MAX_KG,
           is.finite(gestation_weeks), gestation_weeks >= GA_MIN_M, gestation_weeks <= GA_MAX) %>%
    mutate(
        facility     = factor(as.character(facility), levels = c("KCH", "SMCH")),
        sga_status   = classify_sga(bw_kg, gestation_weeks, gender_std),
        svn_category = factor(make_svn_category(sga_status, gestation_weeks), levels = PH),
        bw_cat       = factor(bw_category_from_kg(bw_kg), levels = BWC),
        stillbirth   = as.integer(!is.na(neotreeoutcome) &
                                  toupper(trimws(neotreeoutcome)) %in% STILLBIRTH_CODES)
    )
cat(sprintf("All-births base: %d (KCH %d, SMCH %d); stillbirths %d\n",
            nrow(base), sum(base$facility == "KCH"), sum(base$facility == "SMCH"),
            sum(base$stillbirth)))

# ---- (1) SVN phenotype distribution by site (all births) --------------------
svn_site <- function(df, label) {
    tot <- nrow(df)
    df %>% filter(!is.na(svn_category)) %>% count(svn_category, .drop = FALSE) %>%
        mutate(site = label, total = tot, pct = round(100 * n / tot, 1))
}
svn_by_site <- bind_rows(
    svn_site(base %>% filter(facility == "KCH"),  "KCH"),
    svn_site(base %>% filter(facility == "SMCH"), "SMCH"),
    svn_site(base,                                "Combined")
) %>%
    mutate(n_pct = sprintf("%d (%.1f%%)", n, pct)) %>%
    select(site, svn_category, n, total, pct, n_pct) %>%
    pivot_wider(id_cols = svn_category, names_from = site, values_from = n_pct)
write.csv(svn_by_site, file.path(OUT2, "helper_table01_svn_by_site_allbirths.csv"), row.names = FALSE)
cat("\nSVN phenotype by site (all births):\n"); print(as.data.frame(svn_by_site))

# ---- (2) combined stillbirth OR vs reference (pooled sites) ------------------
or_table <- function(df, var, ref) {
    d <- df %>% filter(!is.na(.data[[var]]))
    d[[var]] <- relevel(factor(d[[var]]), ref = ref)
    m <- glm(stillbirth ~ d[[var]], data = d, family = binomial())
    ci <- suppressMessages(confint.default(m))
    broom::tidy(m) %>%
        filter(term != "(Intercept)") %>%
        mutate(category  = sub("^d\\[\\[var\\]\\]", "", term),
               reference = ref,
               OR        = round(exp(estimate), 2),
               CI_low    = round(exp(ci[-1, 1]), 2),
               CI_high   = round(exp(ci[-1, 2]), 2),
               p_value   = signif(p.value, 3),
               p_fmt     = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)),
               OR_CI_p   = sprintf("%.1f (%.1f-%.1f), p %s", OR, CI_low, CI_high, p_fmt)) %>%
        select(category, reference, OR, CI_low, CI_high, p_value, p_fmt, OR_CI_p)
}
or_svn <- or_table(base, "svn_category", "Term-AGA")
or_bw  <- or_table(base, "bw_cat",       "NBW")
write.csv(or_svn, file.path(OUT2, "helper_table01_stillbirth_or_svn.csv"), row.names = FALSE)
write.csv(or_bw,  file.path(OUT2, "helper_table01_stillbirth_or_bw.csv"),  row.names = FALSE)
cat("\nStillbirth OR vs Term-AGA (combined) - Table 2:\n");  print(as.data.frame(or_svn))
cat("\nStillbirth OR vs NBW (combined) - Table 2b:\n");      print(as.data.frame(or_bw))
cat("\nSaved helper_table01_svn_by_site_allbirths.csv, helper_table01_stillbirth_or_svn.csv, helper_table01_stillbirth_or_bw.csv\n")
cat("55_table1_svn_by_site_and_stillbirth_OR.R complete.\n")
