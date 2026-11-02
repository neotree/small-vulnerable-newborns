# ==============================================================================
# Purpose: Table 11 -- case fatality with vs without antenatal corticosteroids
# among eligible preterm neonates, restricted to the clinically eligible
# window (24+0 to 34+6 weeks gestation).
#
# Definitions: the `ansteroids` variable is grouped into Steroids Used (Y) /
# Not Used (N) / Unsure (U) / Missing. Reports CFR (deaths/n) for "Steroids
# used" vs "not used", by facility and by LBW/SGA category, with Fisher's
# exact test p-values.
#
# Output: table_11_full_*.csv, table_11_ga24to34wk_*.csv, table_11_cfr_by_steroids.csv
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

# Steroids Groups
STEROIDS_GROUPS <- c(
    "Steroids Used (Y)",
    "Steroids Not Used (N)",
    "Unsure (U)",
    "Missing"
)

# Normalize outcomes
normalise_outcome <- function(x) {
    y <- tolower(trimws(as.character(x)))
    y <- gsub("[^a-z0-9]+", " ", y)
    y <- gsub("\\s+", " ", y)

    dplyr::case_when(
        grepl("\\bdc\\b|discharged|\\bdcr\\b|\\bdpc\\b|\\bdckmc\\b|kangaroo mother care", y) ~ "Alive",
        grepl("\\babs\\b|abscond|\\bdama\\b|against medical advice", y) ~ "Alive",
        grepl("\\btrh\\b|\\btro\\b|transfer", y) ~ "Alive",
        grepl("\\bnnd\\b|neonatal death|\\bdda\\b|died during admission", y) ~ "Neonatal Death",
        grepl("\\bbid\\b|brought in dead", y) ~ "Neonatal Death",
        grepl("\\bstb\\b|stillbirth|still born", y) ~ "Stillbirth",
        y == "" ~ NA_character_,
        TRUE ~ "Unknown/Other"
    )
}

# Normalize Antenatal Steroids
# Values: Y, N, U / Not Sure, or Missing
normalise_steroids <- function(x) {
    y <- toupper(trimws(as.character(x)))
    dplyr::case_when(
        y == "Y" ~ "Steroids Used (Y)",
        y == "N" ~ "Steroids Not Used (N)",
        y == "U" ~ "Unsure (U)",
        # "Not Sure" is a second raw coding of the same answer (219 records in
        # the 24-34+6 wk window); previously it fell through to "Missing".
        y %in% c("NOT SURE", "NOTSURE") ~ "Unsure (U)",
        is.na(y) | y == "" | y == "NA" ~ "Missing",
        TRUE ~ "Missing" # Treat unexpected values as missing or check them
    )
}

# ===== LOAD AND PROCESS DATA =====

cat("\n=== Loading NNU Data ===\n")
all_admissions <- load_nnu_data()
recorded_sex <- if ("gender_std" %in% names(all_admissions)) {
    all_admissions$gender_std
} else if ("gender" %in% names(all_admissions)) {
    all_admissions$gender
} else {
    rep(NA_character_, nrow(all_admissions))
}

# Extract antenatal steroids column
all_admissions <- all_admissions %>%
    mutate(steroids_raw = ansteroids)

# Process
cat("Processing variables...\n")
analysis_df <- all_admissions %>%
    mutate(
        bw_kg = bw_kg_canonical(.),
        gestation_weeks = suppressWarnings(as.numeric(gestation)),
        gender_std = recorded_sex,

        # Outcome
        outcome_norm = normalise_outcome(neotreeoutcome),

        # Steroids Categorization
        steroids_cat = normalise_steroids(steroids_raw),

        # Validations
        is_valid_bw = is.finite(bw_kg) & bw_kg >= BW_MIN_KG & bw_kg <= BW_MAX_KG,
        is_invalid_ga = is.na(gestation_weeks) | gestation_weeks < GA_MIN | gestation_weeks > GA_MAX,
        sga_status = dplyr::if_else(gestation_weeks >= GA_MIN & !is.na(gestation_weeks), classify_sga(bw_kg, gestation_weeks, gender_std), NA_character_),
        sga_cat_temp = sga_category_from_gestational_age(sga_status, gestation_weeks),
        is_valid_sga = is_valid_bw & !is.na(sga_cat_temp),

        # Categorizations including "Discarded"
        bw_category = if_else(is_valid_bw, bw_category_from_kg(bw_kg), "Discarded"),
        lbw_explanation = if_else(is_valid_bw, "Included in analysis", "Missing or invalid birth weight"),
        sga_explanation = case_when(
            is_valid_sga ~ "Included in analysis",
            !is_valid_bw & is_invalid_ga ~ "Missing both BW and GA",
            !is_valid_bw ~ "Missing BW",
            is_invalid_ga ~ "Missing GA",
            TRUE ~ "Missing GA"
        ),
        sga_category = if_else(is_valid_sga, sga_cat_temp, sga_explanation)
    )

# Filter for valid data (Valid Outcome) - Birthweight mappings managed securely within explicit matrices natively gracefully natively accurately explicitly
valid_df <- analysis_df %>%
    filter(outcome_norm %in% c("Alive", "Neonatal Death"))

cat("Total Analysis Population (Valid Outcome):", nrow(valid_df), "\n")

# ===== ANALYSIS FUNCTION =====

run_q21_analysis <- function(valid_df, analysis_df, file_prefix, label) {
    cat(sprintf("\n=== Running Analysis: %s ===\n", label))
    cat("Total Analysis Population (Valid Outcome):", nrow(valid_df), "\n")

    # ----- CALC STATISTICS -----
    cat("Calculating Statistics...\n")

    # Helper for group stats
    calc_group_stats <- function(data, group_name, filter_expr, group_cols, group_levels) {
        category_col <- group_cols[1]
        stats <- data %>%
            filter(eval(parse(text = filter_expr))) %>%
            filter(!is.na(.data[[category_col]])) %>%
            group_by(facility, category = .data[[category_col]]) %>%
            summarise(
                NND = sum(outcome_norm == "Neonatal Death"),
                DC = sum(outcome_norm == "Alive"),
                .groups = "drop"
            ) %>%
            mutate(steroids_group = group_name)

        totals <- data %>%
            filter(eval(parse(text = filter_expr))) %>%
            group_by(facility) %>%
            summarise(
                NND = sum(outcome_norm == "Neonatal Death"),
                DC = sum(outcome_norm == "Alive"),
                .groups = "drop"
            ) %>%
            mutate(category = "Total", steroids_group = group_name)
        bind_rows(stats, totals)
    }

    # Helper for wide outcome
    calc_wide_outcome <- function(data, group_cols, group_levels) {
        defs <- list(
            list(name = "Steroids Used (Y)", expr = "steroids_cat == 'Steroids Used (Y)'"),
            list(name = "Steroids Not Used (N)", expr = "steroids_cat == 'Steroids Not Used (N)'"),
            list(name = "Unsure (U)", expr = "steroids_cat == 'Unsure (U)'"),
            list(name = "Missing", expr = "steroids_cat == 'Missing'")
        )
        all_stats <- list()
        for (d in defs) {
            all_stats[[d$name]] <- calc_group_stats(data, d$name, d$expr, group_cols, group_levels)
        }
        summary_long <- bind_rows(all_stats)
        complete_grid <- expand_grid(
            facility = unique(data$facility),
            category = c(group_levels, "Discarded", "Missing BW", "Missing GA", "Missing both BW and GA", "Total"),
            steroids_group = STEROIDS_GROUPS
        )
        summary_complete <- complete_grid %>%
            left_join(summary_long, by = c("facility", "category", "steroids_group")) %>%
            replace_na(list(NND = 0, DC = 0))

        summary_wide <- summary_complete %>%
            pivot_wider(
                names_from = category,
                values_from = c(NND, DC),
                names_glue = "{category}_{.value}"
            )

        col_order <- c("facility", "steroids_group")
        for (cat in c(group_levels, "Discarded", "Missing BW", "Missing GA", "Missing both BW and GA", "Total")) {
            col_name_nnd <- paste0(cat, "_NND")
            col_name_dc <- paste0(cat, "_DC")
            if (col_name_nnd %in% names(summary_wide)) {
                col_order <- c(col_order, col_name_nnd, col_name_dc)
            }
        }
        summary_wide %>% select(any_of(col_order)) %>% arrange(facility, match(steroids_group, STEROIDS_GROUPS))
    }

    # LBW Categories
    lbw_levels <- c("ELBW", "VLBW", "LBW", "NBW", "HBW")
    lbw_wide <- calc_wide_outcome(valid_df, c("bw_category", "lbw_explanation"), lbw_levels) %>%
        mutate(classification = "Birthweight")

    # SGA Categories
    sga_levels <- c("Term-AGA", "Term-SGA", "Preterm-AGA", "Preterm-SGA")
    sga_wide <- calc_wide_outcome(valid_df, c("sga_category", "sga_explanation"), sga_levels) %>%
        mutate(classification = "SGA")

    # ----- SAVE RESULTS -----
    cat("Saving Results...\n")

    bw_note <- sprintf("is.finite(bw_kg) & bw_kg >= %.1f kg & bw_kg <= %.1f kg", BW_MIN_KG, BW_MAX_KG)
    ga_note <- sprintf("!is.na(gestation) & gestation >= %d wks & gestation <= %d wks", GA_MIN, GA_MAX)
    sga_note <- paste0(bw_note, " AND ", ga_note, " AND !is.na(sga_category)")
    outcome_note <- "normalise_outcome(neotreeoutcome) %in% c(\"Alive\", \"Neonatal Death\")"
    
    n_total <- nrow(analysis_df)
    sample_size_summary <- tibble(
        description = c(
            "Total admissions",
            "Analysis population (valid outcome)",
            "Excluded: missing/invalid outcome",
            "Analysis population - Steroids Used (Y)",
            "Analysis population - Steroids Not Used (N)",
            "Analysis population - Unsure (U)",
            "Analysis population - Missing steroids record",
            "Analysis population (KCH) - Steroids Used (Y)",
            "Analysis population (KCH) - Steroids Not Used (N)",
            "Analysis population (KCH) - Unsure (U)",
            "Analysis population (KCH) - Missing steroids record",
            "Analysis population (SMCH) - Steroids Used (Y)",
            "Analysis population (SMCH) - Steroids Not Used (N)",
            "Analysis population (SMCH) - Unsure (U)",
            "Analysis population (SMCH) - Missing steroids record",
            "Analysis population with valid birthweight (LBW sample)",
            "Discarded for LBW (missing/invalid birthweight)",
            "Analysis population with valid SGA classification",
            "Discarded for SGA (Missing BW)",
            "Discarded for SGA (Missing GA)",
            "Discarded for SGA (Missing both BW and GA)"
        ),
        n = c(
            n_total,
            nrow(valid_df),
            sum(!analysis_df$outcome_norm %in% c("Alive", "Neonatal Death")),
            sum(valid_df$steroids_cat == "Steroids Used (Y)"),
            sum(valid_df$steroids_cat == "Steroids Not Used (N)"),
            sum(valid_df$steroids_cat == "Unsure (U)"),
            sum(valid_df$steroids_cat == "Missing"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Steroids Used (Y)"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Steroids Not Used (N)"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Unsure (U)"),
            sum(valid_df$facility == "KCH" & valid_df$steroids_cat == "Missing"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Steroids Used (Y)"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Steroids Not Used (N)"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Unsure (U)"),
            sum(valid_df$facility == "SMCH" & valid_df$steroids_cat == "Missing"),
            sum(valid_df$is_valid_bw),
            sum(!valid_df$is_valid_bw),
            sum(valid_df$is_valid_sga),
            sum(valid_df$sga_explanation == "Missing BW"),
            sum(valid_df$sga_explanation == "Missing GA"),
            sum(valid_df$sga_explanation == "Missing both BW and GA")
        ),
        pct_of_total = round(100 * n / n_total, 1),
        Filter_notes = c(
            NA_character_, outcome_note, paste0("NOT (", outcome_note, ")"),
            paste0(outcome_note, " AND steroids == 'Y'"), paste0(outcome_note, " AND steroids == 'N'"),
            paste0(outcome_note, " AND steroids == 'U'"), paste0(outcome_note, " AND steroids is NA/blank/other"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids == 'Y'"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids == 'N'"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids == 'U'"),
            paste0(outcome_note, " AND facility == 'KCH' AND steroids is Missing"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids == 'Y'"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids == 'N'"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids == 'U'"),
            paste0(outcome_note, " AND facility == 'SMCH' AND steroids is Missing"),
            paste0(outcome_note, " AND ", bw_note),
            paste0(outcome_note, " AND NOT (", bw_note, ")"),
            paste0(outcome_note, " AND ", sga_note),
            "Invalid/missing birth weight", "Invalid/missing gestational age", "Both BW and GA invalid/missing"
        )
    )

    write_csv(sample_size_summary, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_00_sample_size_summary.csv")))
    write_csv(lbw_wide, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_01_steroids_outcomes_LBW_wide.csv")))
    write_csv(sga_wide, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_02_steroids_outcomes_SGA_wide.csv")))

    # ----- CFR TESTS -----
    cat("Running CFR Chi-Square and Pairwise Tests...\n")
    
    calc_cfr_tests <- function(data, category_col, valid_categories, group_col, group_levels, group_label) {
        d <- data %>% filter(.data[[group_col]] %in% group_levels, .data[[category_col]] %in% valid_categories)
        facilities <- unique(d$facility)
        results <- list()
        for (fac in facilities) {
            for (cat in valid_categories) {
                sub <- d %>% filter(facility == fac, .data[[category_col]] == cat)
                ct_df <- sub %>% group_by(group = .data[[group_col]]) %>%
                    summarise(NND = sum(outcome_norm == "Neonatal Death"), DC  = sum(outcome_norm == "Alive"), .groups = "drop") %>%
                    filter((NND + DC) > 0)
                if (nrow(ct_df) < 2) next
                ct_mat <- as.matrix(ct_df[, c("NND", "DC")]); rownames(ct_mat) <- ct_df$group
                cfr_vals <- ct_df %>% mutate(total = NND + DC, CFR = round(100 * NND / total, 1))
                chisq_res <- tryCatch(suppressWarnings(chisq.test(ct_mat, correct = FALSE)), error = function(e) NULL)
                overall_p <- if (!is.null(chisq_res)) chisq_res$p.value else NA_real_
                overall_stat <- if (!is.null(chisq_res)) chisq_res$statistic else NA_real_
                overall_df   <- if (!is.null(chisq_res)) chisq_res$parameter else NA_real_
                expected_warn <- if (!is.null(chisq_res)) any(chisq_res$expected < 5) else NA
                grp_names <- rownames(ct_mat); pairs <- combn(grp_names, 2, simplify = FALSE); n_pairs <- length(pairs)
                pairwise_rows <- map_dfr(pairs, function(pair) {
                    g1 <- pair[1]; g2 <- pair[2]; sub_mat <- ct_mat[c(g1, g2), , drop = FALSE]
                    fish_p <- tryCatch(fisher.test(sub_mat)$p.value, error = function(e) NA_real_)
                    tibble(facility = fac, category = cat, comparison = paste0(g1, " vs. ", g2), group1 = g1, group2 = g2,
                           CFR_group1_pct = cfr_vals$CFR[cfr_vals$group == g1], n_group1 = cfr_vals$total[cfr_vals$group == g1], NND_group1 = cfr_vals$NND[cfr_vals$group == g1],
                           CFR_group2_pct = cfr_vals$CFR[cfr_vals$group == g2], n_group2 = cfr_vals$total[cfr_vals$group == g2], NND_group2 = cfr_vals$NND[cfr_vals$group == g2],
                           fisher_p_raw = fish_p, n_pairs_tested = n_pairs)
                })
                pairwise_rows <- pairwise_rows %>% mutate(
                    fisher_p_bonferroni = pmin(fisher_p_raw * n_pairs_tested, 1),
                    fisher_p_bonferroni = round(fisher_p_bonferroni, 6),
                    fisher_p_raw = round(fisher_p_raw, 6),
                    chisq_statistic = round(as.numeric(overall_stat), 4),
                    chisq_df = as.integer(overall_df),
                    chisq_p_value = round(overall_p, 6),
                    low_expected_count_warning = expected_warn,
                    exposure_variable = group_label,
                    classification = if (category_col == "bw_category") "Birthweight" else "SGA"
                )
                results[[paste(fac, cat, sep = "_")]] <- pairwise_rows
            }
        }
        bind_rows(results)
    }

    testable_steroids <- c("Steroids Used (Y)", "Steroids Not Used (N)", "Unsure (U)")
    cfr_tests_lbw <- calc_cfr_tests(valid_df, "bw_category", lbw_levels, "steroids_cat", testable_steroids, "Antenatal Steroids")
    cfr_tests_sga <- calc_cfr_tests(valid_df, "sga_category", sga_levels, "steroids_cat", testable_steroids, "Antenatal Steroids")
    cfr_tests_all <- bind_rows(cfr_tests_lbw, cfr_tests_sga) %>% arrange(facility, classification, category, comparison)

    write_csv(cfr_tests_all, file.path(OUTPUT_DIR_NNU, paste0(file_prefix, "_03_steroids_CFR_chisquare_pairwise.csv")))
    cat(sprintf("Saved: 03_steroids_CFR_chisquare_pairwise.csv (Rows: %d)\n", nrow(cfr_tests_all)))
    invisible(cfr_tests_all)
}

# ===== EXECUTE ANALYSES =====

# 1. Full Sample Analysis
run_q21_analysis(valid_df, analysis_df, "table_11_full", "Full Sample")

# 2. Gestation Subsample Analysis (24 to 34 weeks)
cat("\n=== Creating Gestation Subsample (24 to 34 weeks) ===\n")
sub_analysis_df <- analysis_df %>%
    filter(gestation_weeks >= 24 & gestation_weeks < 35)

sub_valid_df <- sub_analysis_df %>%
    filter(outcome_norm %in% c("Alive", "Neonatal Death"))

res_sub <- run_q21_analysis(sub_valid_df, sub_analysis_df, "table_11_ga24to34wk",
                            "Gestation Subsample (24 to 34 weeks)")

# ===== MANUSCRIPT TABLE 11 =====
# CFR expressed as deaths per 1,000 admissions, with deaths/n shown, No vs Yes
# steroids, for the preterm SVN phenotypes in the eligible GA 24-34+6 wk window.
#   p (No vs Yes) = Fisher exact on the two shown groups (matches a 2-column table)
#   p (3-group)   = 3-group (N/Y/U) chi-square, reported for reference
fmtp    <- function(p) ifelse(is.na(p), "-", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
comma   <- function(x) formatC(x, format = "d", big.mark = ",")
cfr1000 <- function(nnd, n) round(1000 * nnd / n)
cell    <- function(nnd, n) sprintf("%d (%s/%s)", cfr1000(nnd, n), comma(nnd), comma(n))

table13 <- res_sub %>%
    filter(classification == "SGA",
           comparison == "Steroids Not Used (N) vs. Steroids Used (Y)",
           category %in% c("Preterm-AGA", "Preterm-SGA")) %>%
    mutate(.fac = factor(facility, levels = c("KCH", "SMCH")),
           .cat = factor(category, levels = c("Preterm-AGA", "Preterm-SGA"))) %>%
    arrange(.fac, .cat) %>%
    transmute(
        Group = paste(facility, category),
        `No steroids, CFR per 1,000 (deaths/n)` = cell(NND_group1, n_group1),
        `Steroids, CFR per 1,000 (deaths/n)`    = cell(NND_group2, n_group2),
        `p (No vs Yes, Fisher)`  = fmtp(fisher_p_raw),
        `p (3-group chi-square)` = fmtp(chisq_p_value))

write_csv(table13, file.path(OUTPUT_DIR_NNU, "table_11_cfr_by_steroids.csv"))
cat("\n=== Manuscript Table 11: CFR (per 1,000) by antenatal steroids, eligible preterm ===\n")
print(as.data.frame(table13), row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(OUTPUT_DIR_NNU, "table_11_cfr_by_steroids.csv")))

# ==============================================================================
# POOLED AND SITE-ADJUSTED STEROID ANALYSIS (eligible window, 24+0 to 34+6 wk)
# ------------------------------------------------------------------------------
# Added in response to reviewer comment 88 (September 2026): the site-specific
# Preterm-SGA comparisons are small, so a both-sites analysis is needed before
# any claim that steroid benefit differs by growth status.
#
# For each preterm phenotype (Preterm-AGA, Preterm-SGA), restricted to
# neonates with steroids recorded as Y or N:
#   * site-specific OR (Wald 95% CI, univariable logistic model)
#   * crude pooled OR (sites simply combined)
#   * site-adjusted OR: Mantel-Haenszel (Robins-Breslow-Greenland CI) and
#     logistic model death ~ steroids + facility
#   * heterogeneity between sites: Breslow-Day (Tarone-corrected) test and the
#     likelihood-ratio test of a steroids x facility interaction
#   * sensitivity: logistic model additionally adjusted for gestational age
#     (weeks+days converted to decimal weeks, natural spline, 3 df)
#   * Mantel-Haenszel risk difference (Greenland-Robins variance) and the
#     corresponding number needed to treat
# Effect modification by growth status is tested in both phenotypes together:
#   death ~ steroids * phenotype + facility * phenotype   (LRT on steroids:phenotype)
#
# These are observational associations; steroid receipt is not randomised and
# is likely to be confounded (antenatal care, referral pathway, time to birth).
# ==============================================================================

cat("\n=== Pooled and site-adjusted steroid analysis (24+0 to 34+6 wk) ===\n")

PRETERM_PHENOTYPES <- c("Preterm-AGA", "Preterm-SGA")

# Gestation is recorded in weeks+days notation (e.g. 33.4 = 33 weeks 4 days);
# convert to decimal weeks for use as a continuous covariate.
ga_weeks_days_to_decimal <- function(x) {
    wk <- floor(x)
    dy <- round(10 * (x - wk))
    wk + pmin(dy, 6) / 7
}

ster_df <- sub_valid_df %>%
    filter(sga_category %in% PRETERM_PHENOTYPES,
           steroids_cat %in% c("Steroids Used (Y)", "Steroids Not Used (N)")) %>%
    mutate(
        death     = as.integer(outcome_norm == "Neonatal Death"),
        steroids  = factor(if_else(steroids_cat == "Steroids Used (Y)", "Yes", "No"),
                           levels = c("No", "Yes")),
        facility  = factor(facility, levels = c("KCH", "SMCH")),
        phenotype = factor(sga_category, levels = PRETERM_PHENOTYPES),
        ga_dec    = ga_weeks_days_to_decimal(gestation_weeks)
    )

cat(sprintf("Steroid analysis population (Y or N, preterm AGA/SGA): %d\n", nrow(ster_df)))

wald_or <- function(fit, term = "steroidsYes") {
    est <- coef(summary(fit))[term, ]
    c(OR = exp(est[["Estimate"]]),
      lo = exp(est[["Estimate"]] - qnorm(0.975) * est[["Std. Error"]]),
      hi = exp(est[["Estimate"]] + qnorm(0.975) * est[["Std. Error"]]),
      p  = est[["Pr(>|z|)"]])
}

# Stratum tables in the layout [steroids (Yes, No)] x [death (1, 0)] x [facility],
# so that mantelhaen.test() returns the OR for steroids Yes vs No.
strata_array <- function(d) {
    table(steroids = factor(d$steroids, levels = c("Yes", "No")),
          death    = factor(d$death, levels = c(1, 0)),
          facility = droplevels(d$facility))
}

# Breslow-Day test for homogeneity of ORs, with Tarone's correction
breslow_day_test <- function(tab, or_mh) {
    K <- dim(tab)[3]
    stat <- 0; tarone_num <- 0; tarone_var <- 0
    for (k in seq_len(K)) {
        a  <- tab[1, 1, k]                  # exposed deaths
        m1 <- sum(tab[, 1, k])              # total deaths
        n1 <- sum(tab[1, , k])              # exposed total
        n0 <- sum(tab[2, , k])              # unexposed total
        lo <- max(0, m1 - n0); hi <- min(n1, m1)
        # expected exposed deaths under common OR (quadratic root in range)
        f <- function(x) x * (n0 - m1 + x) - or_mh * (n1 - x) * (m1 - x)
        ea <- uniroot(f, c(lo, hi), tol = 1e-10)$root
        va <- 1 / (1 / ea + 1 / (n1 - ea) + 1 / (m1 - ea) + 1 / (n0 - m1 + ea))
        stat <- stat + (a - ea)^2 / va
        tarone_num <- tarone_num + (a - ea)
        tarone_var <- tarone_var + va
    }
    stat_t <- stat - tarone_num^2 / tarone_var
    c(chisq = stat_t, df = K - 1, p = pchisq(stat_t, df = K - 1, lower.tail = FALSE))
}

# Mantel-Haenszel risk difference with Greenland-Robins variance
mh_risk_difference <- function(tab) {
    num <- 0; den <- 0; var_num <- 0
    for (k in seq_len(dim(tab)[3])) {
        a <- tab[1, 1, k]; b <- tab[1, 2, k]   # exposed: deaths, survivors
        c <- tab[2, 1, k]; d <- tab[2, 2, k]   # unexposed: deaths, survivors
        n1 <- a + b; n0 <- c + d; N <- n1 + n0
        if (n1 == 0 || n0 == 0) next
        num <- num + (a * n0 - c * n1) / N
        den <- den + n1 * n0 / N
        var_num <- var_num + (a * b * n0^3 + c * d * n1^3) / (n1 * n0 * N^2)
    }
    rd <- num / den; se <- sqrt(var_num) / den
    c(RD = rd, lo = rd - qnorm(0.975) * se, hi = rd + qnorm(0.975) * se, se = se)
}

# NNT from a risk difference (steroids minus no steroids; negative = benefit)
nnt_label <- function(rd, lo, hi) {
    if (rd >= 0) return("Not applicable (no benefit)")
    if (hi < 0) return(sprintf("%.0f (%.0f to %.0f)", -1 / rd, -1 / lo, -1 / hi))
    sprintf("%.0f (NNT %.0f to benefit, NNH %.0f to harm)", -1 / rd, -1 / lo, 1 / hi)
}

cell_counts <- function(d) {
    no  <- d %>% filter(steroids == "No");  yes <- d %>% filter(steroids == "Yes")
    tibble(deaths_no = sum(no$death),  n_no  = nrow(no),
           deaths_yes = sum(yes$death), n_yes = nrow(yes))
}

pooled_rows <- list()
het_rows    <- list()

for (ph in PRETERM_PHENOTYPES) {
    d_ph <- ster_df %>% filter(phenotype == ph)

    # --- site-specific ---
    for (fac in levels(ster_df$facility)) {
        d_s <- d_ph %>% filter(facility == fac)
        fit <- glm(death ~ steroids, family = binomial, data = d_s)
        o   <- wald_or(fit)
        rd  <- mh_risk_difference(strata_array(d_s))
        pooled_rows[[paste(ph, fac)]] <- bind_cols(
            tibble(phenotype = ph, analysis = paste0(fac, " only"), model = "Univariable logistic"),
            cell_counts(d_s),
            tibble(OR = o[["OR"]], OR_lo = o[["lo"]], OR_hi = o[["hi"]], p = o[["p"]],
                   fisher_p = fisher.test(table(d_s$steroids, d_s$death))$p.value,
                   RD_per_1000 = 1000 * rd[["RD"]], RD_lo = 1000 * rd[["lo"]], RD_hi = 1000 * rd[["hi"]],
                   NNT = nnt_label(rd[["RD"]], rd[["lo"]], rd[["hi"]])))
    }

    # --- crude pool ---
    fit_crude <- glm(death ~ steroids, family = binomial, data = d_ph)
    o_crude   <- wald_or(fit_crude)
    tab_crude <- table(steroids = factor(d_ph$steroids, levels = c("Yes", "No")),
                       death = factor(d_ph$death, levels = c(1, 0)))
    rd_crude  <- mh_risk_difference(array(tab_crude, dim = c(2, 2, 1)))
    pooled_rows[[paste(ph, "crude")]] <- bind_cols(
        tibble(phenotype = ph, analysis = "Both sites, crude pool", model = "Univariable logistic"),
        cell_counts(d_ph),
        tibble(OR = o_crude[["OR"]], OR_lo = o_crude[["lo"]], OR_hi = o_crude[["hi"]], p = o_crude[["p"]],
               fisher_p = fisher.test(tab_crude)$p.value,
               RD_per_1000 = 1000 * rd_crude[["RD"]], RD_lo = 1000 * rd_crude[["lo"]], RD_hi = 1000 * rd_crude[["hi"]],
               NNT = nnt_label(rd_crude[["RD"]], rd_crude[["lo"]], rd_crude[["hi"]])))

    # --- site-adjusted: Mantel-Haenszel ---
    tab <- strata_array(d_ph)
    mh  <- mantelhaen.test(tab, correct = FALSE)
    rd  <- mh_risk_difference(tab)
    pooled_rows[[paste(ph, "MH")]] <- bind_cols(
        tibble(phenotype = ph, analysis = "Both sites, site-stratified", model = "Mantel-Haenszel"),
        cell_counts(d_ph),
        tibble(OR = unname(mh$estimate), OR_lo = mh$conf.int[1], OR_hi = mh$conf.int[2],
               p = mh$p.value, fisher_p = NA_real_,
               RD_per_1000 = 1000 * rd[["RD"]], RD_lo = 1000 * rd[["lo"]], RD_hi = 1000 * rd[["hi"]],
               NNT = nnt_label(rd[["RD"]], rd[["lo"]], rd[["hi"]])))

    # --- site-adjusted: logistic ---
    fit_adj <- glm(death ~ steroids + facility, family = binomial, data = d_ph)
    o_adj   <- wald_or(fit_adj)
    stopifnot(sign(log(unname(mh$estimate))) == sign(log(o_adj[["OR"]])))
    pooled_rows[[paste(ph, "adj")]] <- bind_cols(
        tibble(phenotype = ph, analysis = "Both sites, site-adjusted", model = "Logistic: steroids + facility"),
        cell_counts(d_ph),
        tibble(OR = o_adj[["OR"]], OR_lo = o_adj[["lo"]], OR_hi = o_adj[["hi"]], p = o_adj[["p"]],
               fisher_p = NA_real_, RD_per_1000 = NA_real_, RD_lo = NA_real_, RD_hi = NA_real_, NNT = NA_character_))

    # --- sensitivity: + gestational age ---
    d_ga <- d_ph %>% filter(is.finite(ga_dec))
    fit_ga <- glm(death ~ steroids + facility + splines::ns(ga_dec, 3), family = binomial, data = d_ga)
    o_ga   <- wald_or(fit_ga)
    pooled_rows[[paste(ph, "ga")]] <- bind_cols(
        tibble(phenotype = ph, analysis = "Both sites, site- and GA-adjusted (sensitivity)",
               model = "Logistic: steroids + facility + ns(GA, 3)"),
        cell_counts(d_ga),
        tibble(OR = o_ga[["OR"]], OR_lo = o_ga[["lo"]], OR_hi = o_ga[["hi"]], p = o_ga[["p"]],
               fisher_p = NA_real_, RD_per_1000 = NA_real_, RD_lo = NA_real_, RD_hi = NA_real_, NNT = NA_character_))

    # --- heterogeneity between sites ---
    bd      <- breslow_day_test(tab, unname(mh$estimate))
    fit_int <- glm(death ~ steroids * facility, family = binomial, data = d_ph)
    lrt     <- anova(fit_adj, fit_int, test = "LRT")
    or_k    <- sapply(levels(ster_df$facility), function(f) {
        wald_or(glm(death ~ steroids, family = binomial, data = d_ph %>% filter(facility == f)))[["OR"]]
    })
    het_rows[[ph]] <- tibble(
        phenotype = ph,
        test = c("Breslow-Day (Tarone) homogeneity of OR across sites",
                 "LRT steroids x facility interaction"),
        statistic = c(bd[["chisq"]], lrt$Deviance[2]),
        df = c(bd[["df"]], lrt$Df[2]),
        p = c(bd[["p"]], lrt$`Pr(>Chi)`[2]),
        note = sprintf("Site ORs: KCH %.2f, SMCH %.2f; ratio SMCH/KCH %.2f",
                       or_k[["KCH"]], or_k[["SMCH"]], or_k[["SMCH"]] / or_k[["KCH"]])
    )
}

# --- effect modification by growth status (both phenotypes together) ---
fit_main <- glm(death ~ steroids + phenotype * facility, family = binomial, data = ster_df)
fit_emod <- glm(death ~ steroids * phenotype + phenotype * facility, family = binomial, data = ster_df)
lrt_emod <- anova(fit_main, fit_emod, test = "LRT")
rom      <- coef(summary(fit_emod))["steroidsYes:phenotypePreterm-SGA", ]

# difference in MH risk differences (SGA minus AGA)
rd_by_ph <- lapply(PRETERM_PHENOTYPES, function(ph) mh_risk_difference(strata_array(ster_df %>% filter(phenotype == ph))))
names(rd_by_ph) <- PRETERM_PHENOTYPES
rd_diff    <- rd_by_ph[["Preterm-SGA"]][["RD"]] - rd_by_ph[["Preterm-AGA"]][["RD"]]
rd_diff_se <- sqrt(rd_by_ph[["Preterm-SGA"]][["se"]]^2 + rd_by_ph[["Preterm-AGA"]][["se"]]^2)

het_rows[["growth"]] <- tibble(
    phenotype = "Preterm-AGA vs Preterm-SGA",
    test = c("LRT steroids x phenotype interaction (site-adjusted)",
             "Difference in MH risk difference, SGA minus AGA (per 1,000)"),
    statistic = c(lrt_emod$Deviance[2], 1000 * rd_diff),
    df = c(lrt_emod$Df[2], NA_real_),
    p = c(lrt_emod$`Pr(>Chi)`[2], 2 * pnorm(-abs(rd_diff / rd_diff_se))),
    note = c(sprintf("Ratio of ORs (SGA/AGA) %.2f (95%% CI %.2f to %.2f)",
                     exp(rom[["Estimate"]]),
                     exp(rom[["Estimate"]] - qnorm(0.975) * rom[["Std. Error"]]),
                     exp(rom[["Estimate"]] + qnorm(0.975) * rom[["Std. Error"]])),
             sprintf("95%% CI %.0f to %.0f per 1,000",
                     1000 * (rd_diff - qnorm(0.975) * rd_diff_se),
                     1000 * (rd_diff + qnorm(0.975) * rd_diff_se)))
)

pooled_tbl <- bind_rows(pooled_rows) %>%
    mutate(CFR_no_per_1000 = round(1000 * deaths_no / n_no),
           CFR_yes_per_1000 = round(1000 * deaths_yes / n_yes),
           .after = n_yes) %>%
    mutate(across(c(OR, OR_lo, OR_hi), ~ round(.x, 3)),
           across(c(p, fisher_p), ~ signif(.x, 4)),
           across(c(RD_per_1000, RD_lo, RD_hi), ~ round(.x, 1)))
het_tbl <- bind_rows(het_rows) %>%
    mutate(statistic = round(statistic, 3), p = signif(p, 4))

write_csv(pooled_tbl, file.path(OUTPUT_DIR_NNU, "table_11_ga24to34wk_04_steroids_pooled_site_adjusted.csv"))
write_csv(het_tbl,    file.path(OUTPUT_DIR_NNU, "table_11_ga24to34wk_05_steroids_heterogeneity.csv"))

# --- completeness of the steroid field, by raw coding ---
STEROID_RECORD_LEVELS <- c("Yes (Y)", "No (N)", "Unsure (U)", "Not sure", "Missing (blank)")
label_steroid_record <- function(x) {
    y <- toupper(trimws(as.character(x)))
    dplyr::case_when(
        is.na(y) | y == "" ~ "Missing (blank)",
        y == "Y" ~ "Yes (Y)",
        y == "N" ~ "No (N)",
        y == "U" ~ "Unsure (U)",
        y == "NOT SURE" ~ "Not sure",
        TRUE ~ paste0("Other: ", x))
}

rec_df <- sub_valid_df %>%
    mutate(steroid_record = factor(label_steroid_record(steroids_raw),
                                   levels = STEROID_RECORD_LEVELS))
stopifnot(!anyNA(rec_df$steroid_record))   # fails loudly if a new raw coding appears

completeness <- bind_rows(
        rec_df %>% count(facility, steroid_record, .drop = FALSE),
        rec_df %>% count(steroid_record, .drop = FALSE) %>% mutate(facility = "Both")) %>%
    group_by(facility) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup() %>%
    arrange(factor(facility, levels = c("KCH", "SMCH", "Both")), steroid_record)

write_csv(completeness, file.path(OUTPUT_DIR_NNU, "table_11_ga24to34wk_06_steroid_completeness.csv"))

# --- Supplementary Table S7 (manuscript layout, with pooled rows) ---
fmt_or <- function(or, lo, hi) sprintf("%.2f (%.2f-%.2f)", or, lo, hi)
s7 <- pooled_tbl %>%
    filter(analysis %in% c("KCH only", "SMCH only", "Both sites, site-stratified")) %>%
    mutate(Group = case_when(
               analysis == "KCH only"  ~ paste("KCH", phenotype),
               analysis == "SMCH only" ~ paste("SMCH", phenotype),
               TRUE ~ paste("Both sites,", phenotype, "(site-stratified)")),
           .ord = match(analysis, c("KCH only", "SMCH only", "Both sites, site-stratified"))) %>%
    arrange(phenotype, .ord) %>%
    transmute(
        Group,
        `No steroids, CFR per 1,000 (deaths/n)` = cell(deaths_no, n_no),
        `Steroids, CFR per 1,000 (deaths/n)`    = cell(deaths_yes, n_yes),
        `OR (95% CI)` = fmt_or(OR, OR_lo, OR_hi),
        `p` = fmtp(if_else(is.na(fisher_p), p, fisher_p)),
        `Test` = if_else(is.na(fisher_p), "Mantel-Haenszel", "Fisher exact"))

write_csv(s7, file.path(OUTPUT_DIR_NNU, "supplementary_table_s7_steroids_with_pooled.csv"))

cat("\n--- Pooled / site-adjusted steroid ORs ---\n")
print(as.data.frame(pooled_tbl %>% select(phenotype, analysis, n_no, CFR_no_per_1000, n_yes,
                                          CFR_yes_per_1000, OR, OR_lo, OR_hi, p, RD_per_1000, NNT)),
      row.names = FALSE)
cat("\n--- Heterogeneity ---\n")
print(as.data.frame(het_tbl), row.names = FALSE)
cat("\n--- Steroid field completeness ---\n")
print(as.data.frame(completeness), row.names = FALSE)
cat("\n--- Supplementary Table S7 ---\n")
print(as.data.frame(s7), row.names = FALSE)

cat("\nAll analyses complete!\n")
