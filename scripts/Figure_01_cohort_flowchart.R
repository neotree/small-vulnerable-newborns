# ==============================================================================
# Purpose: Figure 1 -- CONSORT-style cohort-selection flowchart for the NNU
# (neonatal) dataset.
#
# Definitions -- selection flow (all counts computed from the data):
#   Neonatal unit admissions (KCH + SMCH, Jan 2022 - Dec 2025)
#     -> Excluded: no matched discharge / outcome
#   Matched admission-discharge records
#     -> Excluded: missing/invalid birthweight
#   After birthweight exclusion
#     -> Excluded: missing/invalid gestational age
#   Analysis population (SVN-classified, determinate outcome)
#
# Output: figure_01_cohort_flowchart.png/.pdf (+ _kch, _smch),
#         figure_01_flowchart_counts.csv (+ _kch, _smch)
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
cat("\n=== Figure_01_NNU_cohort_flowchart_v2.R ===\n")

fmt <- function(n) formatC(n, format = "d", big.mark = ",")
pctf <- function(num, den) {
    if (is.na(den) || den == 0) return("0.0%")
    sprintf("%.1f%%", 100 * num / den)
}

# ------------------------------------------------------------------------------
# STEP 0: TOTAL ADMISSIONS (PRE-MATCH DENOMINATOR) FROM MASTER_JOINED FILES
# ------------------------------------------------------------------------------
cat("Reading master_joined files (all admissions)...\n")

read_master_mini <- function(path) {
    suppressMessages(read_csv(
        path,
        col_select = any_of(c("facility", "adm_date_parsed",
                              "match_type", "methodestgest")),
        col_types  = cols(.default = "c"),
        guess_max  = 50000
    ))
}

master <- bind_rows(
    read_master_mini(KCH_NNU_MASTER_FILE),
    read_master_mini(SMCH_NNU_MASTER_FILE)
)

master <- master %>% filter(facility %in% VALID_FACILITIES)
master <- filter_study_period(master, date_col = "adm_date_parsed")

if (POPULATION == "USS") {
    n_before_uss <- nrow(master)
    master <- master %>%
        filter(!is.na(methodestgest) &
               toupper(trimws(methodestgest)) %in% c("USS", "EUSS"))
    cat(sprintf("  USS filter (master): %d -> %d rows\n",
                n_before_uss, nrow(master)))
}

if (!"match_type" %in% names(master)) {
    stop("master_joined files have no match_type column; cannot derive the ",
         "no-matched-outcome step. Check the data files in ../Data.")
}

n_total_adm    <- nrow(master)
n_total_kch    <- sum(master$facility == "KCH",  na.rm = TRUE)
n_total_smch   <- sum(master$facility == "SMCH", na.rm = TRUE)
n_unmatched_kch  <- sum(master$facility == "KCH"  & master$match_type == "unmatched", na.rm = TRUE)
n_unmatched_smch <- sum(master$facility == "SMCH" & master$match_type == "unmatched", na.rm = TRUE)

cat(sprintf("  Total admissions (study sites, period): %d (KCH %d, SMCH %d)\n",
            n_total_adm, n_total_kch, n_total_smch))

# ------------------------------------------------------------------------------
# STEP 1: MATCHED STUDY POPULATION (joined files)
# ------------------------------------------------------------------------------
cat("Loading matched study population via load_nnu_data()...\n")
raw <- load_nnu_data()

n_matched      <- nrow(raw)
n_matched_kch  <- sum(raw$facility == "KCH",  na.rm = TRUE)
n_matched_smch <- sum(raw$facility == "SMCH", na.rm = TRUE)

n_excl_nomatch      <- n_total_adm  - n_matched
n_excl_nomatch_kch  <- n_total_kch  - n_matched_kch
n_excl_nomatch_smch <- n_total_smch - n_matched_smch

if (abs(n_excl_nomatch - (n_unmatched_kch + n_unmatched_smch)) > 5) {
    cat(sprintf(paste0("  NOTE: derived no-match (%d) differs from master ",
                       "unmatched (%d); using derived (total - matched).\n"),
                n_excl_nomatch, n_unmatched_kch + n_unmatched_smch))
}

# ------------------------------------------------------------------------------
# STEP 2-5: BIRTHWEIGHT, GESTATION, SVN CLASSIFICATION, OUTCOME
# ------------------------------------------------------------------------------
proc <- raw %>%
    mutate(
        bw_kg           = bw_kg_canonical(.),
        gestation_weeks = suppressWarnings(as.numeric(gestation)),
        gender_std      = if ("gender_std" %in% names(.)) gender_std
                          else if ("gender" %in% names(.)) gender
                          else NA_character_,
        outcome_norm    = normalise_outcome(neotreeoutcome),
        outcome_bin     = binary_outcome(outcome_norm)
    )

step_bw <- proc %>%
    filter(is.finite(bw_kg), bw_kg >= BW_MIN_KG, bw_kg <= BW_MAX_KG)
n_excl_bw  <- n_matched - nrow(step_bw)

step_ga <- step_bw %>%
    filter(is.finite(gestation_weeks),
           gestation_weeks >= GA_MIN_NNU,
           gestation_weeks <= GA_MAX)
n_excl_ga  <- nrow(step_bw) - nrow(step_ga)

step_svn <- step_ga %>%
    mutate(
        sga_status   = classify_sga(bw_kg, gestation_weeks, gender_std),
        svn_category = make_svn_category(sga_status, gestation_weeks)
    ) %>%
    filter(!is.na(svn_category))
n_excl_svn <- nrow(step_ga) - nrow(step_svn)

step_final <- step_svn %>% filter(!is.na(outcome_bin))
n_excl_out <- nrow(step_svn) - nrow(step_final)

n_excl_tail <- n_excl_svn + n_excl_out
n_after_bw  <- nrow(step_bw)
n_after_ga  <- nrow(step_ga)
n_final     <- nrow(step_final)
n_final_kch  <- sum(step_final$facility == "KCH",  na.rm = TRUE)
n_final_smch <- sum(step_final$facility == "SMCH", na.rm = TRUE)

# Per-site counts for the intermediate steps (so every combined box can show
# its KCH / SMCH split). Pure additions; do not affect any existing counts.
n_after_bw_kch  <- sum(step_bw$facility == "KCH",  na.rm = TRUE)
n_after_bw_smch <- sum(step_bw$facility == "SMCH", na.rm = TRUE)
n_after_ga_kch  <- sum(step_ga$facility == "KCH",  na.rm = TRUE)
n_after_ga_smch <- sum(step_ga$facility == "SMCH", na.rm = TRUE)
n_excl_bw_kch   <- n_matched_kch  - n_after_bw_kch
n_excl_bw_smch  <- n_matched_smch - n_after_bw_smch
n_excl_ga_kch   <- n_after_bw_kch  - n_after_ga_kch
n_excl_ga_smch  <- n_after_bw_smch - n_after_ga_smch

cat("\n--- NNU FLOWCHART SUMMARY (combined) ---\n")
cat(sprintf("Total admissions (KCH + SMCH, period): %d (KCH %d, SMCH %d)\n",
            n_total_adm, n_total_kch, n_total_smch))
cat(sprintf("Excluded - no matched outcome:         %d (%s)\n",
            n_excl_nomatch, pctf(n_excl_nomatch, n_total_adm)))
cat(sprintf("Matched admission-discharge records:   %d (KCH %d, SMCH %d)\n",
            n_matched, n_matched_kch, n_matched_smch))
cat(sprintf("Excluded - missing/invalid BW:         %d\n", n_excl_bw))
cat(sprintf("Excluded - missing/invalid GA:         %d\n", n_excl_ga))
cat(sprintf("Excluded - SGA/indeterminate outcome:  %d\n", n_excl_tail))
cat(sprintf("Analysis population:                   %d (KCH %d, SMCH %d)\n",
            n_final, n_final_kch, n_final_smch))

# ------------------------------------------------------------------------------
# SAVE COUNTS TABLE (combined)
# ------------------------------------------------------------------------------
counts_df <- data.frame(
    Step = c(
        "1. NNU admissions (KCH + SMCH, Jan 2022 - Dec 2025)",
        "1a. KCH (Malawi)",
        "1b. SMCH (Zimbabwe)",
        "2. Excluded: no matched discharge/outcome",
        "3. Matched admission-discharge records",
        "3a. KCH (Malawi)",
        "3b. SMCH (Zimbabwe)",
        "4. Excluded: missing/invalid birthweight",
        "5. After birthweight exclusion",
        "6. Excluded: missing/invalid gestational age",
        "7. After gestational-age exclusion",
        "8. Excluded: SGA not classifiable",
        "9. Excluded: indeterminate outcome",
        "10. Analysis population (SVN-classified, determinate outcome)",
        "10a. KCH (Malawi)",
        "10b. SMCH (Zimbabwe)"
    ),
    n = c(
        n_total_adm, n_total_kch, n_total_smch,
        n_excl_nomatch,
        n_matched, n_matched_kch, n_matched_smch,
        n_excl_bw, n_after_bw,
        n_excl_ga, n_after_ga,
        n_excl_svn, n_excl_out,
        n_final, n_final_kch, n_final_smch
    )
)
write.csv(counts_df, file.path(OUTPUT_DIR_NNU, "figure_01_flowchart_counts.csv"),
          row.names = FALSE)
cat("  Saved: figure_01_flowchart_counts.csv\n")

# ==============================================================================
# GENERIC CONSORT DRAWING FUNCTION (redesigned for clarity)
# nodes : character vector of main-box labels (top -> bottom)
# fills : character vector, same length, values in c("main","final")
# excls : character vector of exclusion-box labels; excls[i] sits in the gap
#         AFTER node i (length = length(nodes) - 1). Use NA for "no box".
# ==============================================================================
# No title / subtitle / caption are drawn - those are added in the manuscript.
# Every box (main and exclusion) carries its own KCH / SMCH split from the
# labels supplied by build_labels().
draw_consort <- function(nodes, fills, excls, name, width = 9.5, height = 11) {

    M  <- length(nodes)
    ym <- seq(0.90, 0.10, length.out = M)   # full canvas (no title/caption)
    h_m <- 0.188                            # main box height
    h_e <- 0.14                             # exclusion box height (4 lines)
    x_spine <- 0.31                         # spine x-centre
    w_m <- 0.55                             # main box width
    x_excl <- 0.80                          # exclusion box x-centre
    w_e <- 0.38                             # exclusion box width

    main_df <- data.frame(
        xc = rep(x_spine, M), yc = ym, w = rep(w_m, M), h = rep(h_m, M),
        lbl = nodes, fill = fills, stringsAsFactors = FALSE
    )

    arrow_df <- data.frame(
        x    = rep(x_spine, M - 1),
        xend = rep(x_spine, M - 1),
        y    = ym[seq_len(M - 1)]     - h_m / 2,
        yend = ym[seq_len(M - 1) + 1] + h_m / 2
    )

    excl_rows  <- which(!is.na(excls))
    excl_df    <- NULL
    side_df    <- NULL
    if (length(excl_rows) > 0) {
        mid_y <- (ym[excl_rows] + ym[excl_rows + 1]) / 2
        excl_df <- data.frame(
            xc = rep(x_excl, length(excl_rows)), yc = mid_y,
            w  = rep(w_e, length(excl_rows)),    h  = rep(h_e, length(excl_rows)),
            lbl = excls[excl_rows], fill = rep("excl", length(excl_rows)),
            stringsAsFactors = FALSE
        )
        # dashed connectors branch off the vertical stem (longer line)
        side_df <- data.frame(
            x = rep(x_spine, length(excl_rows)),
            xend = rep(x_excl - w_e / 2, length(excl_rows)),
            y = mid_y, yend = mid_y
        )
    }

    box_df <- rbind(main_df, excl_df)

    fill_map <- c("main" = "#DDEEFF", "final" = "#C8E6C9", "excl" = "#FFE5E5")
    main_face <- ifelse(main_df$fill == "final", "bold", "plain")

    p <- ggplot() +
        geom_rect(data = box_df,
                  aes(xmin = xc - w / 2, xmax = xc + w / 2,
                      ymin = yc - h / 2, ymax = yc + h / 2, fill = fill),
                  colour = "grey20", linewidth = 0.7) +
        scale_fill_manual(values = fill_map, guide = "none") +
        geom_text(data = main_df,
                  aes(x = xc, y = yc, label = lbl),
                  size = 5.4, lineheight = 1.2, colour = "grey5",
                  fontface = main_face) +
        geom_segment(data = arrow_df,
                     aes(x = x, xend = xend, y = y, yend = yend),
                     arrow = arrow(length = unit(0.26, "cm"), type = "closed"),
                     colour = "grey20", linewidth = 0.7)

    if (!is.null(excl_df)) {
        p <- p +
            geom_text(data = excl_df,
                      aes(x = xc, y = yc, label = lbl),
                      size = 5.4, lineheight = 1.15, colour = "grey15") +
            geom_segment(data = side_df,
                         aes(x = x, xend = xend, y = y, yend = yend),
                         arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
                         colour = "#C0392B", linewidth = 0.6, linetype = "dashed")
    }

    p <- p +
        scale_x_continuous(limits = c(0, 1.0), expand = c(0, 0)) +
        scale_y_continuous(limits = c(0, 1),   expand = c(0, 0)) +
        theme_void() +
        theme(plot.margin = margin(10, 10, 10, 10))

    save_figure(p, name, width = width, height = height)
}

# ------------------------------------------------------------------------------
# BUILD AND DRAW COMBINED FLOWCHART
# ------------------------------------------------------------------------------
build_labels <- function(total, total_kch, total_smch,
                         excl_nomatch, excl_nomatch_kch, excl_nomatch_smch,
                         matched, matched_kch, matched_smch,
                         excl_bw, excl_bw_kch, excl_bw_smch,
                         after_bw, after_bw_kch, after_bw_smch,
                         excl_ga, excl_ga_kch, excl_ga_smch,
                         after_ga, final, final_kch, final_smch,
                         excl_tail, site_split = TRUE) {

    # site line appended to every box; blank for single-site charts (kch = NA)
    sl <- function(kch, smch) {
        if (site_split && !is.na(kch))
            paste0("\n(KCH ", fmt(kch), "; SMCH ", fmt(smch), ")") else ""
    }

    node_total <- paste0("Neonatal unit admissions\n",
                         "KCH + SMCH, Jan 2022 - Dec 2025\n",
                         "N = ", fmt(total), sl(total_kch, total_smch))
    node_match <- paste0("Matched admission-discharge records\n",
                         "n = ", fmt(matched), sl(matched_kch, matched_smch))
    node_bw    <- paste0("After birthweight exclusion\n",
                         "n = ", fmt(after_bw), sl(after_bw_kch, after_bw_smch))
    node_ga    <- paste0("After gestational-age exclusion\n",
                         "n = ", fmt(after_ga), sl(final_kch, final_smch))
    node_final <- paste0("Analysis population\n",
                         "(valid birthweight + gestation, SVN-classified)\n",
                         "n = ", fmt(final), sl(final_kch, final_smch))
    # When tail == 0, the after-GA box IS the analysis population
    node_final_at_ga <- paste0("Analysis population\n",
                         "(valid birthweight + gestation, SVN-classified)\n",
                         "n = ", fmt(after_ga), sl(final_kch, final_smch))

    e_nomatch <- paste0("Excluded: no matched\ndischarge/outcome\n",
                        "n = ", fmt(excl_nomatch), sl(excl_nomatch_kch, excl_nomatch_smch))
    e_bw      <- paste0("Excluded: missing/invalid\nbirthweight\n",
                        "n = ", fmt(excl_bw), sl(excl_bw_kch, excl_bw_smch))
    e_ga      <- paste0("Excluded: missing/invalid\ngestational age\n",
                        "n = ", fmt(excl_ga), sl(excl_ga_kch, excl_ga_smch))
    e_tail    <- paste0("Excluded: SGA not classifiable\nor indeterminate outcome\n",
                        "n = ", fmt(excl_tail))

    if (excl_tail > 0) {
        list(
            nodes = c(node_total, node_match, node_bw, node_ga, node_final),
            fills = c("main", "main", "main", "main", "final"),
            excls = c(e_nomatch, e_bw, e_ga, e_tail)
        )
    } else {
        list(
            nodes = c(node_total, node_match, node_bw, node_final_at_ga),
            fills = c("main", "main", "main", "final"),
            excls = c(e_nomatch, e_bw, e_ga)
        )
    }
}

lab <- build_labels(
    total = n_total_adm, total_kch = n_total_kch, total_smch = n_total_smch,
    excl_nomatch = n_excl_nomatch,
    excl_nomatch_kch = n_excl_nomatch_kch, excl_nomatch_smch = n_excl_nomatch_smch,
    matched = n_matched, matched_kch = n_matched_kch, matched_smch = n_matched_smch,
    excl_bw = n_excl_bw, excl_bw_kch = n_excl_bw_kch, excl_bw_smch = n_excl_bw_smch,
    after_bw = n_after_bw, after_bw_kch = n_after_bw_kch, after_bw_smch = n_after_bw_smch,
    excl_ga = n_excl_ga, excl_ga_kch = n_excl_ga_kch, excl_ga_smch = n_excl_ga_smch,
    after_ga = n_after_ga, final = n_final,
    final_kch = n_final_kch, final_smch = n_final_smch,
    excl_tail = n_excl_tail, site_split = TRUE
)

draw_consort(
    nodes = lab$nodes, fills = lab$fills, excls = lab$excls,
    name = "figure_01_cohort_flowchart",
    width = 9.5, height = 11
)

# ------------------------------------------------------------------------------
# PER-SITE FLOWCHARTS (KCH and SMCH)
# ------------------------------------------------------------------------------
for (site in c("KCH", "SMCH")) {
    site_long <- if (site == "KCH") "KCH (Kamuzu Central Hospital, Malawi)"
                 else "SMCH (Sally Mugabe Central Hospital, Zimbabwe)"
    site_fn <- tolower(site)

    s_total   <- if (site == "KCH") n_total_kch   else n_total_smch
    s_matched <- if (site == "KCH") n_matched_kch else n_matched_smch
    s_nomatch <- s_total - s_matched

    sp_bw  <- step_bw   %>% filter(facility == site)
    sp_ga  <- step_ga   %>% filter(facility == site)
    sp_svn <- step_svn  %>% filter(facility == site)
    sp_fin <- step_final %>% filter(facility == site)

    s_excl_bw   <- s_matched   - nrow(sp_bw)
    s_after_bw  <- nrow(sp_bw)
    s_excl_ga   <- nrow(sp_bw) - nrow(sp_ga)
    s_after_ga  <- nrow(sp_ga)
    s_excl_svn  <- nrow(sp_ga) - nrow(sp_svn)
    s_excl_out  <- nrow(sp_svn) - nrow(sp_fin)
    s_excl_tail <- s_excl_svn + s_excl_out
    s_final     <- nrow(sp_fin)

    counts_s <- data.frame(
        Step = c(
            paste0("1. NNU admissions (", site, ", Jan 2022 - Dec 2025)"),
            "2. Excluded: no matched discharge/outcome",
            "3. Matched admission-discharge records",
            "4. Excluded: missing/invalid birthweight",
            "5. After birthweight exclusion",
            "6. Excluded: missing/invalid gestational age",
            "7. After gestational-age exclusion",
            "8. Excluded: SGA not classifiable",
            "9. Excluded: indeterminate outcome",
            "10. Analysis population"
        ),
        n = c(s_total, s_nomatch, s_matched,
              s_excl_bw, s_after_bw, s_excl_ga, s_after_ga,
              s_excl_svn, s_excl_out, s_final)
    )
    write.csv(counts_s,
              file.path(OUTPUT_DIR_NNU,
                        paste0("figure_01_flowchart_counts_", site_fn, ".csv")),
              row.names = FALSE)
    cat(sprintf("  Saved: figure_01_flowchart_counts_%s.csv\n", site_fn))

    lab_s <- build_labels(
        total = s_total, total_kch = NA, total_smch = NA,
        excl_nomatch = s_nomatch, excl_nomatch_kch = NA, excl_nomatch_smch = NA,
        matched = s_matched, matched_kch = NA, matched_smch = NA,
        excl_bw = s_excl_bw, excl_bw_kch = NA, excl_bw_smch = NA,
        after_bw = s_after_bw, after_bw_kch = NA, after_bw_smch = NA,
        excl_ga = s_excl_ga, excl_ga_kch = NA, excl_ga_smch = NA,
        after_ga = s_after_ga, final = s_final,
        final_kch = NA, final_smch = NA,
        excl_tail = s_excl_tail, site_split = FALSE
    )

    draw_consort(
        nodes = lab_s$nodes, fills = lab_s$fills, excls = lab_s$excls,
        name = paste0("figure_01_cohort_flowchart_", site_fn),
        width = 9.5, height = 10
    )
}

cat("Figure_01_NNU_cohort_flowchart_v2.R complete.\n")
