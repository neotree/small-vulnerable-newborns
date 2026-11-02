# ==============================================================================
# Purpose: shared setup sourced by every script in this folder -- data
# loading, SVN/SGA classification, birthweight handling, and plotting theme.
#
# Definitions: POPULATION is "all" (default, every newborn with valid GA) or
# "USS" (only newborns with GA estimated by ultrasound scan); set it before
# sourcing this file, or leave unset for "all". Study period 1 January 2022
# to 31 December 2025; facilities KCH (Malawi) and SMCH (Zimbabwe).
#
# Expects, relative to this file's own folder:
#   ../input/   the six source CSVs, birthweight-corrected
#   ../output/  all/{tables,graphs}/ and uss/{tables,graphs}/
#
# Data files (see input/README.md for the full list):
#   NNU, matched only:   zim-matched-nnu-data-*.csv, mwi-matched-nnu-data-*.csv
#   NNU, all admissions: zim-master-data-*.csv, mwi-master-data-*.csv
#                        (used only by Figure_01_cohort_flowchart.R, for the
#                        pre-match denominator)
#   Maternal:            zim-maternal-data-*.csv, mwi-maternal-data-*.csv
# ==============================================================================

suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(tidyr)
    library(ggplot2)
    library(scales)
    library(purrr)
    library(openxlsx)
    library(broom)
    library(stringr)
    library(lubridate)
})

# ==============================================================================
# POPULATION CONFIG
# ==============================================================================
if (!exists("POPULATION")) {
    .pp <- Sys.getenv("POPULATION"); POPULATION <- if (nzchar(.pp)) .pp else "all"
}
stopifnot(POPULATION %in% c("all", "USS"))
POP_FOLDER <- if (POPULATION == "USS") "uss" else "all"   # output subfolder name

cat(sprintf("=== POPULATION: %s (outputs -> outputs/%s/) ===\n", POPULATION, POP_FOLDER))

# ==============================================================================
# PATHS (all relative to this script's own directory, so the repository can
# be moved or copied anywhere and everything still resolves correctly)
# ==============================================================================
# SCRIPT_DIR is set by the small header at the top of every script in this
# folder, BEFORE it sources this file. Fall back to getwd() if sourced some
# other way (e.g. interactively with setwd() already done).
SCRIPT_DIR <- if (exists("SCRIPT_DIR")) SCRIPT_DIR else getwd()

# Data folder: sibling of scripts/. Use the locked, already corrected exports listed in
# config/locked_inputs.csv. This project does not alter the input files.
PROJECT_ROOT <- normalizePath(file.path(SCRIPT_DIR, ".."), mustWork = TRUE)
resolve_project_dir <- function(value, fallback) {
    if (!nzchar(value)) return(file.path(PROJECT_ROOT, fallback))
    value <- path.expand(value)
    if (!grepl("^(/|[A-Za-z]:[/\\\\])", value)) value <- file.path(PROJECT_ROOT, value)
    value
}
DATA_DIR <- resolve_project_dir(Sys.getenv("LBW_INPUT_DIR"), "input")
if (!dir.exists(DATA_DIR))
    stop("Cannot find the data folder at ", DATA_DIR,
         ". Expected input/ next to scripts/ -- see input/README.md for the ",
         "exact files required.")

# --- Data source files ---
SMCH_NNU_FILE       <- file.path(DATA_DIR, "zim-matched-nnu-data-20220101-to-20251231.csv")
KCH_NNU_FILE        <- file.path(DATA_DIR, "mwi-matched-nnu-data-20220101-to-20251231.csv")
SMCH_MATERNAL_FILE  <- file.path(DATA_DIR, "zim-maternal-data-20220101-to-20251231.csv")
KCH_MATERNAL_FILE   <- file.path(DATA_DIR, "mwi-maternal-data-20220101-to-20251231.csv")
SMCH_NNU_MASTER_FILE <- file.path(DATA_DIR, "zim-master-data-20220101-to-20251231.csv")
KCH_NNU_MASTER_FILE  <- file.path(DATA_DIR, "mwi-master-data-20220101-to-20251231.csv")

for (.f in c(SMCH_NNU_FILE, KCH_NNU_FILE, SMCH_MATERNAL_FILE, KCH_MATERNAL_FILE)) {
    if (!file.exists(.f)) stop("Required data file not found: ", .f)
}
rm(.f)

# --- Output directories --------------------------------------------------
# output/<all|uss>/tables/ and output/<all|uss>/graphs/. OUTPUT_DIR_MATERNAL
# and OUTPUT_DIR_NNU are aliases for OUTPUT_DIR (some scripts use one name,
# some the other; all three resolve to the same folder).
OUTPUT_ROOT <- resolve_project_dir(Sys.getenv("LBW_OUTPUT_DIR"), "output")
OUTPUT_DIR_POP      <- file.path(OUTPUT_ROOT, POP_FOLDER)
OUTPUT_DIR          <- file.path(OUTPUT_DIR_POP, "tables")
OUTPUT_DIR_MATERNAL <- OUTPUT_DIR
OUTPUT_DIR_NNU      <- OUTPUT_DIR
OUTPUT_DIR_GRAPHS   <- file.path(OUTPUT_DIR_POP, "graphs")
for (.d in c(OUTPUT_DIR, OUTPUT_DIR_GRAPHS)) {
    if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}
rm(.d)

# ==============================================================================
# STUDY PERIOD AND FACILITY FILTERS
# ==============================================================================
STUDY_START <- as.Date("2022-01-01")
STUDY_END   <- as.Date("2025-12-31")

VALID_FACILITIES <- c("KCH", "SMCH")

# ==============================================================================
# VALIDITY THRESHOLDS
# ==============================================================================
BW_MIN_KG  <- 0.3    # 300 g
BW_MAX_KG  <- 7.0    # 7000 g
GA_MIN_NNU <- 24     # minimum GA for NNU SGA analysis (weeks)
GA_MIN_MAT <- 22     # minimum GA for maternal analysis (weeks)
GA_MAX     <- 42

# ------------------------------------------------------------------------------
# BIRTHWEIGHT-FOR-GESTATION PLAUSIBILITY CHECK (NNU) -- AUDIT ONLY
# Updated 17 September 2026.
#
# No record is excluded on this criterion. Implausible birthweights were
# corrected in the data (input/apply_birthweight_corrections_2026-09-17.R,
# 21 records, corrected to the discharge-form birthweight after manual review);
# the remaining 4 were reviewed and kept as recorded.
#
# process_nnu_data() runs the criterion as an audit: it counts the records that
# meet it and warns about any UID not listed in input/bw_ga_reviewed_keep_uids.txt
# (one UID per line, confidential, not committed), so that a new extract with
# fresh entry errors is noticed.
# ------------------------------------------------------------------------------
IMPLAUSIBLE_BW_TERM_GA_MIN_WEEKS <- 37    # check applies at GA >= 37 weeks
IMPLAUSIBLE_BW_TERM_MAX_KG       <- 1.0   # ... and birthweight < 1,000 g
.reviewed_uid_file <- file.path(DATA_DIR, "bw_ga_reviewed_keep_uids.txt")
BW_GA_REVIEWED_KEEP_UIDS <- if (file.exists(.reviewed_uid_file))
    trimws(readLines(.reviewed_uid_file, warn = FALSE)) else character(0)

is_implausible_bw_for_ga <- function(bw_kg, ga_weeks) {
    is.finite(bw_kg) & is.finite(ga_weeks) &
        ga_weeks >= IMPLAUSIBLE_BW_TERM_GA_MIN_WEEKS &
        bw_kg    <  IMPLAUSIBLE_BW_TERM_MAX_KG
}

# ==============================================================================
# INTERGROWTH-21st 10TH CENTILE REFERENCE (grams)
# Birthweight: published boys/girls values for 24+0 through 42+6.
# Source: official INTERGROWTH-21st Oxford newborn/very-preterm weight tables;
# full day-level transcription in references/official_intergrowth_p10_birthweight_24_42.csv.
# Weeks 22-23 use the published exact-week fetal EFW 10th centiles as
# a study-defined proxy for measured birthweight in the maternal analysis.
# Weeks 14-21 have no lookup because they are outside the study cohort.
# There is no official combined-sex 10th-centile table: the unknown-sex
# proxy is the arithmetic midpoint of the official boys/girls thresholds.
# Do not describe that proxy as an official unisex centile.
# ==============================================================================
PCT10_UNISEX <- c("22"=481, "23"=538)
PCT10_BOYS   <- c("22"=481, "23"=538)
PCT10_GIRLS  <- c("22"=481, "23"=538)

INTERGROWTH_REFERENCE <- read.csv(file.path(SCRIPT_DIR, "..", "references",
    "official_intergrowth_p10_birthweight_24_42.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(INTERGROWTH_REFERENCE) == 266L,
          !anyDuplicated(INTERGROWTH_REFERENCE[c("sex", "week", "day")]))
# Gestation encodes completed weeks and days: 33.4 means 33+4; 33 means 33+0.
# Do not round this notation to the nearest whole week.

# ==============================================================================
# CATEGORY LEVELS (for ordered factors)
# ==============================================================================
BW_LEVELS  <- c("ELBW","VLBW","LBW","NBW","HBW")
SVN_LEVELS <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
SVN_CATS   <- SVN_LEVELS   # alias used by graph scripts

GA_MIN <- GA_MIN_NNU

SVN_COLOURS <- c(
    "Term-AGA"    = "#2166AC",
    "Term-SGA"    = "#74ADD1",
    "Preterm-AGA" = "#D73027",
    "Preterm-SGA" = "#FDAE61"
)
SVN_COLOURS_LABELLED <- c(
    "Term AGA"    = "#2166AC",
    "Term SGA"    = "#74ADD1",
    "Preterm AGA" = "#D73027",
    "Preterm SGA" = "#FDAE61"
)
FACILITY_COLOURS <- c("KCH" = "#1B7837", "SMCH" = "#762A83")
FACILITY_LABELS  <- c(KCH = "KCH (Malawi)", SMCH = "SMCH (Zimbabwe)")

SVN_LABELS <- c(
    "Term-AGA"    = "Term AGA",
    "Term-SGA"    = "Term SGA",
    "Preterm-AGA" = "Preterm AGA",
    "Preterm-SGA" = "Preterm SGA"
)

svn_factor <- function(x) {
    factor(x, levels = SVN_LEVELS, labels = unname(SVN_LABELS))
}

# ==============================================================================
# HELPER FUNCTIONS - BIRTHWEIGHT
# ==============================================================================
# Canonical birthweight is produced UPSTREAM by the Neotree cleaning pipeline
# as a single column, birthweight_g, in GRAMS (coalesce of the legacy
# era/site columns birthweight, bwt, bwtdis, with a kg->g guard). Read
# birthweight_g directly; the fallback below is for pre-fix inputs only.
BW_CANONICAL_COL <- "birthweight_g"

to_kg <- function(x) {
    x_num <- suppressWarnings(as.numeric(x))
    pos   <- x_num[is.finite(x_num) & x_num > 0]
    if (length(pos) > 0) {
        med <- median(pos, na.rm = TRUE)
        if (is.finite(med) && med > 20 && med <= 7000) return(x_num / 1000)
    }
    x_num
}

get_bw_column <- function(df) {
    if (BW_CANONICAL_COL %in% names(df)) return(as.character(df[[BW_CANONICAL_COL]]))
    result <- rep(NA_character_, nrow(df))
    for (col in c("bwt","birthweight","bwtdis")) {
        if (col %in% names(df))
            result <- ifelse(is.na(result) & !is.na(df[[col]]), df[[col]], result)
    }
    result
}

bw_kg_canonical <- function(df) {
    if (BW_CANONICAL_COL %in% names(df)) {
        return(suppressWarnings(as.numeric(df[[BW_CANONICAL_COL]])) / 1000)
    }
    to_kg(get_bw_column(df))
}

# ==============================================================================
# HELPER FUNCTIONS - SGA / SVN CLASSIFICATION
# ==============================================================================

get_10th_pct <- function(ga_weeks, gender_val = NA_character_) {
    if (is.na(ga_weeks) || !is.finite(ga_weeks)) return(NA_real_)
    g <- toupper(trimws(as.character(gender_val)))
    sex <- if (!is.na(g) && g %in% c("M", "MALE", "BOY")) "B" else
           if (!is.na(g) && g %in% c("F", "FEMALE", "GIRL")) "G" else "U"
    week <- floor(ga_weeks)
    day <- round(10 * (ga_weeks - week))
    if (day < 0 || day > 6 || abs(ga_weeks - (week + day/10)) > 1e-6)
        return(NA_real_) # .7-.9 cannot be interpreted as weeks+days
    # The pre-agreed upper bound is 42+0, not the end of week 42.
    if (ga_weeks < 22 || ga_weeks > GA_MAX) return(NA_real_)
    if (week < 24) {
        # The published fetal EFW table supplies exact-week rows, used only
        # as a clearly labelled proxy for maternal deliveries at 22+0/23+0.
        if (day != 0L) return(NA_real_)
        p <- switch(sex, B=PCT10_BOYS, G=PCT10_GIRLS, U=PCT10_UNISEX)
        return(as.numeric(p[as.character(week)]) / 1000)
    }
    ref <- INTERGROWTH_REFERENCE
    boy <- ref$p10_g[ref$sex == "B" & ref$week == week & ref$day == day]
    girl <- ref$p10_g[ref$sex == "G" & ref$week == week & ref$day == day]
    if (length(boy) != 1L || length(girl) != 1L) return(NA_real_)
    (switch(sex, B=boy, G=girl, U=(boy+girl)/2)) / 1000
}

classify_sga <- function(bw_kg, ga_weeks, gender_val = NA_character_) {
    pct10 <- mapply(get_10th_pct, ga_weeks, gender_val, USE.NAMES = FALSE)
    ifelse(is.na(bw_kg) | is.na(pct10), NA_character_,
           ifelse(bw_kg < pct10, "SGA", "AGA"))
}

# Resuscitation at birth (Neotree `resus`): coded dropdown, not free text.
#   Norm / NONE        no resuscitation (reference level)
#   Stim, Suc, O2      stimulation / suction / oxygen, alone or in combination
#   BVM, CPR           bag-valve-mask; chest compressions (highest level wins)
#   UNK or blank       "Unknown" -- an answer option on the form (Malawi), kept
#                      as its own level and never folded into "None"
# Stops on any code it does not recognise, so nothing can fall into a catch-all.
resus_category <- function(resus) {
    x <- toupper(trimws(as.character(resus)))
    unknown <- is.na(resus) | x %in% c("", "NA", "NR", "U", "UNK", "UNKNOWN",
                                       "-6", "-7", "-8", "-9")
    cat <- dplyr::case_when(
        unknown                                   ~ "Unknown",
        x %in% c("NONE", "NORM", "NO", "N", "FALSE", "0") ~ "None",
        grepl("CPR|CARDIO", x)                    ~ "CPR",
        grepl("BVM|BAG.*VALVE|MASK", x)           ~ "BVM",
        grepl("STIM|SUC|O2|OXYGEN", x)            ~ "Stimulation/Suction/O2",
        TRUE                                      ~ NA_character_)
    if (anyNA(cat))
        stop("resus_category(): unrecognised code(s): ",
             paste(unique(x[is.na(cat)]), collapse = ", "))
    factor(cat, levels = c("None", "Stimulation/Suction/O2", "BVM", "CPR", "Unknown"))
}

make_svn_category <- function(sga, ga_weeks) {
    dplyr::case_when(
        is.na(sga) | is.na(ga_weeks)          ~ NA_character_,
        ga_weeks >= 37 & sga == "SGA"         ~ "Term-SGA",
        ga_weeks >= 37 & sga == "AGA"         ~ "Term-AGA",
        ga_weeks <  37 & sga == "SGA"         ~ "Preterm-SGA",
        ga_weeks <  37 & sga == "AGA"         ~ "Preterm-AGA",
        TRUE                                   ~ NA_character_
    )
}

bw_category <- function(bw_kg) {
    dplyr::case_when(
        is.na(bw_kg)          ~ NA_character_,
        bw_kg <  1.000        ~ "ELBW",
        bw_kg <  1.500        ~ "VLBW",
        bw_kg <  2.500        ~ "LBW",
        # HBW is >= 4,000 g (manuscript: NBW 2,500-3,999 g, HBW >= 4,000 g;
        # Neotree bwgroup puts 4,000 g in HBW). Previously "<= 4.000" put
        # babies of exactly 4,000 g in NBW (fix, 17 Sep 2026).
        bw_kg <  4.000        ~ "NBW",
        TRUE                  ~ "HBW"
    )
}

bw_category_from_kg              <- bw_category
sga_category_from_gestational_age <- make_svn_category

# Shared ANC count/category rules, including recovered historical labels.
source(file.path(PROJECT_ROOT, "config", "anc_rules.R"))

# ==============================================================================
# HELPER FUNCTIONS - OUTCOMES
# ==============================================================================

normalise_outcome <- function(x) {
    y <- tolower(trimws(as.character(x)))
    y <- gsub("[^a-z0-9]+", " ", y)
    y <- trimws(gsub("\\s+", " ", y))
    dplyr::case_when(
        grepl("\\bdc\\b|discharged|\\bdcr\\b|\\bdpc\\b|\\bdckmc\\b|kangaroo mother care", y) ~ "Alive",
        grepl("\\babs\\b|abscond|\\bdama\\b|against medical advice", y) ~ "Alive",
        grepl("\\btrh\\b|\\btro\\b|transfer", y) ~ "Alive",
        y == "lb" ~ "Alive",
        grepl("\\bnnd\\b|neonatal death|\\bdda\\b|died during admission", y) ~ "Neonatal Death (in-facility)",
        grepl("\\bbid\\b|brought in dead", y) ~ "Neonatal Death (pre-facility)",
        grepl("\\bennd\\b|early neonatal death", y) ~ "Neonatal Death (in-facility)",
        grepl("\\bstb\\b|\\bsbf\\b|\\bsbm\\b|\\bstbf\\b|\\bstbm\\b|stillbirth|still born", y) ~ "Stillbirth",
        y == "" ~ NA_character_,
        TRUE    ~ "Unknown/Other"
    )
}

binary_outcome <- function(outcome_norm) {
    dplyr::case_when(
        outcome_norm %in% c("Neonatal Death (in-facility)",
                            "Neonatal Death (pre-facility)") ~ 1L,
        outcome_norm == "Alive"                              ~ 0L,
        TRUE                                                  ~ NA_integer_
    )
}

# ==============================================================================
# HELPER FUNCTION - LENGTH OF STAY
# ==============================================================================
MAX_LOS_DAYS <- 180

get_los_days <- function(df) {
    if (all(c("datetimeadmission","datetimedischarge") %in% names(df))) {
        los <- suppressWarnings(as.numeric(difftime(
            lubridate::parse_date_time(df$datetimedischarge,
                                       orders = c("Ymd HMS", "ymd HMS", "Ymd", "ymd"),
                                       quiet = TRUE),
            lubridate::parse_date_time(df$datetimeadmission,
                                       orders = c("Ymd HMS", "ymd HMS", "Ymd", "ymd"),
                                       quiet = TRUE),
            units = "days"
        )))
        los <- ifelse(is.finite(los) & los >= 0 & los <= MAX_LOS_DAYS, los, NA_real_)
        return(los)
    }
    if ("timespent" %in% names(df)) {
        ts <- suppressWarnings(as.numeric(df$timespent))
        los <- ts / 24
        los <- ifelse(is.finite(los) & los >= 0 & los <= MAX_LOS_DAYS, los, NA_real_)
        return(los)
    }
    rep(NA_real_, nrow(df))
}

# ==============================================================================
# HELPER FUNCTION - DATE FILTERING
# ==============================================================================

filter_study_period <- function(df, date_col = "adm_date_parsed") {
    if (!date_col %in% names(df)) {
        fallbacks <- c("dateadmission", "datetimeadmission", "startedat")
        for (fb in fallbacks) {
            if (fb %in% names(df)) {
                date_col <- fb
                break
            }
        }
    }
    if (!date_col %in% names(df)) {
        warning("No date column found for study period filter. Returning all rows.")
        return(df)
    }

    df <- df %>%
        mutate(.adm_date = suppressWarnings(
            as.Date(lubridate::parse_date_time(
                .data[[date_col]],
                orders = c("Ymd HMS", "ymd HMS", "Ymd", "ymd"),
                quiet = TRUE
            ))
        ))

    n_before <- nrow(df)
    df <- df %>%
        filter(!is.na(.adm_date),
               .adm_date >= STUDY_START,
               .adm_date <= STUDY_END)
    n_after <- nrow(df)

    cat(sprintf("  Study period filter (%s to %s): %d -> %d rows (%d excluded)\n",
                STUDY_START, STUDY_END, n_before, n_after, n_before - n_after))
    df
}

# ==============================================================================
# DATA LOADING FUNCTIONS
# ==============================================================================

load_nnu_data <- function(sites = c("KCH","SMCH")) {
    datasets <- list()

    if ("KCH" %in% sites) {
        cat("  Reading KCH NNU data...\n")
        d <- suppressMessages(read_csv(KCH_NNU_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        if (!"facility" %in% names(d)) d <- d %>% mutate(facility = "KCH")
        d <- d %>% filter(facility %in% "KCH")
        datasets[["KCH"]] <- d
    }

    if ("SMCH" %in% sites) {
        cat("  Reading SMCH NNU data...\n")
        d <- suppressMessages(read_csv(SMCH_NNU_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        d <- d %>% filter(facility %in% "SMCH")
        datasets[["SMCH"]] <- d
    }

    combined <- bind_rows(datasets)

    if (!"antenatalcare_category" %in% names(combined)) combined$antenatalcare_category <- NA_character_
    combined$antenatalcare_category <- anc_band(combined$antenatalcare, combined$antenatalcare_category)
    combined <- filter_study_period(combined, date_col = "adm_date_parsed")

    if (POPULATION == "USS") {
        n_before_uss <- nrow(combined)
        combined <- combined %>%
            filter(!is.na(methodestgest) &
                   toupper(trimws(methodestgest)) %in% c("USS", "EUSS"))
        cat(sprintf("  USS filter: %d -> %d rows\n",
                    n_before_uss, nrow(combined)))
    }

    cat(sprintf("  NNU data loaded: %d rows (KCH: %d, SMCH: %d)\n",
                nrow(combined),
                sum(combined$facility == "KCH"),
                sum(combined$facility == "SMCH")))

    combined
}

load_maternal_data <- function(sites = c("KCH","SMCH")) {
    datasets <- list()

    if ("KCH" %in% sites) {
        cat("  Reading KCH maternal data...\n")
        d <- suppressMessages(read_csv(KCH_MATERNAL_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        if (!"facility" %in% names(d)) d <- d %>% mutate(facility = "KCH")
        d <- d %>% filter(facility %in% "KCH")
        datasets[["KCH"]] <- d
    }

    if ("SMCH" %in% sites) {
        cat("  Reading SMCH maternal data...\n")
        d <- suppressMessages(read_csv(SMCH_MATERNAL_FILE, guess_max = 50000,
                                       col_types = cols(.default = "c")))
        d <- d %>% filter(facility %in% "SMCH")
        datasets[["SMCH"]] <- d
    }

    combined <- bind_rows(datasets)

    combined <- filter_study_period(combined, date_col = "dateadmission")

    cat(sprintf("  Maternal data loaded: %d rows (KCH: %d, SMCH: %d)\n",
                nrow(combined),
                sum(combined$facility == "KCH"),
                sum(combined$facility == "SMCH")))

    combined
}

# ==============================================================================
# CORE PROCESSING - NNU DATA
# ==============================================================================

process_nnu_data <- function(raw) {
    gender_col <- if ("gender_std" %in% names(raw)) "gender_std"
                  else if ("gender" %in% names(raw)) "gender"
                  else NULL

    processed <- raw %>%
        mutate(
            bw_kg          = bw_kg_canonical(.),
            gestation_weeks= suppressWarnings(as.numeric(gestation)),
            gender_std     = if (!is.null(gender_col)) .data[[gender_col]] else NA_character_,
            outcome_norm   = normalise_outcome(neotreeoutcome),
            outcome_bin    = binary_outcome(outcome_norm)
        )

    lbw <- processed %>%
        filter(is.finite(bw_kg), bw_kg >= BW_MIN_KG, bw_kg <= BW_MAX_KG) %>%
        mutate(bw_cat = bw_category(bw_kg))

    ga_valid <- lbw %>%
        filter(is.finite(gestation_weeks),
               gestation_weeks >= GA_MIN_NNU,
               gestation_weeks <= GA_MAX)

    # Audit only (see IMPLAUSIBLE_BW_* above): no record is excluded.
    flag <- is_implausible_bw_for_ga(ga_valid$bw_kg, ga_valid$gestation_weeks)
    flagged_bw_ga <- ga_valid[flag, , drop = FALSE]
    unreviewed <- setdiff(flagged_bw_ga$uid, BW_GA_REVIEWED_KEEP_UIDS)
    cat(sprintf(paste0("  Plausibility check (BW < %.0f g at GA >= %d wk, not excluded): ",
                       "%d record(s) retained after review (KCH %d, SMCH %d)\n"),
                IMPLAUSIBLE_BW_TERM_MAX_KG * 1000,
                IMPLAUSIBLE_BW_TERM_GA_MIN_WEEKS,
                nrow(flagged_bw_ga),
                sum(flagged_bw_ga$facility == "KCH"),
                sum(flagged_bw_ga$facility == "SMCH")))
    if (length(unreviewed) > 0)
        warning("Birthweight < 1,000 g at GA >= 37 wk for UID(s) not yet reviewed: ",
                paste(unreviewed, collapse = ", "), call. = FALSE)

    svn <- ga_valid %>%
        mutate(
            sga_status   = classify_sga(bw_kg, gestation_weeks, gender_std),
            svn_category = make_svn_category(sga_status, gestation_weeks)
        ) %>%
        filter(!is.na(svn_category))

    analysis <- svn %>%
        filter(!is.na(outcome_bin))

    list(raw = raw, processed = processed,
         lbw_filtered = lbw, ga_filtered = ga_valid,
         flagged_bw_ga_retained = flagged_bw_ga,
         svn_filtered = svn, analysis = analysis)
}

# ==============================================================================
# CORE PROCESSING - MATERNAL DATA
# ==============================================================================

process_maternal_data <- function(raw) {
    gender_col <- if ("sexdis" %in% names(raw)) "sexdis"
                  else if ("gender_std" %in% names(raw)) "gender_std"
                  else if ("gender" %in% names(raw)) "gender"
                  else NULL
    raw %>%
        mutate(
            bw_kg          = bw_kg_canonical(.),
            gestation_weeks= suppressWarnings(as.numeric(gestation)),
            gender_std     = if (!is.null(gender_col)) .data[[gender_col]] else NA_character_,
            outcome_norm   = normalise_outcome(neotreeoutcome),
            matageyrs_num  = suppressWarnings(as.numeric(matageyrs))
        )
}

# ==============================================================================
# GGPLOT2 THEME - publication-ready
# ==============================================================================
theme_svn <- function(base_size = 11, legend_position = "right") {
    theme_classic(base_size = base_size) %+replace%
        theme(
            axis.line   = element_line(colour = "grey30", linewidth = 0.5),
            axis.ticks  = element_line(colour = "grey50", linewidth = 0.4),
            axis.text   = element_text(colour = "grey20"),
            axis.title  = element_text(colour = "grey10", face = "bold"),
            panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", colour = NA),
            legend.position  = legend_position,
            legend.title     = element_text(face = "bold", size = rel(0.9)),
            legend.text      = element_text(size = rel(0.85)),
            legend.key.size  = unit(0.45, "cm"),
            legend.background = element_rect(fill = "white", colour = "grey80"),
            strip.background = element_rect(fill = "grey94", colour = "grey70"),
            strip.text       = element_text(face = "bold", colour = "grey15"),
            plot.title       = element_text(face = "bold", hjust = 0,
                                            size = rel(1.15), colour = "grey10",
                                            margin = margin(b = 4)),
            plot.subtitle    = element_text(hjust = 0, colour = "grey30",
                                            size = rel(0.95), margin = margin(b = 6)),
            plot.caption     = element_text(hjust = 1, colour = "grey50",
                                            size = rel(0.75), margin = margin(t = 6)),
            plot.margin      = margin(10, 14, 8, 10)
        )
}

theme_set(theme_svn())

# ==============================================================================
# FIGURE SAVING
# ==============================================================================

save_figure <- function(plot, name, width = 10, height = 7, dpi = 300,
                        outdir = OUTPUT_DIR_GRAPHS) {
    pdf_file <- file.path(outdir, paste0(name, ".pdf"))
    png_file <- file.path(outdir, paste0(name, ".png"))
    ggsave(pdf_file, plot = plot, width = width, height = height, device = "pdf")
    ggsave(png_file, plot = plot, width = width, height = height, dpi = dpi, device = "png")
    cat(sprintf("  Saved: %s (PDF + PNG) [%.0f x %.0f in]\n", name, width, height))
}

# ==============================================================================
# WILSON 95% CI FOR A PROPORTION
# ==============================================================================

prop_ci_wilson <- function(x, n, conf = 0.95) {
    if (is.na(n) || n == 0) return(c(lower = NA_real_, upper = NA_real_))
    z <- qnorm(1 - (1 - conf) / 2)
    p <- x / n
    d <- 1 + z^2 / n
    ctr <- (p + z^2 / (2 * n)) / d
    marg <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
    c(lower = max(0, ctr - marg), upper = min(1, ctr + marg))
}

prop_ci <- function(x, n, conf = 0.95) {
    res <- prop_ci_wilson(x, n, conf)
    list(lower = unname(res["lower"]), upper = unname(res["upper"]))
}

add_prop_ci <- function(df, n_col = "deaths", total_col = "total_admissions",
                        rate_col = "mortality_rate") {
    df %>%
        rowwise() %>%
        mutate(
            ci_res   = list(prop_ci_wilson(.data[[n_col]], .data[[total_col]])),
            ci_lower = ci_res[["lower"]] * 100,
            ci_upper = ci_res[["upper"]] * 100
        ) %>%
        ungroup() %>%
        select(-ci_res)
}

# ==============================================================================
# STANDARD CAPTIONS
# ==============================================================================

svn_caption <- function(extra = NULL) {
    base <- paste0(
        "KCH = Kamuzu Central Hospital, Malawi; ",
        "SMCH = Sally Mugabe Central Hospital, Zimbabwe. 2022-2025.")
    pop_note <- if (POPULATION == "USS") " [USS-confirmed GA subset]" else ""
    full <- paste0(base, pop_note)
    if (!is.null(extra)) paste0(full, "\n", extra) else full
}

# Suppress default Rplots.pdf generation
if (capabilities("cairo")) options(bitmapType = "cairo")
pdf(file = NULL)
options(device = function() { pdf(file = NULL) })

cat(sprintf("00a_setup_shared.R loaded. Data: %s | Output: %s\n", DATA_DIR, OUTPUT_DIR_POP))
