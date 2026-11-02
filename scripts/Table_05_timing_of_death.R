# ==============================================================================
# Purpose: Table 5 -- timing of neonatal death by SVN phenotype (main) and by
# site (Supplementary Table S2). Intervals: 0-1, 1-3, 3-7, 7-14, 14-21,
# >21 days, plus Unknown age. Reported as n and % of deaths, and as CFR per
# 1,000 admissions in that interval.
#
# Definitions: age at death = (age at admission in hours, capped at 672) / 24
# + (time from admission to death, in days). Population = svn_filtered
# deaths (n=5,344).
#
# Output: table_05_timing_by_site.csv, table_05_timing_by_phenotype.csv,
#         table_05_timing_by_svn_vs_termaga.csv, table_05_timing_by_bw.csv,
#         and a matching *_cfr.csv for each
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
suppressPackageStartupMessages({library(dplyr); library(tidyr); library(lubridate)})
cat("\n=== 46_timing_of_death_by_site.R ===\n")

AGE_MAX_HOURS <- 672
classify_age_at_death <- function(days){
    dplyr::case_when(
        is.na(days) | days < 0 ~ NA_character_,
        days <= 1.0 ~ "0-1 days",
        days > 1.0 & days <= 3.0 ~ "1-3 days",
        days > 3.0 & days <= 7.0 ~ "3-7 days",
        days > 7.0 & days <= 14.0 ~ "7-14 days",
        days > 14.0 & days <= 21.0 ~ "14-21 days",
        days > 21.0 ~ ">21 days")
}
INTERVALS <- c("0-1 days","1-3 days","3-7 days","7-14 days","14-21 days",">21 days","Unknown age")

nnu <- process_nnu_data(load_nnu_data())
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
# Cohort = all admissions in the four SVN phenotypes (denominator for CFR).
cohort <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(
        facility = as.character(facility),
        svn_category = factor(as.character(svn_category), levels=PH),
        svn_group = ifelse(as.character(svn_category)=="Term-AGA","Term-AGA","All-SVN"),
        bw_kg = bw_kg_canonical(.),
        bw_category = factor(ifelse(is.finite(bw_kg) & bw_kg>=BW_MIN_KG & bw_kg<=BW_MAX_KG,
                                    bw_category_from_kg(bw_kg), NA_character_),
                             levels=c("ELBW","VLBW","LBW","NBW","HBW")))

# Deaths = the fatal subset, with postnatal age at death classified into intervals.
deaths <- cohort %>%
    filter(outcome_bin == 1) %>%
    mutate(
        age_adm_h = suppressWarnings(as.numeric(age)),
        age_adm_h = ifelse(is.finite(age_adm_h) & age_adm_h <= AGE_MAX_HOURS, age_adm_h, NA_real_),
        age_adm_d = age_adm_h/24,
        adm_ts = ymd_hms(datetimeadmission, quiet=TRUE),
        dth_ts = ymd_hms(datetimedeath, quiet=TRUE),
        dur_d  = as.numeric(difftime(dth_ts, adm_ts, units="days")),
        age_death_d = age_adm_d + dur_d,
        interval = classify_age_at_death(age_death_d),
        interval = ifelse(is.na(interval), "Unknown age", interval),
        interval = factor(interval, levels=INTERVALS))

# Admissions denominators for CFR (deaths in interval / admissions in group * 1000).
adm_site <- c(KCH = sum(cohort$facility=="KCH"), SMCH = sum(cohort$facility=="SMCH"))
adm_phen <- setNames(as.integer(table(cohort$svn_category)), levels(cohort$svn_category))
adm_grp2 <- c("All-SVN" = sum(cohort$svn_group=="All-SVN"),
              "Term-AGA" = sum(cohort$svn_group=="Term-AGA"))
adm_bw   <- setNames(as.integer(table(cohort$bw_category)), levels(cohort$bw_category))
N_ADM    <- nrow(cohort)

# generic: counts + % (+ CFR per 1,000 admissions when adm is supplied) by a
# grouping column, wide over intervals. adm = named integer vector of admissions
# per group level; CFR = deaths in interval / admissions in group * 1000.
timing_by <- function(df, gcol, adm=NULL){
    tot <- df %>% group_by(.data[[gcol]]) %>% summarise(total=n(), .groups="drop")
    cnt <- df %>% count(.data[[gcol]], interval, name="n") %>%
        complete(.data[[gcol]], interval=INTERVALS, fill=list(n=0)) %>%
        left_join(tot, by=gcol) %>%
        mutate(pct=round(100*n/total,1), n_pct=sprintf("%d (%.1f%%)",n,pct))
    res <- list(
        n_wide   = cnt %>% select(all_of(gcol),interval,n)   %>% pivot_wider(names_from=interval,values_from=n),
        pct_wide = cnt %>% select(all_of(gcol),interval,n_pct) %>% pivot_wider(names_from=interval,values_from=n_pct),
        totals   = tot)
    if (!is.null(adm)) {
        cnt <- cnt %>% mutate(adm_n = as.integer(adm[as.character(.data[[gcol]])]),
                              cfr   = round(1000 * n / adm_n))
        res$cfr_wide <- cnt %>% select(all_of(gcol),interval,cfr) %>%
            pivot_wider(names_from=interval, values_from=cfr)
    }
    res
}

site <- timing_by(deaths %>% mutate(grp=factor(facility,levels=c("KCH","SMCH"))), "grp", adm=adm_site)
both <- deaths %>% count(interval, name="n") %>%
    complete(interval=INTERVALS, fill=list(n=0)) %>%
    mutate(total=nrow(deaths), pct=round(100*n/total,1), grp="Both sites",
           n_pct=sprintf("%d (%.1f%%)",n,pct), cfr=round(1000*n/N_ADM))
site_pct <- bind_rows(
    site$pct_wide %>% rename(group=grp),
    both %>% select(group=grp,interval,n_pct) %>% pivot_wider(names_from=interval,values_from=n_pct))
site_cfr <- bind_rows(
    site$cfr_wide %>% rename(group=grp),
    both %>% select(group=grp,interval,cfr) %>% pivot_wider(names_from=interval,values_from=cfr))
write.csv(site_pct, file.path(OUT2,"table_05_timing_by_site.csv"),     row.names=FALSE)
write.csv(site_cfr, file.path(OUT2,"table_05_timing_by_site_cfr.csv"), row.names=FALSE)

phen <- timing_by(deaths %>% mutate(grp=svn_category), "grp", adm=adm_phen)
write.csv(phen$pct_wide %>% rename(phenotype=grp),
          file.path(OUT2,"table_05_timing_by_phenotype.csv"),     row.names=FALSE)
write.csv(phen$cfr_wide %>% rename(phenotype=grp),
          file.path(OUT2,"table_05_timing_by_phenotype_cfr.csv"), row.names=FALSE)

grp2 <- timing_by(deaths %>% mutate(grp=factor(svn_group,levels=c("All-SVN","Term-AGA"))), "grp", adm=adm_grp2)
write.csv(grp2$pct_wide %>% rename(group=grp),
          file.path(OUT2,"table_05_timing_by_svn_vs_termaga.csv"),     row.names=FALSE)
write.csv(grp2$cfr_wide %>% rename(group=grp),
          file.path(OUT2,"table_05_timing_by_svn_vs_termaga_cfr.csv"), row.names=FALSE)

# Timing by birth-weight category (appendix / discrimination section: same
# distribution shown by birth weight as by SVN phenotype). Drops the few deaths
# with missing/invalid birth weight from this breakdown only.
bw <- timing_by(deaths %>% filter(!is.na(bw_category)) %>% mutate(grp=bw_category), "grp", adm=adm_bw)
write.csv(bw$pct_wide %>% rename(bw_category=grp),
          file.path(OUT2,"table_05_timing_by_bw.csv"),     row.names=FALSE)
write.csv(bw$cfr_wide %>% rename(bw_category=grp),
          file.path(OUT2,"table_05_timing_by_bw_cfr.csv"), row.names=FALSE)

cat(sprintf("Deaths with timing: %d (of %d); Unknown age %d\n",
            sum(deaths$interval!="Unknown age"), nrow(deaths), sum(deaths$interval=="Unknown age")))
cat("Timing by site - % of deaths:\n");              print(as.data.frame(site_pct))
cat("Timing by site - CFR per 1,000 admissions:\n"); print(as.data.frame(site_cfr))
cat("\nSaved *_timing_by_{site,phenotype,bw,svn_vs_termaga}.csv plus matching *_cfr.csv\n")
cat("46_timing_of_death_by_site.R complete.\n")
