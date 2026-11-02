# ==============================================================================
# Purpose: Table 3 -- case fatality (deaths per 1,000 admissions) by SVN
# phenotype and by birth-weight category, by site and combined, with 95% CI.
# Table 3 also carries a pooled "All SVN" row (Term-SGA + Preterm-AGA +
# Preterm-SGA) in the by-phenotype CFR.
# Table 4 -- distribution of admissions by SVN phenotype x birth-weight
# category (both sites combined), n and % of all admissions.
#
# Definitions: population = svn_filtered. Respects POPULATION ("all" default,
# "USS" for the ultrasound subset).
#
# Output: table_04_distribution_crosstab.csv (counts, + row/col totals),
#         table_03_cfr_by_svn.csv, table_03_cfr_by_bw.csv (CFR per 1,000)
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

suppressPackageStartupMessages({library(dplyr); library(tidyr)})
cat("\n=== 53_distribution_and_cfr.R ===\n")
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
BW <- c("ELBW","VLBW","LBW","NBW","HBW")

wilson <- function(x, n){
    if (n == 0) return(c(NA, NA))
    p <- x/n; z <- 1.96; d <- 1 + z^2/n
    lo <- (p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d
    hi <- (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d
    round(c(1000*lo, 1000*hi), 1)
}

nnu <- process_nnu_data(load_nnu_data())
d <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(facility = factor(as.character(facility), levels=c("KCH","SMCH")),
           svn_category = factor(as.character(svn_category), levels=PH),
           bw_kg = bw_kg_canonical(.),
           bw_category = factor(ifelse(is.finite(bw_kg)&bw_kg>=BW_MIN_KG&bw_kg<=BW_MAX_KG,
                                       bw_category_from_kg(bw_kg), NA_character_), levels=BW))
GRAND <- nrow(d)

# ---- Table 4: distribution crosstab (both sites), counts + row/col totals ----
ct <- d %>% filter(!is.na(bw_category)) %>%
    count(bw_category, svn_category) %>%
    complete(bw_category=BW, svn_category=PH, fill=list(n=0)) %>%
    pivot_wider(names_from=svn_category, values_from=n, values_fill=0) %>%
    arrange(factor(bw_category, levels=BW))
ct <- ct %>% mutate(`Row total` = rowSums(across(all_of(PH))))
col_tot <- ct %>% summarise(across(all_of(c(PH,"Row total")), sum)) %>%
    mutate(bw_category="Column total") %>% select(bw_category, everything())
ct_out <- bind_rows(ct %>% mutate(bw_category=as.character(bw_category)), col_tot)
write.csv(ct_out, file.path(OUT2,"table_04_distribution_crosstab.csv"), row.names=FALSE)

# ---- Table 5 / 5b: CFR per 1,000 by category and site ------------------------
cfr_by <- function(df, gcol, levels){
    f <- function(dd, site) dd %>% group_by(category=.data[[gcol]]) %>%
        summarise(n=n(), deaths=sum(outcome_bin==1), .groups="drop") %>%
        rowwise() %>%
        mutate(site=site, per1000=round(1000*deaths/n,1),
               ci=list(wilson(deaths,n)), ci95=sprintf("%.1f-%.1f", ci[1], ci[2])) %>%
        ungroup() %>% select(site, category, n, deaths, per1000, ci95)
    bind_rows(f(filter(df,facility=="KCH"),"KCH"),
              f(filter(df,facility=="SMCH"),"SMCH"),
              f(df,"Combined")) %>%
        filter(!is.na(category)) %>%
        mutate(category=factor(category, levels=levels)) %>% arrange(category, site)
}
cfr_svn <- cfr_by(d, "svn_category", PH)
cfr_pooled <- cfr_by(d %>% filter(svn_category != "Term-AGA") %>% mutate(grp = "All SVN"),
                     "grp", "All SVN")
cfr_svn <- bind_rows(mutate(cfr_svn, category = as.character(category)),
                     mutate(cfr_pooled, category = as.character(category)))
cfr_bw  <- cfr_by(filter(d, !is.na(bw_category)), "bw_category", BW)
write.csv(cfr_svn, file.path(OUT2,"table_03_cfr_by_svn.csv"), row.names=FALSE)
write.csv(cfr_bw,  file.path(OUT2,"table_03_cfr_by_bw.csv"),  row.names=FALSE)

cat(sprintf("Population %d; deaths %d (%.0f per 1,000)\n", GRAND, sum(d$outcome_bin==1),
            1000*mean(d$outcome_bin==1)))
cat("\nDistribution crosstab (BW x phenotype):\n"); print(as.data.frame(ct_out), row.names=FALSE)
cat("\nCFR per 1,000 by SVN phenotype:\n"); print(as.data.frame(cfr_svn), row.names=FALSE)
cat("\nCFR per 1,000 by birth weight:\n"); print(as.data.frame(cfr_bw), row.names=FALSE)
cat("\nSaved table_04_distribution_crosstab.csv, table_03_cfr_by_svn.csv, table_03_cfr_by_bw.csv\n")
cat("53_distribution_and_cfr.R complete.\n")
