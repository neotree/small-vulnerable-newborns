# ==============================================================================
# Purpose: Table 1 -- stillbirth rate (per 1,000 total deliveries) by SVN
# phenotype and by birth-weight category, by site and combined, with 95%
# Wilson CI.
#
# Definitions: base = maternal analysis population (valid birth weight +
# gestation). Stillbirth = neotreeoutcome in {SBF, SBM, STBF, STBM}. Rate =
# stillbirths / all deliveries in the category x 1,000.
#
# Output: table_01b_stillbirth_rate_by_svn.csv, table_01b_stillbirth_rate_by_bw.csv
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
cat("\n=== 49_stillbirth_by_category.R ===\n")

STILLBIRTH_CODES <- c("SBF","SBM","STBF","STBM")
GA_MIN_M <- if (exists("GA_MIN_MAT")) GA_MIN_MAT else 22

wilson <- function(x, n){
    if (n == 0) return(c(NA, NA))
    p <- x/n; z <- 1.96; d <- 1 + z^2/n
    c(1000*((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d),
      1000*((p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2)))/d))
}

proc <- process_maternal_data(load_maternal_data())
base <- proc %>%
    filter(is.finite(bw_kg), bw_kg>=BW_MIN_KG, bw_kg<=BW_MAX_KG,
           is.finite(gestation_weeks), gestation_weeks>=GA_MIN_M, gestation_weeks<=GA_MAX) %>%
    mutate(facility = factor(as.character(facility), levels=c("KCH","SMCH")),
           is_sb = !is.na(neotreeoutcome) & toupper(trimws(neotreeoutcome)) %in% STILLBIRTH_CODES,
           sga_status = classify_sga(bw_kg, gestation_weeks, gender_std),
           svn_category = factor(make_svn_category(sga_status, gestation_weeks),
                                 levels=c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")),
           bw_category = factor(bw_category_from_kg(bw_kg), levels=c("ELBW","VLBW","LBW","NBW","HBW")))

rate_by <- function(df, gcol){
    df %>% group_by(site=facility, cat=.data[[gcol]]) %>%
        summarise(n_births=n(), n_sb=sum(is_sb), .groups="drop") %>%
        bind_rows(df %>% group_by(cat=.data[[gcol]]) %>%
                  summarise(site=factor("Combined"), n_births=n(), n_sb=sum(is_sb), .groups="drop")) %>%
        rowwise() %>%
        mutate(rate_per1000 = round(1000*n_sb/n_births,1),
               ci = list(round(wilson(n_sb,n_births),1)),
               ci_txt = sprintf("%.1f-%.1f", ci[1], ci[2])) %>%
        ungroup() %>%
        transmute(site, category=cat, n_stillbirths=n_sb, n_deliveries=n_births,
                  stillbirth_per_1000=rate_per1000, ci95=ci_txt) %>%
        arrange(category, site)
}

svn_tab <- rate_by(base, "svn_category")
bw_tab  <- rate_by(base, "bw_category")
write.csv(svn_tab, file.path(OUT2,"table_01b_stillbirth_rate_by_svn.csv"), row.names=FALSE)
write.csv(bw_tab,  file.path(OUT2,"table_01b_stillbirth_rate_by_bw.csv"),  row.names=FALSE)

cat(sprintf("Deliveries %d; stillbirths %d (%.1f per 1,000)\n",
            nrow(base), sum(base$is_sb), 1000*mean(base$is_sb)))
cat("Stillbirth per 1,000 by SVN phenotype:\n"); print(as.data.frame(svn_tab), row.names=FALSE)
cat("\nStillbirth per 1,000 by birth weight:\n"); print(as.data.frame(bw_tab), row.names=FALSE)
cat("\nSaved table_01b_stillbirth_rate_by_svn.csv, table_01b_stillbirth_rate_by_bw.csv\n")
cat("49_stillbirth_by_category.R complete.\n")
