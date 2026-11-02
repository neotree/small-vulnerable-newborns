# ==============================================================================
# Purpose: Table 9 -- categorical admission-temperature odds ratios for
# neonatal death, reference Normal (36.5-37.5 C), overall and by SVN
# phenotype.
#
# Definitions: univariable (death ~ temp_cat), pooled sites. Respects
# POPULATION ("all" default, "USS" for the ultrasound subset).
#
# Output: table_09_temp_cat_or_overall.csv, table_09_temp_cat_or_by_phenotype.csv,
#         table_09_temp_cat_counts.csv
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

suppressPackageStartupMessages({library(dplyr); library(broom)})
cat("\n=== 52_temperature_categorical_OR.R ===\n")
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
TLV <- c("Normal (36.5-37.5)","Mild hypothermia (36.0-36.4)","Moderate hypothermia (32-35.9)",
         "Severe hypothermia (<32)","Hyperthermia (>37.5)")

nnu <- process_nnu_data(load_nnu_data())
d <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(svn_category = factor(as.character(svn_category), levels=PH),
           tv = suppressWarnings(as.numeric(temperature)),
           tv = ifelse(tv>=30 & tv<=45, tv, NA_real_),
           temp_cat = factor(case_when(
               is.na(tv) ~ NA_character_,
               tv < 32.0 ~ "Severe hypothermia (<32)",
               tv < 36.0 ~ "Moderate hypothermia (32-35.9)",
               tv < 36.5 ~ "Mild hypothermia (36.0-36.4)",
               tv <= 37.5 ~ "Normal (36.5-37.5)",
               TRUE ~ "Hyperthermia (>37.5)"), levels=TLV))

or_temp <- function(df, label){
    df2 <- df %>% filter(!is.na(temp_cat)) %>% mutate(temp_cat=droplevels(temp_cat))
    if (nlevels(df2$temp_cat) < 2 || length(unique(df2$outcome_bin)) < 2)
        return(tibble(group=label, temp_cat=NA, OR=NA, ci=NA, p=NA, note="insufficient"))
    fit <- glm(outcome_bin ~ temp_cat, data=df2, family=binomial())
    res <- broom::tidy(fit, conf.int=TRUE, exponentiate=TRUE) %>% filter(term!="(Intercept)")
    ref <- levels(df2$temp_cat)[1]
    # Categories with no deaths or no survivors cannot be estimated (separation)
    sep <- df2 %>% group_by(temp_cat) %>%
        summarise(ev=sum(outcome_bin==1), ne=sum(outcome_bin==0), .groups="drop") %>%
        filter(ev==0 | ne==0) %>% pull(temp_cat) %>% as.character()
    bind_rows(tibble(term=paste0("temp_cat",ref), estimate=1, conf.low=NA, conf.high=NA, p.value=NA), res) %>%
        transmute(group=label,
                  temp_cat=sub("^temp_cat","",term),
                  not_est=temp_cat %in% sep,
                  OR=ifelse(not_est, NA, round(estimate,2)),
                  ci=ifelse(not_est, "not estimable",
                            ifelse(is.na(conf.low),"ref",sprintf("%.2f-%.2f",conf.low,conf.high))),
                  p=ifelse(not_est | is.na(p.value),NA,signif(p.value,3)),
                  note=ifelse(not_est, "no deaths or no survivors in category", "")) %>%
        select(-not_est)
}

overall <- or_temp(d, "All phenotypes")
byph <- bind_rows(lapply(PH, function(ph) or_temp(filter(d, svn_category==ph), ph)))
counts <- d %>% filter(!is.na(temp_cat)) %>%
    group_by(temp_cat) %>% summarise(n=n(), deaths=sum(outcome_bin==1),
                                     cfr_per1000=round(1000*mean(outcome_bin==1),1), .groups="drop")
write.csv(overall, file.path(OUT2,"table_09_temp_cat_or_overall.csv"), row.names=FALSE)
write.csv(byph,    file.path(OUT2,"table_09_temp_cat_or_by_phenotype.csv"), row.names=FALSE)
write.csv(counts,  file.path(OUT2,"table_09_temp_cat_counts.csv"), row.names=FALSE)

cat("Temperature category ORs vs Normal (all phenotypes):\n"); print(as.data.frame(overall), row.names=FALSE)
cat("\nSaved table_09_temp_cat_or_overall.csv, table_09_temp_cat_or_by_phenotype.csv, table_09_temp_cat_counts.csv\n")
cat("52_temperature_categorical_OR.R complete.\n")
