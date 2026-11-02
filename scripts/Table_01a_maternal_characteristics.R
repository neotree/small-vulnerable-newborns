# ==============================================================================
# Purpose: Table 1 -- maternal delivery characteristics by site.
#
# Definitions: the maternal (delivery-outcome) datasets are thinner than the
# NNU datasets. Columns present across both sites are limited to maternal
# age, mode of delivery, type of pregnancy, birth weight and gestation
# (parity, gravidity and antenatal-care counts are not in the maternal
# export; maternal HIV and pregnancy conditions are recorded for KCH/Malawi
# only). This script reports the reliably cross-site variables, plus SVN
# phenotype among live births, as n (%) or median (IQR). Base = maternal
# analysis population (valid birth weight + gestation).
#
# Output: table_01a_maternal_characteristics.csv
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
cat("\n=== 47_table1_maternal_expanded.R ===\n")

STILLBIRTH_CODES <- c("SBF","SBM","STBF","STBM")
GA_MIN_M <- if (exists("GA_MIN_MAT")) GA_MIN_MAT else 22

proc <- process_maternal_data(load_maternal_data())
base <- proc %>%
    filter(is.finite(bw_kg), bw_kg>=BW_MIN_KG, bw_kg<=BW_MAX_KG,
           is.finite(gestation_weeks), gestation_weeks>=GA_MIN_M, gestation_weeks<=GA_MAX) %>%
    mutate(facility = factor(as.character(facility), levels=c("KCH","SMCH")),
           is_stillbirth = !is.na(neotreeoutcome) & toupper(trimws(neotreeoutcome)) %in% STILLBIRTH_CODES,
           matage = matageyrs_num)
has <- function(c) c %in% names(base)

# --- derived factors (only when source column exists) -------------------------
base$matage_cat <- factor(dplyr::case_when(
    is.na(base$matage) ~ NA_character_, base$matage<20 ~ "<20 years",
    base$matage<35 ~ "20-34 years", TRUE ~ ">=35 years"),
    levels=c("<20 years","20-34 years",">=35 years"))

if (has("modedelivery")) {
    m <- suppressWarnings(as.numeric(base$modedelivery))
    base$mode_del <- factor(dplyr::case_when(
        is.na(m) ~ NA_character_, m==1 ~ "SVD", m %in% c(2,3,6,7) ~ "Assisted vaginal",
        m==4 ~ "Elective CS", m %in% c(5,8) ~ "Emergency CS", TRUE ~ NA_character_),
        levels=c("SVD","Assisted vaginal","Elective CS","Emergency CS"))
}
if (has("typebirth")) {
    tb <- toupper(trimws(as.character(base$typebirth)))
    base$typebirth_cat <- factor(dplyr::case_when(
        is.na(base$typebirth) | tb=="" ~ NA_character_,
        tb %in% c("S","SINGLETON") ~ "Singleton",
        tb %in% c("TW1","TW2","TWIN1","TWIN2","TWIN") ~ "Twin",
        tb %in% c("TR1","TR2","TR3","TRIPLET") ~ "Triplet", TRUE ~ NA_character_),
        levels=c("Singleton","Twin","Triplet"))
}
# SVN phenotype among live births (classification as in 02b)
base$sga_status  <- classify_sga(base$bw_kg, base$gestation_weeks, base$gender_std)
base$svn_category<- factor(make_svn_category(base$sga_status, base$gestation_weeks),
                           levels=c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA"))

TK <- sum(base$facility=="KCH"); TS <- sum(base$facility=="SMCH"); TA <- nrow(base)
LK <- sum(!base$is_stillbirth & base$facility=="KCH")   # live births per site (SVN denom)
LS <- sum(!base$is_stillbirth & base$facility=="SMCH"); LA <- sum(!base$is_stillbirth)
fmt_np <- function(n,tot) sprintf("%d (%.1f%%)", n, ifelse(tot>0,100*n/tot,0))
row_hdr <- function(name,k,s,a) data.frame(Variable=name,Category="",KCH=k,SMCH=s,Combined=a,stringsAsFactors=FALSE)

cat_block <- function(var,varname,denomK=TK,denomS=TS,denomA=TA,subset=NULL,show_na=TRUE){
    dat <- base; if(!is.null(subset)) dat <- base[subset,]
    lv <- levels(dat[[var]]); out <- list()
    tk<-table(factor(dat[[var]][dat$facility=="KCH"],levels=lv))
    ts<-table(factor(dat[[var]][dat$facility=="SMCH"],levels=lv))
    for (l in lv) out[[length(out)+1]] <- data.frame(Variable=varname,Category=l,
        KCH=fmt_np(as.integer(tk[l]),denomK),SMCH=fmt_np(as.integer(ts[l]),denomS),
        Combined=fmt_np(as.integer(tk[l])+as.integer(ts[l]),denomA),stringsAsFactors=FALSE)
    if (show_na){nak<-sum(is.na(dat[[var]])&dat$facility=="KCH");nas<-sum(is.na(dat[[var]])&dat$facility=="SMCH")
        if(nak+nas>0) out[[length(out)+1]] <- data.frame(Variable=varname,Category="Not recorded",
            KCH=fmt_np(nak,denomK),SMCH=fmt_np(nas,denomS),Combined=fmt_np(nak+nas,denomA),stringsAsFactors=FALSE)}
    do.call(rbind,out)
}
med_iqr <- function(x){x<-x[!is.na(x)]; if(length(x)==0) return("-");
    sprintf("%.1f (%.1f-%.1f)",median(x),quantile(x,.25),quantile(x,.75))}
cont_block <- function(var,varname) data.frame(Variable=varname,Category="median (IQR)",
    KCH=med_iqr(base[[var]][base$facility=="KCH"]),SMCH=med_iqr(base[[var]][base$facility=="SMCH"]),
    Combined=med_iqr(base[[var]]),stringsAsFactors=FALSE)
mean_sd <- function(x){x<-x[!is.na(x)]; sprintf("%.2f (%.2f)",mean(x),sd(x))}

blocks <- list(
    row_hdr("Deliveries (valid BW + GA), n", as.character(TK), as.character(TS), as.character(TA)),
    row_hdr("Live births, n", as.character(LK), as.character(LS), as.character(LA)),
    row_hdr("Stillbirths, n (per 1,000 deliveries)",
        sprintf("%d (%.0f)",sum(base$is_stillbirth&base$facility=="KCH"),1000*mean(base$is_stillbirth[base$facility=="KCH"])),
        sprintf("%d (%.0f)",sum(base$is_stillbirth&base$facility=="SMCH"),1000*mean(base$is_stillbirth[base$facility=="SMCH"])),
        sprintf("%d (%.0f)",sum(base$is_stillbirth),1000*mean(base$is_stillbirth))),
    cont_block("matage","Maternal age (years)"),
    cat_block("matage_cat","Maternal age group"))
if (has("modedelivery")) blocks <- c(blocks, list(cat_block("mode_del","Mode of delivery")))
if (has("typebirth"))    blocks <- c(blocks, list(cat_block("typebirth_cat","Type of pregnancy")))
blocks <- c(blocks, list(
    row_hdr("Birth weight kg, mean (SD)", mean_sd(base$bw_kg[base$facility=="KCH"]),
        mean_sd(base$bw_kg[base$facility=="SMCH"]), mean_sd(base$bw_kg)),
    row_hdr("Low birth weight (<2.5 kg)", fmt_np(sum(base$bw_kg<2.5&base$facility=="KCH"),TK),
        fmt_np(sum(base$bw_kg<2.5&base$facility=="SMCH"),TS), fmt_np(sum(base$bw_kg<2.5),TA)),
    cont_block("gestation_weeks","Gestation (weeks)"),
    row_hdr("Preterm (<37 weeks)", fmt_np(sum(base$gestation_weeks<37&base$facility=="KCH"),TK),
        fmt_np(sum(base$gestation_weeks<37&base$facility=="SMCH"),TS), fmt_np(sum(base$gestation_weeks<37),TA)),
    # SVN phenotype among LIVE BIRTHS (denominator = live births)
    cat_block("svn_category","SVN phenotype (live births)", denomK=LK,denomS=LS,denomA=LA,
              subset=!base$is_stillbirth, show_na=FALSE)))

out <- do.call(rbind, blocks)
write.csv(out, file.path(OUT2,"table_01a_maternal_characteristics.csv"), row.names=FALSE)
cat(sprintf("Maternal base: %d (KCH %d, SMCH %d); live births %d\n", TA, TK, TS, LA))
cat("NOTE: parity, gravidity, ANC not in maternal export; HIV/pregconditions KCH-only (omitted for cross-site consistency).\n")
cat("Saved table_01a_maternal_characteristics.csv (", nrow(out)," rows)\n", sep="")
print(out, row.names=FALSE)
cat("47_table1_maternal_expanded.R complete.\n")
