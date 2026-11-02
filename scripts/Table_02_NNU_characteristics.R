# ==============================================================================
# Purpose: Table 2 -- neonatal-unit population characteristics by site, with
# n (%), explicit "Not recorded" rows for variables with missingness, and
# column totals.
#
# Definitions: population = svn_filtered with a determinate outcome (current locked n=30,432).
#
# Output: table_02_nnu_characteristics.csv
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
cat("\n=== 48_table3_nnu_characteristics.R ===\n")

nnu <- process_nnu_data(load_nnu_data())
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
d <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(
        facility = factor(as.character(facility), levels=c("KCH","SMCH")),
        svn_category = factor(as.character(svn_category), levels=PH),
        bw_kg = bw_kg_canonical(.),
        bw_category = factor(ifelse(is.finite(bw_kg) & bw_kg>=BW_MIN_KG & bw_kg<=BW_MAX_KG,
                                    bw_category_from_kg(bw_kg), NA_character_),
                             levels=c("ELBW","VLBW","LBW","NBW","HBW")),
        sex_cat = factor(case_when(
            toupper(trimws(as.character(gender_std))) %in% c("M","MALE","BOY") ~ "Male",
            toupper(trimws(as.character(gender_std))) %in% c("F","FEMALE","GIRL") ~ "Female",
            TRUE ~ NA_character_), levels=c("Male","Female")),
        inout_cat = factor(case_when(
            toupper(trimws(as.character(inorout))) %in% c("TRUE","INBORN","1","IN","YES") ~ "Inborn",
            toupper(trimws(as.character(inorout))) %in% c("FALSE","OUTBORN","0","OUT","NO") ~ "Outborn",
            !is.na(admittedfrom) & toupper(trimws(as.character(admittedfrom))) %in% c("LW","LABOUR WARD") ~ "Inborn",
            !is.na(admittedfrom) ~ "Outborn", TRUE ~ NA_character_), levels=c("Inborn","Outborn")),
        temp_val = suppressWarnings(as.numeric(temperature)),
        temp_val = ifelse(temp_val>=30 & temp_val<=45, temp_val, NA_real_),
        temp_cat = factor(case_when(
            is.na(temp_val) ~ NA_character_,
            temp_val < 32.0 ~ "Severe hypothermia (<32)",
            temp_val < 36.0 ~ "Moderate hypothermia (32-35.9)",
            temp_val < 36.5 ~ "Mild hypothermia (36.0-36.4)",
            temp_val <= 37.5 ~ "Normal (36.5-37.5)",
            TRUE ~ "Hyperthermia (>37.5)"),
            levels=c("Normal (36.5-37.5)","Mild hypothermia (36.0-36.4)","Moderate hypothermia (32-35.9)",
                     "Severe hypothermia (<32)","Hyperthermia (>37.5)")),
        resus_cat = resus_category(resus),
        apgar5_num = suppressWarnings(as.numeric(apgar5)),
        apgar5_num = ifelse(apgar5_num>=0 & apgar5_num<=10, apgar5_num, NA_real_),
        apgar5_cat = factor(case_when(
            is.na(apgar5_num) ~ NA_character_,
            apgar5_num <= 3 ~ "Severely depressed (0-3)",
            apgar5_num <= 6 ~ "Moderately depressed (4-6)",
            TRUE ~ "Normal (7-10)"),
            levels=c("Normal (7-10)","Moderately depressed (4-6)","Severely depressed (0-3)")),
        hiv_status = factor(case_when(
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y") &
                toupper(trimws(as.character(hivtestresult)))=="R" ~ "HIV positive",
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y") &
                toupper(trimws(as.character(hivtestresult)))=="NR" ~ "HIV negative",
            toupper(trimws(as.character(mathivtest))) %in% c("TRUE","Y") ~ "Tested, result unknown",
            toupper(trimws(as.character(mathivtest))) %in% c("FALSE","N") ~ "Not tested",
            TRUE ~ NA_character_),
            levels=c("HIV negative","HIV positive","Tested, result unknown","Not tested")),
        age_adm_h = suppressWarnings(as.numeric(age)),
        age_adm_h = ifelse(is.finite(age_adm_h) & age_adm_h<=672, age_adm_h, NA_real_),
        # Maternal factors captured on the NNU admission form (available for the
        # admitted cohort; NOT in the maternal delivery export).
        matage_num = suppressWarnings(as.numeric(matageyrs)),
        matage_cat = factor(case_when(
            is.na(matage_num) ~ NA_character_, matage_num<16 ~ "<16 years",
            matage_num<20 ~ "16-19 years", matage_num<25 ~ "20-24 years",
            matage_num<35 ~ "25-34 years", TRUE ~ ">=35 years"),
            levels=c("<16 years","16-19 years","20-24 years","25-34 years",">=35 years")),
        par_num = suppressWarnings(as.numeric(par)),
        parity = factor(case_when(
            is.na(par_num) ~ NA_character_, par_num==1 ~ "Primiparous (1st)",
            par_num<=4 ~ "Multiparous (2-4)", TRUE ~ "Grand multipara (5+)"),
            levels=c("Primiparous (1st)","Multiparous (2-4)","Grand multipara (5+)")),
        anc_num = suppressWarnings(as.numeric(antenatalcare)),
        anc = anc_factor(antenatalcare, antenatalcare_category, short=TRUE),
        typebirth_cat = factor(case_when(
            is.na(typebirth) | trimws(typebirth)=="" ~ NA_character_,
            toupper(trimws(typebirth)) %in% c("S","SINGLETON") ~ "Singleton",
            toupper(trimws(typebirth)) %in% c("TW1","TW2","TWIN1","TWIN2","TWIN") ~ "Twin",
            toupper(trimws(typebirth)) %in% c("TR1","TR2","TR3","TRIPLET") ~ "Triplet",
            TRUE ~ NA_character_), levels=c("Singleton","Twin","Triplet")),
        # Method of gestational-age estimation.
        gadating_clean = toupper(trimws(as.character(methodestgest))),
        gadating_cat = factor(case_when(
            is.na(methodestgest) |
                gadating_clean %in% c("","NA","U","UNK","UNKNOWN","-6","-7","-8","-9") ~ NA_character_,
            gadating_clean %in% c("USS","EUSS") | grepl("ULTRA|SCAN", gadating_clean) ~ "Ultrasound",
            gadating_clean == "LMP" | grepl("MENSTR", gadating_clean)                 ~ "Last menstrual period",
            grepl("SFH|FUNDAL|FUND|SYMPH", gadating_clean)                            ~ "Fundal height",
            TRUE                                                                       ~ "Other"),
            levels=c("Ultrasound","Last menstrual period","Fundal height","Other")))

fmt_np <- function(n,tot) sprintf("%d (%.1f%%)", n, ifelse(tot>0,100*n/tot,0))

# categorical block with an explicit Not recorded row when any missing
cat_block <- function(var, varname, show_na=TRUE){
    rows <- list()
    lv <- levels(d[[var]])
    for (site in c("KCH","SMCH")) NULL
    tab_kch <- table(factor(d[[var]][d$facility=="KCH"], levels=lv))
    tab_smc <- table(factor(d[[var]][d$facility=="SMCH"], levels=lv))
    tot_kch <- sum(d$facility=="KCH"); tot_smc <- sum(d$facility=="SMCH"); tot_all <- nrow(d)
    for (l in lv){
        nk <- as.integer(tab_kch[l]); ns <- as.integer(tab_smc[l]); nc <- nk+ns
        rows[[length(rows)+1]] <- data.frame(Variable=varname, Category=l,
            KCH=fmt_np(nk,tot_kch), SMCH=fmt_np(ns,tot_smc), Combined=fmt_np(nc,tot_all),
            stringsAsFactors=FALSE)
    }
    if (show_na){
        nak <- sum(is.na(d[[var]]) & d$facility=="KCH")
        nas <- sum(is.na(d[[var]]) & d$facility=="SMCH")
        if (nak+nas > 0)
            rows[[length(rows)+1]] <- data.frame(Variable=varname, Category="Not recorded",
                KCH=fmt_np(nak,tot_kch), SMCH=fmt_np(nas,tot_smc), Combined=fmt_np(nak+nas,tot_all),
                stringsAsFactors=FALSE)
    }
    do.call(rbind, rows)
}
cont_block <- function(var, varname){
    med <- function(x){x<-x[!is.na(x)]; if(length(x)==0) return("-");
        sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x,.25), quantile(x,.75))}
    data.frame(Variable=varname, Category="median (IQR)",
        KCH=med(d[[var]][d$facility=="KCH"]),
        SMCH=med(d[[var]][d$facility=="SMCH"]),
        Combined=med(d[[var]]), stringsAsFactors=FALSE)
}

hdr <- data.frame(Variable="N admissions", Category="",
    KCH=as.character(sum(d$facility=="KCH")),
    SMCH=as.character(sum(d$facility=="SMCH")),
    Combined=as.character(nrow(d)), stringsAsFactors=FALSE)
deaths_row <- data.frame(Variable="Deaths, n (per 1,000)", Category="",
    KCH=sprintf("%d (%.0f)", sum(d$outcome_bin[d$facility=="KCH"]==1), 1000*mean(d$outcome_bin[d$facility=="KCH"]==1)),
    SMCH=sprintf("%d (%.0f)", sum(d$outcome_bin[d$facility=="SMCH"]==1), 1000*mean(d$outcome_bin[d$facility=="SMCH"]==1)),
    Combined=sprintf("%d (%.0f)", sum(d$outcome_bin==1), 1000*mean(d$outcome_bin==1)), stringsAsFactors=FALSE)

out <- rbind(
    hdr, deaths_row,
    cat_block("svn_category","SVN phenotype", show_na=FALSE),
    cat_block("bw_category","Birth-weight category"),
    cat_block("gadating_cat","Method of gestational-age estimation"),
    cat_block("sex_cat","Sex"),
    cat_block("inout_cat","Place of birth"),
    cont_block("age_adm_h","Age at admission (hours)"),
    cat_block("temp_cat","Admission temperature"),
    cat_block("resus_cat","Resuscitation at birth"),
    cat_block("apgar5_cat","Apgar (5 min)"),
    cat_block("hiv_status","Maternal HIV status"),
    cat_block("matage_cat","Maternal age (NNU admission form)"),
    cat_block("parity","Parity (NNU admission form)"),
    cat_block("anc","Antenatal care visits (NNU admission form)"),
    cat_block("typebirth_cat","Type of pregnancy"))
write.csv(out, file.path(OUT2,"table_02_nnu_characteristics.csv"), row.names=FALSE)
cat("Saved table_02_nnu_characteristics.csv (", nrow(out)," rows)\n", sep="")
print(out, row.names=FALSE)
cat("48_table3_nnu_characteristics.R complete.\n")
