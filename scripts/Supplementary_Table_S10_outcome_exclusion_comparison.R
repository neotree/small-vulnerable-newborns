# ==============================================================================
# Supplementary Table S10: admission characteristics by outcome linkage.
#
# Place this file beside 00a_setup_shared.R in the corrected scripts/ folder.
# Run: Rscript Supplementary_Table_S10_outcome_exclusion_comparison.R
# It also runs from any working directory when given the full script path.
#
# Dependencies: the existing corrected 00a_setup_shared.R, its reference files,
# and the same confidential input/ CSVs as the manuscript analysis. Packages
# are loaded by the shared setup; no packages or patient data are downloaded.
#
# Outputs: output/all/tables/supplementary_s10/
#   supplementary_outcome_exclusion_comparison.csv
#   Supplementary_Table_outcome_exclusions.xlsx (including source code)
#   R_session_info.txt
#
# The manuscript uses pooled, KCH and SMCH panels. For the site-specific Word
# panels, redundant Site: KCH / Site: SMCH rows are omitted from the CSV.
# No patient-level identifiers are exported. Percentages use explicit available
# denominators. Missing outcomes are not imputed; no equivalence is asserted.
# Cohort assertions intentionally pin this script to the current manuscript
# extract; investigate any assertion failure before updating published values.
# ==============================================================================
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grepl("^--file=", args)]
if (length(file_arg)) {
    SCRIPT_FILE <- normalizePath(sub("^--file=", "", file_arg[1]))
} else {
    source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(f) f$ofile))
    SCRIPT_FILE <- if (length(source_files)) normalizePath(tail(source_files, 1)[[1]])
                   else normalizePath("Supplementary_Table_S10_outcome_exclusion_comparison.R")
}
SCRIPT_DIR <- dirname(SCRIPT_FILE)
POPULATION <- "all"  # S10 compares all admissions; never restrict to USS.
source(file.path(SCRIPT_DIR, "00a_setup_shared.R"))
out <- file.path(OUTPUT_DIR, "supplementary_s10")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
master <- bind_rows(lapply(c(KCH_NNU_MASTER_FILE,SMCH_NNU_MASTER_FILE), function(p)
  suppressMessages(read_csv(p, col_types=cols(.default='c'))))) %>%
  filter(facility %in% VALID_FACILITIES) %>% filter_study_period()
matched <- load_nnu_data()
final <- process_nnu_data(matched)$analysis
key <- function(d) paste(d$facility,d$uid,sep=':')
stopifnot(!anyDuplicated(key(master)), !anyDuplicated(key(matched)),
          all(key(matched) %in% key(master)))
excluded <- master[!key(master) %in% key(matched),]
stopifnot(nrow(master)==34628,nrow(matched)==31223,nrow(excluded)==3405,
          nrow(final)==30432,all(excluded$match_type=='unmatched'))
num <- function(x) suppressWarnings(as.numeric(x))
prepare <- function(d){
 d$bw <- bw_kg_canonical(d)*1000
 d$bw[!is.finite(d$bw)|d$bw<300|d$bw>7000] <- NA
 ga <- num(d$gestation); w <- floor(ga); day <- round((ga-w)*10)
 good <- is.finite(ga)&abs((ga-w)*10-day)<1e-6&day>=0&day<=6&ga>=24&ga<=42
 d$ga <- ifelse(good,w+day/7,NA_real_)
 d$temp <- num(d$temperature); d$temp[d$temp<30|d$temp>45] <- NA
 d$apgar <- num(d$apgar5); d$apgar[d$apgar<0|d$apgar>10] <- NA
 d$male <- ifelse(toupper(d$gender)=='M','Male',ifelse(toupper(d$gender)=='F','Female',NA))
 z <- toupper(trimws(d$inorout)); a <- toupper(trimws(d$admittedfrom))
 d$place <- ifelse(z %in% c('TRUE','INBORN','1','IN','YES'),'Inborn',
  ifelse(z %in% c('FALSE','OUTBORN','0','OUT','NO'),'Outborn',
   ifelse(!is.na(a)&a %in% c('LW','LABOUR WARD'),'Inborn',ifelse(!is.na(a),'Outborn',NA))))
 d$hypo <- ifelse(is.na(d$temp),NA,ifelse(d$temp<36.5,'Hypothermia','Other'))
 d$lowapgar <- ifelse(is.na(d$apgar),NA,ifelse(d$apgar<=6,'Apgar 0–6','Other'))
 d$lbw <- ifelse(is.na(d$bw),NA,ifelse(d$bw<2500,'LBW','Other'))
 d$preterm <- ifelse(is.na(d$ga),NA,ifelse(d$ga<37,'Preterm','Other'))
 d$uss <- ifelse(is.na(d$methodestgest)|d$methodestgest %in% c('','NA','U','UNK','UNKNOWN','-6','-7','-8','-9'),NA,
                 ifelse(toupper(d$methodestgest) %in% c('USS','EUSS'),'Ultrasound','Other'))
 d
}
groups <- lapply(list(Excluded=excluded,Matched=matched,Final=final),prepare)
rows <- list()
add <- function(site,label,stat,values){
 rows[[length(rows)+1]] <<- data.frame(Site=site,Characteristic=label,Statistic=stat,
  Excluded=values[1],Matched=values[2],Final=values[3])
}
for(site in c('All sites','KCH','SMCH')){
 g <- lapply(groups,function(d)if(site=='All sites')d else d[d$facility==site,])
 ns <- vapply(g,nrow,integer(1)); add(site,'Admissions','n',as.character(ns))
 for(s in c('KCH','SMCH'))add(site,paste('Site:',s),'n (% of group)',sapply(g,function(d)sprintf('%d (%.1f%%)',sum(d$facility==s),100*mean(d$facility==s))))
 labels <- c(bw='Birthweight (g)',ga='Gestation (weeks)',temp='Admission temperature (°C)',apgar='Five-minute Apgar')
 for(v in names(labels)){
  add(site,labels[v],'Median [Q1, Q3]',sapply(g,function(d){x<-d[[v]];q<-quantile(x,c(.5,.25,.75),na.rm=TRUE);sprintf('%.1f [%.1f, %.1f]',q[1],q[2],q[3])}))
  add(site,labels[v],'Missing/invalid, n (% of group)',sapply(g,function(d)sprintf('%d (%.1f%%)',sum(is.na(d[[v]])),100*mean(is.na(d[[v]])))))
 }
 cats <- c(male='Male',place='Outborn',hypo='Hypothermia',lowapgar='Apgar 0–6',lbw='LBW',preterm='Preterm',uss='Ultrasound')
 for(v in names(cats)){
  add(site,cats[v],'n / available (% available)',sapply(g,function(d){x<-d[[v]];n<-sum(x==cats[v],na.rm=TRUE);den<-sum(!is.na(x));sprintf('%d / %d (%.1f%%)',n,den,100*n/den)}))
  add(site,cats[v],'Missing/invalid, n (% of group)',sapply(g,function(d)sprintf('%d (%.1f%%)',sum(is.na(d[[v]])),100*mean(is.na(d[[v]])))))
 }
}
tab<-bind_rows(rows); write.csv(tab,file.path(out,'supplementary_outcome_exclusion_comparison.csv'),row.names=FALSE,na='')
print(tab[tab$Site=='All sites',],row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
# Publication-ready editable supplementary table, plus the complete audit code.
wb <- openxlsx::createWorkbook(creator='David De Lorenzo Lopez')
for(site in c('All sites','KCH','SMCH')){
 openxlsx::addWorksheet(wb,site)
 z<-tab[tab$Site==site,c('Characteristic','Statistic','Excluded','Matched','Final')]
 names(z)<-c('Characteristic','Summary','Excluded unmatched outcome','All matched admissions','Final analysis admissions')
 openxlsx::writeData(wb,site,z,withFilter=FALSE)
 openxlsx::addStyle(wb,site,openxlsx::createStyle(textDecoration='bold',fgFill='#D9E5F2',wrapText=TRUE),rows=1,cols=1:5,gridExpand=TRUE)
 openxlsx::setColWidths(wb,site,cols=1:5,widths=c(29,33,30,30,30))
 openxlsx::freezePane(wb,site,firstActiveRow=2,firstActiveCol=3)
}
notes <- c(
 'Supplementary Table S10',
 'Admission characteristics of neonates excluded because no matched discharge outcome was available compared with matched and final analysis cohorts',
 'Neotree KCH and SMCH admissions from 1 January 2022 to 31 December 2025.',
 'Excluded cohort: 3,405 master admissions absent from the matched export. All are labelled unmatched. Linked using facility and admission uid; keys unique in both exports.',
 'Matched cohort: 31,223 matched admissions before birthweight, gestation and classification exclusions. Final cohort: 30,432 admissions using the corrected manuscript classifier.',
 'Excluded and final cohorts are disjoint. Matched and final columns overlap; final is a subset of matched. They are not three independent groups.',
 'Categorical percentages use the available nonmissing denominator, explicitly shown as n / available. Missingness percentages use all admissions in the corresponding group.',
 'Continuous summaries are median [25th, 75th percentiles], using available valid values. No hypothesis tests of equivalence were performed.',
 'Birthweight: canonical birthweight_g, valid 300 to 7000 g. GA: 24+0 to 42+0 with day endings 0 to 6. GA shown in genuine weeks using weeks plus days divided by 7; 33.4 in the source denotes 33 weeks 4 days.',
 'Missing/invalid GA here is assessed independently of birthweight. Its counts therefore differ from the sequential flowchart exclusions.',
 'Temperature: 30 to 45 degrees C; hypothermia below 36.5 degrees C. Apgar 0 to 10; low Apgar 0 to 6. LBW below 2500 g; preterm below 37+0.',
 'Sex, inborn/outborn and availability rules follow Table_02_NNU_characteristics.R. Ultrasound includes USS and EUSS; unknown source codes are missing.',
 'No missing outcomes were imputed. Comparisons describe recorded admission characteristics and cannot establish whether mortality was similar or exclude selection bias.',
 'Proposed Results sentence: Compared with neonates included in the analysis, those excluded because of unmatched outcome data had similar sex distributions and median admission temperatures, but were more often outborn (35.2% versus 31.4%) and preterm (41.4% versus 38.5%; Supplementary Table S10).',
 'Suggested limitation: Outcomes were unavailable for 9.8% of admissions. Excluded admissions were more often outborn and preterm, so selection bias cannot be excluded.'
)
openxlsx::addWorksheet(wb,'Definitions and wording')
openxlsx::writeData(wb,'Definitions and wording',data.frame(Note=notes))
openxlsx::setColWidths(wb,'Definitions and wording',1,110)
openxlsx::addStyle(wb,'Definitions and wording',openxlsx::createStyle(wrapText=TRUE,valign='top'),rows=2:(length(notes)+1),cols=1,gridExpand=TRUE)
openxlsx::setRowHeights(wb,'Definitions and wording',rows=2:(length(notes)+1),heights=55)
openxlsx::addWorksheet(wb,'Reproducible R code')
openxlsx::writeData(wb,'Reproducible R code',data.frame(R_code=readLines(SCRIPT_FILE)))
openxlsx::setColWidths(wb,'Reproducible R code',1,150)
openxlsx::saveWorkbook(wb,file.path(out,'Supplementary_Table_outcome_exclusions.xlsx'),overwrite=TRUE)
