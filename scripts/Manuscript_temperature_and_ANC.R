# Continuous temperature estimates for Table 4 and Results; ANC availability audit.
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
nnu <- process_nnu_data(load_nnu_data())
d <- nnu$svn_filtered %>%
  filter(!is.na(outcome_bin), svn_category %in% c('Term-AGA','Term-SGA','Preterm-AGA','Preterm-SGA')) %>%
  mutate(temp_val=suppressWarnings(as.numeric(temperature)),
         temp_val=ifelse(temp_val>=30 & temp_val<=45,temp_val,NA_real_),
         group=factor(ifelse(svn_category=='Term-AGA','Term-AGA','All-SVN'),levels=c('Term-AGA','All-SVN')),
         anc_num=anc_count(antenatalcare), anc_category=anc_band(antenatalcare,antenatalcare_category))
out <- OUTPUT_DIR
fits <- lapply(c('All neonates','All-SVN','Term-AGA'), function(g) {
  x <- if(g=='All neonates') d else filter(d,as.character(group)==g)
  f <- glm(outcome_bin ~ temp_val, data=x,family=binomial())
  b <- coef(summary(f))['temp_val',]
  ci <- suppressMessages(confint(f,'temp_val'))
  data.frame(population=g,n=nobs(f),OR_decrease=exp(-b[1]),lo_decrease=exp(-ci[2]),hi_decrease=exp(-ci[1]),p=b[4])
})
result <- bind_rows(fits)
f0 <- glm(outcome_bin ~ temp_val + group,data=d,family=binomial())
f1 <- glm(outcome_bin ~ temp_val * group,data=d,family=binomial())
a <- anova(f0,f1,test='LRT')
result$interaction_AllSVN_TermAGA <- a$`Pr(>Chi)`[2]
write.csv(result,file.path(out,'manuscript_temperature_population_check.csv'),row.names=FALSE)
anc <- d %>% group_by(facility) %>% summarise(n=n(),
 exact_recorded=sum(!is.na(anc_num)),exact_missing=sum(is.na(anc_num)),exact_missing_pct=100*mean(is.na(anc_num)),
 category_recorded=sum(!is.na(anc_category)),category_missing=sum(is.na(anc_category)),category_missing_pct=100*mean(is.na(anc_category)),
 category_only=sum(is.na(anc_num)&!is.na(anc_category)),.groups='drop')
anc_year <- d %>% mutate(year=substr(adm_date_parsed,1,4)) %>% group_by(facility,year) %>% summarise(n=n(),
 exact_missing=sum(is.na(anc_num)),category_missing=sum(is.na(anc_category)),category_missing_pct=100*mean(is.na(anc_category)),
 category_only=sum(is.na(anc_num)&!is.na(anc_category)),.groups='drop')
write.csv(anc_year,file.path(out,'manuscript_anc_missingness_by_year.csv'),row.names=FALSE)
write.csv(anc,file.path(out,'manuscript_anc_missingness_check.csv'),row.names=FALSE)
cat('\nTEMPERATURE CHECK\n');print(result)
cat('\nANC CHECK\n');print(anc)
