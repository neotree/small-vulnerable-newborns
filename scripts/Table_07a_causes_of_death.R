# ==============================================================================
# Purpose: Table 7a -- leading causes of neonatal death by SVN phenotype,
# reported as n and % of deaths in each group, both sites combined and by
# site.
#
# Definitions: population = svn_filtered (valid outcome); death totals:
# All-SVN 4,235, Term-AGA 1,109, total 5,344.
#
# Output: table_07a_causes_by_group.csv, table_07a_causes_by_group_and_site.csv,
#         table_07a_causes_top_wide.csv
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
suppressPackageStartupMessages({library(dplyr); library(stringr); library(tidyr)})
cat("\n=== 45_causes_of_death_by_svn.R ===\n")

normalise_cause <- function(code){
    code_up <- toupper(trimws(code))
    dplyr::case_when(
        code_up %in% c("PRRDS","PREMRD","EXPREM","VPREM","PR","RDS","EXLBW","VLBW","LBW",
                       "TERMRD","RDNT","TTN","CLEFTRD") ~ "Prematurity/RDS",
        code_up %in% c("BA","HIE","SHIE","MAS","MA","ASP","LAS","AP","BBA") ~ "Birth Asphyxia/HIE",
        code_up %in% c("NS","SEPS","HYPOGSY") ~ "Neonatal Sepsis",
        code_up == "EONS" ~ "Neonatal Sepsis (Early Onset)",
        code_up == "LONS" ~ "Neonatal Sepsis (Late Onset)",
        code_up %in% c("G","GSCH","GSCHIS","GASTRO") ~ "Gastroschisis",
        code_up %in% c("OMPH","OM","EXOM") ~ "Omphalocele",
        code_up %in% c("CONG","CA","OCA") ~ "Congenital Anomaly",
        code_up %in% c("OTH") ~ "Cause Not Specified",
        TRUE ~ code)
}

nnu <- process_nnu_data(load_nnu_data())
PH <- c("Term-AGA","Term-SGA","Preterm-AGA","Preterm-SGA")
d <- nnu$svn_filtered %>%
    filter(!is.na(outcome_bin), svn_category %in% PH) %>%
    mutate(svn_group = ifelse(as.character(svn_category)=="Term-AGA","Term-AGA","All-SVN"),
           facility  = as.character(facility))
if (!"causedeath" %in% names(d)) stop("Column 'causedeath' not found in svn_filtered.")

deaths <- d %>% filter(outcome_bin == 1) %>%
    mutate(cod = str_trim(as.character(causedeath)),
           cod = ifelse(cod=="" | is.na(cod) | cod=="NA","Unknown",cod),
           cause = normalise_cause(cod))

cause_table <- function(df,label){
    tot <- nrow(df)
    df %>% count(cause,name="n") %>% arrange(desc(n)) %>%
        mutate(group=label, total_deaths=tot, pct=round(100*n/tot,1),
               rank=row_number(), n_pct=sprintf("%d (%.1f%%)",n,pct)) %>%
        select(group,rank,cause,n,pct,n_pct,total_deaths)
}

by_group <- bind_rows(
    cause_table(filter(deaths,svn_group=="All-SVN"),"All-SVN"),
    cause_table(filter(deaths,svn_group=="Term-AGA"),"Term-AGA"),
    cause_table(deaths,"Overall"))
write.csv(by_group, file.path(OUT2,"table_07a_causes_by_group.csv"), row.names=FALSE)

by_group_site <- deaths %>%
    group_split(facility, svn_grp = svn_group) %>%
    lapply(function(g) cause_table(g, paste0(unique(g$facility),"-",unique(g$svn_group)))) %>%
    bind_rows()
# add overall-per-site rows
by_site_overall <- deaths %>% group_split(facility) %>%
    lapply(function(g) cause_table(g, paste0(unique(g$facility),"-Overall"))) %>% bind_rows()
write.csv(bind_rows(by_group_site, by_site_overall),
          file.path(OUT2,"table_07a_causes_by_group_and_site.csv"), row.names=FALSE)

# wide top-8 comparison (both sites combined)
topcauses <- by_group %>% filter(group=="Overall") %>% slice_head(n=8) %>% pull(cause)
wide <- by_group %>%
    filter(cause %in% topcauses) %>%
    mutate(cause=factor(cause,levels=topcauses)) %>%
    select(group,cause,n_pct) %>%
    pivot_wider(names_from=group, values_from=n_pct, values_fill="0 (0.0%)") %>%
    arrange(cause)
write.csv(wide, file.path(OUT2,"table_07a_causes_top_wide.csv"), row.names=FALSE)

cat(sprintf("Deaths: All-SVN %d, Term-AGA %d, Overall %d\n",
            sum(deaths$svn_group=="All-SVN"), sum(deaths$svn_group=="Term-AGA"), nrow(deaths)))
cat("Top causes (All-SVN vs Term-AGA vs Overall):\n"); print(as.data.frame(wide))
cat("\nSaved table_07a_causes_by_group.csv, table_07a_causes_by_group_and_site.csv, table_07a_causes_top_wide.csv\n")
cat("45_causes_of_death_by_svn.R complete.\n")
