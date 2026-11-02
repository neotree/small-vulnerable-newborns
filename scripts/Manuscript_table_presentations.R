# R equivalents of calculations formerly made while formatting the Word tables.
# Dependencies: Table_03_04_distribution_and_CFR.R.
SCRIPT_DIR<-local({
 a<-commandArgs(FALSE);f<-a[grepl('^--file=',a)]
 src<-Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
 if(length(src))dirname(normalizePath(tail(src,1)[[1]]))
 else if(length(f))dirname(normalizePath(sub('^--file=','',f[1]))) else normalizePath('scripts')
})
source(file.path(SCRIPT_DIR,'00a_setup_shared.R'))
results<-list()
for(which in c('svn','bw')) {
 x<-read.csv(file.path(OUTPUT_DIR,paste0('table_03_cfr_by_',which,'.csv')))
 x<-x[x$site=='Combined',]
 reference<-if(which=='svn')'Term-AGA' else 'NBW'
 ref<-x[x$category==reference,]
 stopifnot(nrow(ref)==1)
 x$OR<-x$ci_lower<-x$ci_upper<-x$p<-NA_real_
 for(i in seq_len(nrow(x))) {
  if(x$category[i]==reference){x$OR[i]<-1;next}
  a<-x$deaths[i];b<-x$n[i]-a;c<-ref$deaths;d<-ref$n-c
  if(any(c(a,b,c,d)<=0))next
  odds<-a*d/(b*c);se<-sqrt(1/a+1/b+1/c+1/d)
  x$OR[i]<-odds;x$ci_lower[i]<-exp(log(odds)-1.96*se)
  x$ci_upper[i]<-exp(log(odds)+1.96*se);x$p[i]<-2*pnorm(-abs(log(odds)/se))
 }
 x$classification<-which;results[[which]]<-x
}
write.csv(bind_rows(results),file.path(OUTPUT_DIR,'manuscript_table_3_combined_odds.csv'),row.names=FALSE)
cross<-read.csv(file.path(OUTPUT_DIR,'table_04_distribution_crosstab.csv'),check.names=FALSE)
svn<-results$svn
denom<-sum(svn$n[svn$category %in% SVN_LEVELS])
long<-tidyr::pivot_longer(cross,cols=where(is.numeric),names_to='phenotype',values_to='n')
long$pct_all_admissions<-100*long$n/denom
write.csv(long,file.path(OUTPUT_DIR,'manuscript_supplementary_s1_counts_percentages.csv'),row.names=FALSE)
