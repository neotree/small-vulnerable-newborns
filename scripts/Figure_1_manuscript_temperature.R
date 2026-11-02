# Manuscript Figure 1: temperature odds ratios plus the admission-count panel.
# Dependencies: dep_Table09_temperature_categorical_OR.R and
# Figure_07_temperature_categorical_ORs_by_SVN.R. Uses aggregate CSVs only.
SCRIPT_DIR <- local({
 a<-commandArgs(FALSE);f<-a[grepl('^--file=',a)]
 src<-Filter(Negate(is.null),lapply(sys.frames(),function(x)x$ofile))
 if(length(src))dirname(normalizePath(tail(src,1)[[1]]))
 else if(length(f))dirname(normalizePath(sub('^--file=','',f[1]))) else normalizePath('scripts')
})
source(file.path(SCRIPT_DIR,'00a_setup_shared.R'))
suppressPackageStartupMessages(library(patchwork))
or<-read.csv(file.path(OUTPUT_DIR,'figure_07_temperature_by_svn_table.csv'))
counts<-read.csv(file.path(OUTPUT_DIR,'dep_table09_temp_cat_sample_sizes.csv'))
groups<-c('All-SVN','Term-AGA','Term-SGA','Preterm-AGA','Preterm-SGA')
cats<-c('Severe hypothermia','Hypothermia','Mild hypothermia','Normal (reference)','Hyperthermia')
keys<-c('1_SevereHypothermia','2_Hypothermia','3_LowNormal','4_Normal','5_Hyperthermia')
labels<-c('Severe hypothermia\n<32.0°C','Moderate hypothermia\n32.0–35.9°C',
 'Mild hypothermia\n36.0–36.4°C','Normal (reference)\n36.5–37.5°C','Hyperthermia\n>37.5°C')
cols<-c('All-SVN'='#202020','Term-AGA'='#2166ac','Term-SGA'='#67a9cf','Preterm-AGA'='#d73027','Preterm-SGA'='#fdae61')
or$x<-match(or$temperature_category,cats)
or$group<-factor(or$group,levels=groups)
p<-ggplot(or,aes(x=x,y=OR,colour=group,group=group))+
 geom_hline(yintercept=1,linetype='dashed',colour='grey50')+
 geom_line(linewidth=.8)+geom_point(size=2.4)+
 scale_y_log10(breaks=c(1,3,10,20))+scale_x_continuous(breaks=1:5,labels=labels)+
 scale_colour_manual(values=cols,name=NULL)+
 labs(title='A  Temperature and neonatal mortality by SVN group',x=NULL,
      y='Unadjusted odds ratio (log scale)')+
 theme_minimal(base_size=10)+theme(legend.position='top',panel.grid.minor=element_blank())
ct<-expand.grid(group=groups,temp_cat=keys,stringsAsFactors=FALSE)
ct$n<-mapply(function(g,k){
 ph<-if(g=='All-SVN')groups[3:5] else g
 sum(counts$n[counts$svn_category %in% ph & counts$temp_cat==k])
},ct$group,ct$temp_cat)
ct$x<-match(ct$temp_cat,keys)
ct$y<-length(groups)+1-match(ct$group,groups)
q<-ggplot(ct,aes(x=x,y=y,label=format(n,big.mark=',',trim=TRUE)))+
 geom_text(size=3.2)+scale_x_continuous(breaks=1:5,labels=c('Severe','Moderate','Mild','Normal','Hyperthermia'),position='top')+
 scale_y_continuous(breaks=1:5,labels=rev(groups))+
 labs(title='B  Admissions in each temperature category (n)',x=NULL,y=NULL)+
 theme_minimal(base_size=10)+theme(panel.grid=element_blank(),axis.ticks=element_blank())
fig<-(p/q)+plot_layout(heights=c(3,1.5))+
 plot_annotation(caption=paste('All-SVN pools Term-SGA, Preterm-AGA and Preterm-SGA.',
 'Exact odds ratios and confidence intervals are in Supplementary Table S6.'))
save_figure(fig,'manuscript_figure_1_temperature',width=11.5,height=8.8)
write.csv(ct[c('group','temp_cat','n')],file.path(OUTPUT_DIR,'manuscript_figure_1_temperature_counts.csv'),row.names=FALSE)
