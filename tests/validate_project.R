# Data-free checks used locally and in GitHub Actions.
if (!exists('ROOT')) {
 a<-commandArgs(FALSE); f<-a[grepl('^--file=',a)]
 ROOT<-if(length(f)) dirname(dirname(normalizePath(sub('^--file=','',f[1])))) else normalizePath('.')
}
spec<-read.csv(file.path(ROOT,'config','scripts.csv'),stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(spec$script),all(file.exists(file.path(ROOT,'scripts',spec$script))),
          is.logical(spec$uss_pass),!anyNA(spec$uss_pass))
files<-list.files(ROOT,pattern='[.]R$',recursive=TRUE,full.names=TRUE)
for(f in files) parse(f)
ref<-read.csv(file.path(ROOT,'references','official_intergrowth_p10_birthweight_24_42.csv'))
stopifnot(nrow(ref)==266L,!anyDuplicated(ref[c('sex','week','day')]),all(ref$sex %in% c('B','G')),
          all(ref$day %in% 0:6),all(ref$p10_g>0))
locks<-read.csv(file.path(ROOT,'config','locked_inputs.csv'),stringsAsFactors=FALSE)
stopifnot(nrow(locks)==6L,all(grepl('^[a-f0-9]{32}$',locks$md5)))
for(f in files) {
 if (basename(f)=='validate_project.R') next
 s<-readLines(f,warn=FALSE)
 stopifnot(!any(grepl('/Users/|/Volumes/|ClaudeCowork|Local-intergrowth-final',s)))
}
cat(sprintf('Validated %d R files, %d analysis scripts and 266 reference rows without patient data.\n',length(files),nrow(spec)))
