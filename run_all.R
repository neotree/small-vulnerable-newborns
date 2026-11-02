# Run the locked-data manuscript pipeline on macOS or Linux/DSH.
# Rscript run_all.R --check     : validate inputs, references and packages only
# Rscript run_all.R --validate  : parse/structural checks without patient data
# Rscript run_all.R            : validate, then run each script in a fresh R process,
#                                once for all admissions (output/all) and once for the
#                                ultrasound-dated subset (output/uss) where uss_pass is TRUE
a <- commandArgs(FALSE)
f <- a[grepl('^--file=', a)]
src <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
project_file <- if (length(src)) tail(src, 1)[[1]] else if (length(f)) sub('^--file=', '', f[1]) else 'run_all.R'
ROOT <- dirname(normalizePath(project_file, mustWork=TRUE))
source(file.path(ROOT, 'config', 'project.R'))
args <- commandArgs(trailingOnly=TRUE)
if (any(!args %in% c('--check','--validate'))) stop('Supported options: --check, --validate')
if ('--validate' %in% args) {
    source(file.path(ROOT, 'tests', 'validate_project.R'))
    quit(status=0L)
}
preflight(ROOT)
if ('--check' %in% args) quit(status=0L)
output <- resolve_dir(Sys.getenv('LBW_OUTPUT_DIR'), ROOT, 'output')
dir.create(output, recursive=TRUE, showWarnings=FALSE)
logs <- file.path(output, 'logs')
dir.create(logs, recursive=TRUE, showWarnings=FALSE)
spec <- read.csv(file.path(ROOT,'config','scripts.csv'), stringsAsFactors=FALSE)
Sys.setenv(LBW_INPUT_DIR=resolve_dir(Sys.getenv('LBW_INPUT_DIR'),ROOT,'input'),
           LBW_OUTPUT_DIR=output)
write.csv(data.frame(filename=basename(input_paths(ROOT)),md5=unname(tools::md5sum(input_paths(ROOT)))),
          file.path(logs,'input_checksums.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(logs,'session_info.txt'))
file.copy(file.path(ROOT,'config','scripts.csv'),file.path(logs,'scripts.csv'),overwrite=TRUE)
rs <- file.path(R.home('bin'),'Rscript')
results <- list()
for (pop in c('all','USS')) {
    run <- if (pop == 'USS') spec$script[spec$uss_pass] else spec$script
    pop_logs <- if (pop == 'USS') file.path(logs,'uss') else logs
    dir.create(pop_logs, recursive=TRUE, showWarnings=FALSE)
    Sys.setenv(POPULATION=pop)
    for (i in seq_along(run)) {
        script <- run[i]
        cat(sprintf('[%s %d/%d] %s\n',pop,i,length(run),script))
        start <- Sys.time()
        status <- system2(rs,args=c('--vanilla',shQuote(file.path(ROOT,'scripts',script))),
                          stdout=file.path(pop_logs,paste0(script,'.log')),
                          stderr=file.path(pop_logs,paste0(script,'.log')))
        results[[length(results)+1]] <- data.frame(population=pop,script=script,status=status,
                                                   seconds=as.numeric(difftime(Sys.time(),start,units='secs')))
        write.csv(do.call(rbind,results),file.path(logs,'run_summary.csv'),row.names=FALSE)
        if (status != 0L) stop('Stopped at ',script,' (',pop,'). See its log in ',pop_logs,call.=FALSE)
    }
}
cat(sprintf('Completed %d script runs (all %d, USS %d). Outputs: %s\n',length(results),
            nrow(spec),sum(spec$uss_pass),output))
