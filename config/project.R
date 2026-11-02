# Base-R project configuration and immutable input checks. No analysis is run here.
resolve_dir <- function(value, root, fallback) {
    if (!nzchar(value)) return(file.path(root,fallback))
    value <- path.expand(value)
    if (!grepl('^(/|[A-Za-z]:[/\\\\])',value)) value <- file.path(root,value)
    value
}
input_paths <- function(root) {
    locked <- read.csv(file.path(root,'config','locked_inputs.csv'),stringsAsFactors=FALSE)
    file.path(resolve_dir(Sys.getenv('LBW_INPUT_DIR'),root,'input'),locked$filename)
}
preflight <- function(root) {
    paths <- input_paths(root)
    missing <- paths[!file.exists(paths)]
    if (length(missing)) stop('Missing locked inputs:\n',paste(basename(missing),collapse='\n'),
        '\nPlace the six authorised files in input/ or set LBW_INPUT_DIR. See input/README.md.',call.=FALSE)
    expected <- read.csv(file.path(root,'config','locked_inputs.csv'),stringsAsFactors=FALSE)
    wrong <- expected$filename[unname(tools::md5sum(paths)) != expected$md5]
    if (length(wrong)) stop('Input lock mismatch:\n',paste(wrong,collapse='\n'),
        '\nDo not change the lock to bypass this check. Confirm the intended extract.',call.=FALSE)
    ref <- file.path(root,'references','official_intergrowth_p10_birthweight_24_42.csv')
    reference_lock <- read.csv(file.path(root,'config','reference_lock.csv'),stringsAsFactors=FALSE)
    if (!file.exists(ref) || unname(tools::md5sum(ref)) != reference_lock$md5[1]) stop('INTERGROWTH reference lock mismatch.',call.=FALSE)
    pkgs <- read.csv(file.path(root,'config','packages.csv'),stringsAsFactors=FALSE)$Package
    missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
    if (length(missing)) stop('Missing packages: ',paste(missing,collapse=', '),
        '\nInstall them in your approved R library (see README.md).',call.=FALSE)
    cat('Preflight passed: six locked inputs, reference table and required packages.\n')
}
