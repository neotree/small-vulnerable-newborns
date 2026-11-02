# Explicit optional installation; the analysis runner never installs packages.
a <- commandArgs(FALSE);f <- a[grepl('^--file=',a)]
root <- if (length(f)) dirname(normalizePath(sub('^--file=','',f[1]))) else getwd()
pkgs <- read.csv(file.path(root,'config','packages.csv'),stringsAsFactors=FALSE)$Package
missing <- pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if (length(missing)) install.packages(missing,repos='https://cloud.r-project.org')
if (any(!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE))) stop('Packages remain unavailable.')
cat('All required packages are available. Tested versions are in config/packages.csv.\n')
