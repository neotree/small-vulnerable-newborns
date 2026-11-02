# One ANC definition for recovery, descriptive tables and regression models.
# Visit counts remain missing for labels giving only a range.
anc_count <- function(x) {
    n <- suppressWarnings(as.numeric(trimws(as.character(x))))
    n[!is.finite(n) | n < 0 | n != floor(n)] <- NA_real_
    n
}
anc_band_from_count <- function(x) {
    n <- anc_count(x)
    ifelse(is.na(n), NA_character_, ifelse(n==0, 'none', ifelse(n<=3,'1_3','ge4')))
}
anc_band <- function(count, recovered=NULL) {
    band <- anc_band_from_count(count)
    if (!is.null(recovered)) {
        recovered <- as.character(recovered)
        if (length(recovered)!=length(band)) stop('ANC category length mismatch.')
        bad <- !is.na(recovered) & nzchar(recovered) & !recovered %in% c('none','1_3','ge4')
        if (any(bad)) stop('Unknown recovered ANC category.')
        conflict <- !is.na(band) & !is.na(recovered) & nzchar(recovered) & band!=recovered
        if (any(conflict)) stop('Exact ANC count disagrees with recovered category.')
        fill <- is.na(band) & recovered %in% c('none','1_3','ge4')
        band[fill] <- recovered[fill]
    }
    band
}
anc_factor <- function(count, recovered=NULL, short=FALSE) {
    band <- anc_band(count,recovered)
    labels <- if (short) c(ge4='Adequate (>=4)',none='None (0)',`1_3`='Inadequate (1-3)')
              else c(ge4='Adequate (>=4 visits)',none='None (0 visits)',`1_3`='Inadequate (1-3 visits)')
    x <- unname(labels[band])
    levels <- if (short) unname(labels[c('ge4','1_3','none')]) else unname(labels)
    factor(x,levels=levels)
}
recover_raw_anc <- function(value,label) {
    value_count <- anc_count(value)
    label <- trimws(as.character(label));label[is.na(label)] <- ''
    lower <- tolower(label)
    label_count <- rep(NA_real_,length(label))
    label_count[lower=='none'] <- 0
    exact <- grepl('^[0-9]+[[:space:]]+anc[[:space:]]+visits?$',lower)
    label_count[exact] <- as.numeric(sub('[[:space:]].*$','',lower[exact]))
    numeric_label <- grepl('^[0-9]+([.]0+)?$',lower)
    label_count[numeric_label] <- anc_count(lower[numeric_label])
    conflict <- !is.na(value_count) & !is.na(label_count) & value_count!=label_count
    if (any(conflict)) stop('Raw ANC Value and exact Label conflict; review before applying.')
    n <- value_count
    use_label <- is.na(n) & !is.na(label_count)
    n[use_label] <- label_count[use_label]
    band <- anc_band_from_count(n)
    range_label <- grepl('^more than (4|8)[[:space:]]+anc[[:space:]]+visits$',lower)
    only_band <- is.na(n) & range_label
    band[only_band] <- 'ge4'
    method <- ifelse(!is.na(value_count),'raw_numeric_value',
               ifelse(use_label,'raw_exact_label',ifelse(only_band,'raw_range_label','unresolved')))
    data.frame(count=n,category=band,source=method,label=label,stringsAsFactors=FALSE)
}
