# ==============================================================================
# Purpose: adds the pooled "All-SVN" column to Table 9 -- univariable
# admission-temperature odds ratios for neonatal death in the pooled All-SVN
# group (Term-SGA + Preterm-AGA + Preterm-SGA; Term-AGA excluded).
#
# Definitions: pools the per-phenotype temperature x death counts from
# dep_table09_temp_cat_sample_sizes.csv and computes the 2x2 odds ratio vs
# Normal (36.5-37.5C) directly -- for a single categorical predictor this
# equals the glm(death ~ temp_cat) MLE.
#
# Output: helper_table09_temp_or_allsvn.csv
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
IN  <- file.path(OUTPUT_DIR, "dep_table09_temp_cat_sample_sizes.csv")
OUT <- file.path(OUTPUT_DIR, "helper_table09_temp_or_allsvn.csv")
if (!file.exists(IN)) stop("Missing dep_table09_temp_cat_sample_sizes.csv - run dep_Table09_temperature_categorical_OR.R first")

x <- read.csv(IN, stringsAsFactors = FALSE)
# pooled All-SVN = the three vulnerable phenotypes
svn <- subset(x, svn_category %in% c("Term-SGA", "Preterm-AGA", "Preterm-SGA"))
agg <- aggregate(cbind(n, n_deaths) ~ temp_cat, data = svn, FUN = sum)

ref <- agg[agg$temp_cat == "4_Normal", ]
c_d <- ref$n_deaths; c_a <- ref$n - ref$n_deaths          # normal deaths / alive
lab <- c("1_SevereHypothermia" = "Severe hypothermia (<32C)",
         "2_Hypothermia"       = "Moderate hypothermia (32-35.9C)",
         "3_LowNormal"         = "Mild hypothermia (36.0-36.4C)",
         "5_Hyperthermia"      = "Hyperthermia (>37.5C)")

res <- do.call(rbind, lapply(names(lab), function(k) {
    r <- agg[agg$temp_cat == k, ]
    a_d <- r$n_deaths; a_a <- r$n - r$n_deaths            # category deaths / alive
    or  <- (a_d / a_a) / (c_d / c_a)
    se  <- sqrt(1/a_d + 1/a_a + 1/c_d + 1/c_a)            # SE of log(OR)
    lo  <- exp(log(or) - 1.96 * se); hi <- exp(log(or) + 1.96 * se)
    p   <- 2 * pnorm(-abs(log(or) / se))
    data.frame(category = lab[[k]], n = r$n, deaths = a_d,
               OR = round(or, 2), CI = sprintf("%.2f-%.2f", lo, hi),
               p = signif(p, 3),
               cell = sprintf("%.2f%s", or, ifelse(p < 0.05, "*", "")))
}))
write.csv(res, OUT, row.names = FALSE)
cat("\nAll-SVN admission-temperature ORs (pooled Term-SGA + Preterm-AGA + Preterm-SGA)\n")
cat(sprintf("Reference Normal: n=%d, deaths=%d\n", ref$n, c_d))
print(res, row.names = FALSE)
cat("\nSaved", OUT, "\n56_temp_or_allsvn.R complete.\n")
