# bass_cost_context.R — R BASS 1.3.1's time per sweep at three points, as context for the C# cost qualification (chip
# brief issues/regression/bass/performance-cost-brief.md, "R"). Context, not a comparison: R's sweep is not the same
# work (its fixed 1/3 move probabilities, F-O6, and its own evaluation order).
#
# The points, on the data BassCostQualificationTests.Counts_EveryPoint exports (BASS_COST_OUT/r):
#   Friedman at n = 500, p = 10, R's defaults;
#   the rich surface at n = 5,000, p = 10, maxBasis = 20, maxInt = 2 (the C# base point's interaction cap);
#   the direct functional at n_x = 100, m = 50, maxBasis = 20, maxInt = 2.
# Each runs R's default lengths (nmcmc = 10,000, nburn = 9,000) once; µs per sweep is the elapsed time over nmcmc.
# R BASS runs as written: its truncated gamma draws come from truncdist, which the renv pins.
# sessionInfo() and the BLAS and LAPACK in use go beside the table.
#
# Usage (cwd = r/):  Rscript scripts/bass_cost_context.R <data folder> <output folder> [nmcmc]

source("oracles/_common.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: Rscript scripts/bass_cost_context.R <data folder> <output folder> [nmcmc]")
data_dir <- args[1]
out_dir <- args[2]
nmcmc <- if (length(args) >= 3) as.integer(args[3]) else 10000L
nburn <- as.integer(round(0.9 * nmcmc))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

bass <- bass_source(c("bass.R", "basis_funcs.R", "mcmc_funcs.R", "mcmc_update.R", "misc.R", "rjmcmc_all.R", "rjmcmc_cat.R",
                      "rjmcmc_cat_func.R", "rjmcmc_des.R", "rjmcmc_des_cat.R", "rjmcmc_des_cat_func.R", "rjmcmc_des_func.R",
                      "rjmcmc_des_ridge.R", "generics.R", "predict.R"))

scalar <- function(file) {
  d <- read.csv(file.path(data_dir, file))
  list(x = as.matrix(d[, setdiff(names(d), "y")]), y = d$y)
}

timed <- function(name, call) {
  set.seed(1)
  t <- system.time(fit <- eval(call))[["elapsed"]]
  data.frame(point = name, nmcmc = nmcmc, nburn = nburn, seconds = t, us_per_sweep = 1e6 * t / nmcmc,
             mean_M = tryCatch(mean(fit$nbasis), error = function(e) NA), stringsAsFactors = FALSE)
}

f <- scalar("friedman.csv")
r <- scalar("rich.csv")
fx <- as.matrix(read.csv(file.path(data_dir, "functional-x.csv"))[, -6])
fy <- as.matrix(read.csv(file.path(data_dir, "functional-y.csv")))
ft <- read.csv(file.path(data_dir, "functional-t.csv"))$t

rows <- rbind(
  timed("Friedman n=500 p=10, R defaults", quote(bass$bass(f$x, f$y, nmcmc = nmcmc, nburn = nburn, verbose = FALSE))),
  timed("rich n=5000 p=10 maxBasis=20 maxInt=2",
        quote(bass$bass(r$x, r$y, maxBasis = 20, maxInt = 2, nmcmc = nmcmc, nburn = nburn, verbose = FALSE))),
  timed("functional n_x=100 m=50 p=5 maxBasis=20 maxInt=2",
        quote(bass$bass(fx, fy, xx.func = ft, maxBasis = 20, maxInt = 2, nmcmc = nmcmc, nburn = nburn, verbose = FALSE)))
)
write.csv(rows, file.path(out_dir, "r-context.csv"), row.names = FALSE)
print(rows)
writeLines(c(capture.output(sessionInfo()), "", "extSoftVersion():", capture.output(print(extSoftVersion())), "",
             paste("La_library():", La_library()), paste("BLAS (sessionInfo):", sessionInfo()$BLAS)),
           file.path(out_dir, "r-session.txt"))
