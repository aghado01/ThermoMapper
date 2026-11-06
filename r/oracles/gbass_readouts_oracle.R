# gbass_readouts_oracle.R — GBASS's posterior readouts for one wrapper on one data set (E1's build brief, slice 0). Never traces:
# GBASS's chain and the engine's differ by design (continuous knots, mutation of a whole term, pseudo-counts 1, a ridge with the
# intercept inside), so the engine meets these readouts with each difference stated beside the comparison.
#
# Usage (cwd = r/):  Rscript oracles/gbass_readouts_oracle.R <data.csv> <out.json> <model> <grid.csv> <nmcmc> <nburn> <seed> [key=value ...]
#   data: headerless CSV, inputs in [0, 1] then y in the last column. grid: headerless CSV of inputs, the surface's points.
#   model: tbass, qbass, nwbass or hbass; key=value pairs go to the wrapper as numeric arguments (df, q, s_beta, m_gamma, s_gamma, ...).
#
# Emits { model, args, N, seconds, surface: {mean, q05, q50, q95} at the grid (the draws' f, predict(predictive = FALSE)),
#   shifted_mean: the mean of f + sqrt(w) beta E[v] at the grid (predict(bias_correct = TRUE), qbass and nwbass), w, s2 (the noise variance
#   w(c E[v] + beta^2 Var[v]) per draw, GBASS's own), beta, gamma, tau,
#   M: {mean, sd, q05, q50, q95}, below_fit: the share of y below the posterior mean surface at the data, cnt1, cnt2 (move counts
#   proposed and accepted: birth, death, mutation) }.

source("oracles/_common.R")
suppressPackageStartupMessages(library(GBASS))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("usage: gbass_readouts_oracle.R <data.csv> <out.json> <model> <grid.csv> <nmcmc> <nburn> <seed> [key=value ...]")
D <- read_matrix(args[1]); model <- args[3]; G <- read_matrix(args[4])
nmcmc <- as.integer(args[5]); nburn <- as.integer(args[6]); seed <- as.integer(args[7])
extra <- if (length(args) > 7) {
  kv <- strsplit(args[-(1:7)], "=", fixed = TRUE)
  setNames(lapply(kv, function(p) as.numeric(p[2])), vapply(kv, `[`, "", 1))
} else list()
if (!model %in% c("tbass", "qbass", "nwbass", "hbass")) stop("unknown model ", model)
X <- D[, -ncol(D), drop = FALSE]; y <- D[, ncol(D)]

set.seed(seed)
t0 <- proc.time()[["elapsed"]]
fit <- do.call(model, c(list(X = X, y = y, nmcmc = nmcmc, nburn = nburn, verbose = FALSE), extra))
secs <- proc.time()[["elapsed"]] - t0

summ <- function(x) if (is.null(x) || all(is.na(x))) NULL else
  list(mean = mean(x), sd = sd(x), q05 = unname(quantile(x, 0.05)), q50 = unname(quantile(x, 0.5)), q95 = unname(quantile(x, 0.95)))
f_grid <- predict(fit, newdata = G, predictive = FALSE)
f_data <- predict(fit, newdata = X, predictive = FALSE)
# GBASS warns that the shifted mean is seldom wanted under qbass; here it is a readout, asked for on purpose.
shifted <- if (model %in% c("qbass", "nwbass")) colMeans(suppressWarnings(predict(fit, newdata = G, predictive = FALSE, bias_correct = TRUE))) else NULL
pick <- function(q) unname(apply(f_grid, 2, quantile, probs = q))

emit_json(list(model = model, args = extra, N = length(y), nmcmc = nmcmc, nburn = nburn, seed = seed, seconds = secs,
               surface = list(mean = colMeans(f_grid), q05 = pick(0.05), q50 = pick(0.5), q95 = pick(0.95)),
               shifted_mean = shifted, w = summ(fit$w), s2 = summ(fit$s2), beta = summ(fit$beta), gamma = summ(fit$gamma), tau = summ(fit$tau),
               M = summ(fit$M), below_fit = mean(y < colMeans(f_data)), cnt1 = fit$cnt1, cnt2 = fit$cnt2),
          args[2], null = "null", na = "null")
