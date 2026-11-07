# khaos_readouts_oracle.R — KHAOS's posterior readouts for one coefficient prior on one data set (E19's build brief, slice 0). Never
# traces: khaos's chain and the engine's differ by design (lambda drawn rather than integrated, g0^2 by a Laplace-proposal MH rather
# than a slice step, the corrections slice 0 records), so the engine meets these readouts from slice 1 on, each difference stated beside
# the comparison.
#
# Usage (cwd = r/):  Rscript oracles/khaos_readouts_oracle.R <data.csv> <out.json> <prior> <grid.csv> <nmcmc> <nburn> <seed> [key=value ...]
#   data: headerless CSV, inputs in [0, 1] then y in the last column. grid: headerless CSV of inputs, the surface's points.
#   prior: ridge or gprior (adaptive_khaos's prior_type); key=value pairs go to it as numeric arguments (degree, order, zeta, tau2, ...).
#
# Emits { prior, args, N, p, nmcmc, nburn, seed, seconds, surface: {mean, q05, q50, q95} at the grid (the draws' f, predict as
#   written), M, s2, lambda, g0sq (gprior): {mean, sd, q05, q50, q95}, count_accept, count_propose (birth, death, degree mutation,
#   variable swap), sobol: sobol_khaos's own readout, its indices over the explained variance plus s2 (labels, S_mean, T_mean, T_q05,
#   T_q95, leftover_mean) and T_f_mean, the total indices over the explained variance alone (the variance of f, the engine's
#   convention), terms: {count, repeated (terms that repeat a variable), mean_order, mean_degree} over the kept draws }.

source("oracles/khaos_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) stop("usage: khaos_readouts_oracle.R <data.csv> <out.json> <prior> <grid.csv> <nmcmc> <nburn> <seed> [key=value ...]")
D <- read_matrix(args[1]); prior <- args[3]; G <- read_matrix(args[4])
nmcmc <- as.integer(args[5]); nburn <- as.integer(args[6]); seed <- as.integer(args[7])
extra <- if (length(args) > 7) {
  kv <- strsplit(args[-(1:7)], "=", fixed = TRUE)
  setNames(lapply(kv, function(p) as.numeric(p[2])), vapply(kv, `[`, "", 1))
} else list()
if (!prior %in% c("ridge", "gprior")) stop("unknown prior ", prior)
X <- D[, -ncol(D), drop = FALSE]; y <- D[, ncol(D)]

set.seed(seed)
t0 <- proc.time()[["elapsed"]]
fit <- do.call(adaptive_khaos, c(list(X = X, y = y, prior_type = prior, nmcmc = nmcmc, nburn = nburn, verbose = FALSE), extra))
secs <- proc.time()[["elapsed"]] - t0

summ <- function(x) if (is.null(x) || all(is.na(x))) NULL else
  list(mean = mean(x), sd = sd(x), q05 = unname(quantile(x, 0.05)), q50 = unname(quantile(x, 0.5)), q95 = unname(quantile(x, 0.95)))
f_grid <- predict(fit, newdata = G)
pick <- function(q) unname(apply(f_grid, 2, quantile, probs = q))

sob <- sobol_khaos(fit, plot_it = FALSE)
Tm <- as.matrix(sob$T); explained <- 1 - sob$leftover
held <- if (length(fit$nint)) col(fit$nint) <= fit$nbasis else matrix(FALSE, 0, 0)
orders <- fit$nint[held]; degrees <- fit$dtot[held]
terms <- sum(held)
repeated <- sum(mapply(function(i, m) anyDuplicated(fit$vars[i, m, seq_len(fit$nint[i, m])]) > 0,
                        row(fit$nint)[held], col(fit$nint)[held]))

emit_json(list(prior = prior, args = extra, N = length(y), p = ncol(X), nmcmc = nmcmc, nburn = nburn, seed = seed, seconds = secs,
               surface = list(mean = colMeans(f_grid), q05 = pick(0.05), q50 = pick(0.5), q95 = pick(0.95)),
               M = summ(fit$nbasis), s2 = summ(fit$s2), lambda = summ(fit$lam), g0sq = summ(fit$g2),
               count_accept = unname(fit$count_accept), count_propose = unname(fit$count_propose),
               sobol = list(labels = names(sob$S), S_mean = unname(colMeans(sob$S, na.rm = TRUE)), T_mean = unname(colMeans(Tm, na.rm = TRUE)),
                            T_q05 = unname(apply(Tm, 2, quantile, 0.05, na.rm = TRUE)), T_q95 = unname(apply(Tm, 2, quantile, 0.95, na.rm = TRUE)),
                            leftover_mean = mean(sob$leftover, na.rm = TRUE), T_f_mean = unname(colMeans(Tm / explained, na.rm = TRUE))),
               terms = list(count = terms, repeated = repeated, mean_order = mean(orders), mean_degree = mean(degrees))),
          args[2], null = "null", na = "null")
