# bass_calibration_defects_oracle.R — defect fixtures for R BASS's calibration (port brief, calibration segment, D09;
# C-O1).
#
# Usage (cwd = r/):  Rscript oracles/bass_calibration_defects_oracle.R <basis.csv> <out.json>
#   basis: headerless CSV, one row per output location: y_m, y_s, then the basis columns K_1, K_2 (m x 4, m even),
#     with y_s not all 1.
#
# R's calibrate.bassBasis (calib.R) runs as written for one iteration (nmcmc = 2) on a mock emulator: a bassBasis
# object whose two components predict w_k(theta) = c_k0 + c_k1 theta_1 + c_k2 theta_2 in physical coordinates through
# R's own predict.bassBasis (bassPCA.R), with 50 emulator draws of component noise and seeded truncation residuals.
# theta has bounds [10, 20] x [-1, 1]; y = y_m + y_s * (K w(theta*)) + noise at theta* = (14, 0.3); two variance
# groups, the first and second halves of the outputs; eight temperatures 1.05^(0:7); no decorrelation step or swap
# runs in one iteration. Sections:
#   stale_likelihood: functional branch with a two-column discrepancy basis. After the discrepancy draw R keeps the
#     cached likelihood; when the theta move is rejected, the returned state's cached llik differs from the
#     likelihood of its own (pred, discrep, cov), recomputed here as -(ldet + r' inv r)/2 with r = y - pred - discrep.
#     Per rung: rejected or not, and cached minus recomputed (the error in that iteration's theta acceptance).
#     Corrected: 0.
#   branch_by_shape: R's branch is whichever shape the discrepancy matrix has: type = 1 (functional) with a square
#     m x m discrepancy basis runs the multivariate branch, and type = 2 (multivariate) with a two-column one runs
#     the functional branch. Corrected: the declared type.
#   coordinates: the returned theta against the inputs the emulator was called with; R's starting theta (its
#     theta[1, , ]) against the emulator's first call, which predicts at it. Corrected: physical coordinates,
#     bounds[1] + theta (bounds[2] - bounds[1]).
#   unit_mismatch: each rung's returned log-determinant against dense log-determinants of the covariance in R's
#     units, diag(v_trunc + s2_g) + K diag(e) K', and in physical units, diag(y_s^2 v_trunc + s2_g) +
#     (y_s K) diag(e) (y_s K)', v_trunc being the stored residuals' variances (type 1), e the emulator draw's
#     component variances. Corrected: the physical one.
#   s2_df_switch: the `if (any(s2.df == 0))` statement evaluated as written: with s2.df = (0, 10) the second group's
#     log prior kernel is -log(1 + sigma^2) instead of its IG(5, 5 sd^2) kernel, -(a + 1) log sigma^2 - b / sigma^2.
#   tempering: R's variance move, discrepancy draw and swap at inverse temperature b, on one variance group (the
#     outputs' noise, IG(5, 0.05)) with theta at theta*, the discrepancy at 0 and the covariance at diag(v_trunc +
#     sigma^2) + K diag(0.01, 0.01) K'. Ladder 1, 2, 5 (b = 1, 0.5, 0.2).
#     - The variance move's `alpha`, evaluated as written, is b times (log likelihood + log prior kernel + log
#       Jacobian): R tempers the Jacobian and the prior with the likelihood. D09's likelihood path tempers the
#       likelihood alone. alpha_linear_gap checks the scaling.
#     - On a grid of 61 values of u = log sigma^2 in [log 1e-4, 0], the variance move as a random walk to the
#       neighbouring values, with R's acceptance and with D09's; each kernel's stationary law by GTH, against R's
#       powered target exp(b (llik + ldig + u)) and the likelihood path's exp(b llik + ldig + u).
#     - R's swap, evaluated as written over every pair of grid states at neighbouring rungs, against the log ratio
#       detailed balance needs for R's own rung laws and for the likelihood path's; D09's swap,
#       (b1 - b2)(llik(x2) - llik(x1)), against its kernels' laws. R's swap leaves out the Jacobian and the
#       discrepancy prior, which its own moves temper.
#     - The discrepancy draw's covariance, R's `discrep.S / itl[t]` evaluated as written, against the conditional
#       covariances of the powered target, S / b, and of the likelihood path, (I / tau + b D' Sigma^-1 D)^-1, with
#       S = (I / tau + D' Sigma^-1 D)^-1 and R's tau = 1; and R's swap between states that differ only in the
#       discrepancy v, against the ratio its powered target needs (the gap is the omitted v prior).
#
#   inputs: what the C# parity tests rebuild the setup from: y, the residuals, the emulator's 50 draws of component
#     variances, the coefficients, the bounds, the groups and the discrepancy basis; run: the functional run's state per
#     rung (theta at the start and the end, normalized, the emulator draw, the variances, the component variances, the
#     discrepancy coefficients, the cached log-likelihood and log-determinant); and the tempering grid's u, llik and log
#     prior kernel. Written at 17 significant digits, so every double round-trips.
#
# Emits { m, theta_star, stale_likelihood, branch_by_shape, coordinates, unit_mismatch, s2_df_switch, tempering, inputs,
#         run }.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_calibration_defects_oracle.R <basis.csv> <out.json>")
Bm <- read_matrix(args[1])
y_m <- Bm[, 1]; y_s <- Bm[, 2]; K <- Bm[, 3:4]
m <- nrow(Bm)
if (m %% 2 != 0) stop("m must be even")
bass <- bass_source(c("misc.R", "bassPCA.R", "calib.R"))

coef <- rbind(c(0.5, 0.08, -0.6), c(-0.3, -0.05, 0.9))
seen <- list()
predict.oracle_linear <- function(object, newdata, mcmc.use, ...) {
  seen[[length(seen) + 1]] <<- as.matrix(newdata)
  w <- cbind(1, as.matrix(newdata)) %*% object$coef
  matrix(rep(w, each = length(mcmc.use)), length(mcmc.use))
}
set.seed(31)
n_train <- 10
emulator <- structure(list(
  mod.list = lapply(1:2, function(k) structure(list(coef = coef[k, ], s2 = 0.01 * rgamma(50, 5, 5)), class = "oracle_linear")),
  dat = list(basis = K, trunc.error = matrix(rnorm(m * n_train, sd = 0.05), m, n_train), n.pc = 2, y.m = y_m, y.s = y_s)),
  class = "bassBasis")
theta_star <- c(14, 0.3)
y <- y_m + y_s * (K %*% (coef %*% c(1, theta_star)))[, 1] + rnorm(m, sd = 0.1)
bounds <- cbind(c(10, 20), c(-1, 1))
groups <- rep(1:2, each = m / 2)
D_narrow <- cbind(1, seq(-1, 1, length.out = m))
calibrate <- function(type, discrep.mat, seed) {
  set.seed(seed)
  bass$calibrate.bassBasis(y, emulator, type = type, sd.est = c(0.1, 0.1), s2.df = c(10, 10), s2.ind = groups,
                           meas.error.cor = diag(m), bounds = bounds, discrep.mean = rep(0, m), discrep.mat = discrep.mat,
                           nmcmc = 2, temperature.ladder = 1.05^(0:7), decor.step.every = 1000, verbose = FALSE)
}

# ---- stale likelihood, coordinates, unit mismatch (functional branch) ---------------------------------------------
seen <- list()
run <- calibrate(1, D_narrow, 32)
ntemps <- dim(run$theta)[2]
recomputed <- vapply(seq_len(ntemps), function(t) {
  cur <- run$curr[[t]]
  r <- y - cur$pred - c(cur$discrep)
  c(-0.5 * (cur$cov$ldet + t(r) %*% cur$cov$inv %*% r))
}, 0)
cached <- vapply(run$curr, function(cur) c(cur$llik), 0)
rejected <- vapply(seq_len(ntemps), function(t) all(run$theta[2, t, ] == run$theta[1, t, ]), TRUE)
stale_likelihood <- list(rungs = ntemps, theta_rejected = rejected, cached_minus_recomputed = cached - recomputed,
                         discrep_nonzero = vapply(run$curr, function(cur) any(c(cur$discrep) != 0), TRUE))

# The emulator's first call predicts at the starting theta of every rung, which R returns as theta[1, , ].
first_call <- seen[[1]]
start <- matrix(run$theta[1, , ], ntemps)
corrected_start <- sweep(sweep(start, 2, bounds[2, ] - bounds[1, ], "*"), 2, bounds[1, ], "+")
coordinates <- list(bounds = bounds, r_returned_range = apply(run$theta, 3, range),
                    emulator_input_range = apply(do.call(rbind, seen), 2, range),
                    r_start_vs_first_call = max(abs(start - first_call)),
                    corrected_start_vs_first_call = max(abs(corrected_start - first_call)))

v_trunc <- diag(cov(t(emulator$dat$trunc.error)))
ldet <- function(S) 2 * sum(log(diag(chol(S))))
unit_mismatch <- lapply(seq_len(ntemps), function(t) {
  cur <- run$curr[[t]]
  noise <- cur$s2[groups]
  r_units <- diag(v_trunc + noise) + K %*% diag(cur$emu.vars) %*% t(K)
  physical <- diag(y_s^2 * v_trunc + noise) + (y_s * K) %*% diag(cur$emu.vars) %*% t(y_s * K)
  list(r_ldet = c(cur$cov$ldet), dense_r_units = ldet(r_units), dense_physical = ldet(physical))
})

# ---- branch by shape ------------------------------------------------------------------------------------------------
branch <- function(res) if (is.null(res$dat$D)) "multivariate" else "functional"
branch_by_shape <- list(
  type1_square = list(declared = "functional", r_branch = branch(calibrate(1, diag(m), 33))),
  type2_narrow = list(declared = "multivariate", r_branch = branch(calibrate(2, D_narrow, 34))),
  type1_narrow = list(declared = "functional", r_branch = branch(run)))

# ---- s2.df switch ---------------------------------------------------------------------------------------------------
switch_stmt <- source_calls(bass$calibrate.bassBasis, function(e) identical(e[[1]], as.name("if")) &&
                              identical(deparse(e[[2]]), "any(s2.df == 0)"))[[1]]
kernel_for <- function(s2.df) { e <- new.env(); assign("s2.df", s2.df, e); eval(switch_stmt, e); e$ldig.kern }
sd_est <- c(0.1, 0.1)
s2_grid <- c(0.001, 0.01, 0.1, 1)
a <- c(0, 10) / 2; b <- a * sd_est^2
s2_df_switch <- list(s2_df = c(0, 10), sigma2 = s2_grid,
                     r_group2 = kernel_for(c(0, 10))(s2_grid, a[2], b[2]),
                     corrected_group2 = -(a[2] + 1) * log(s2_grid) - b[2] / s2_grid,
                     r_group2_when_all_positive = kernel_for(c(10, 10))(s2_grid, a[2], b[2]))

# ---- tempering ------------------------------------------------------------------------------------------------------
fn <- bass$calibrate.bassBasis
alpha_calls <- source_calls(fn, function(e) as.character(e[[1]])[1] %in% c("<-", "=") && identical(e[[2]], as.name("alpha")))
alpha_text <- vapply(alpha_calls, function(e) paste(deparse(e[[3]]), collapse = ""), "")
variance_alpha <- alpha_calls[[which(grepl("ldig.kern(cand$s2", alpha_text, fixed = TRUE))[1]]][[3]]
swap_alpha <- alpha_calls[[which(grepl("sw[2]", alpha_text, fixed = TRUE))[1]]][[3]]
discrep_cov <- source_calls(fn, function(e) identical(e[[1]], as.name("rmnorm")) &&
                              grepl("discrep.S", paste(deparse(e), collapse = ""), fixed = TRUE))[[1]][[3]]
itl <- 1 / c(1, 2, 5)
ldig <- kernel_for(10)
a1 <- 5; b1 <- 5 * 0.1^2
e_emu <- c(0.01, 0.01)
resid0 <- y - y_m - y_s * c(K %*% (coef %*% c(1, theta_star)))
Sigma_at <- function(s2) diag(v_trunc + s2) + K %*% diag(e_emu) %*% t(K)
llik_at <- function(s2, v = c(0, 0)) {
  r <- resid0 - c(D_narrow %*% v); L <- chol(Sigma_at(s2))
  -sum(log(diag(L))) - 0.5 * sum(forwardsolve(t(L), r)^2)
}
u <- seq(log(1e-4), 0, length.out = 61)
ll <- vapply(exp(u), llik_at, 0)
lk <- ldig(exp(u), a1, b1)
state <- function(i) list(llik = ll[i], s2 = exp(u[i]))
r_variance <- function(t, i, j) {
  curr <- vector("list", 3); curr[[t]] <- state(i)
  c(eval(variance_alpha, list2env(list(itl = itl, t = t, cand = state(j), curr = curr, ldig.kern = ldig, a = a1, b = b1))))
}
walk_kernel <- function(accept) {                     # accept(i, j): log acceptance of the move i -> j
  n <- length(u); P <- matrix(0, n, n)
  for (i in seq_len(n)) for (j in c(i - 1, i + 1)) if (j >= 1 && j <= n) P[i, j] <- 0.5 * min(1, exp(accept(i, j)))
  diag(P) <- 1 - rowSums(P)
  P
}
normalize <- function(l) { w <- exp(l - max(l)); w / sum(w) }
tv <- function(p, q) 0.5 * sum(abs(p - q))
rungs <- lapply(1:3, function(t) {
  b <- itl[t]
  s_r <- stationary_gth(walk_kernel(function(i, j) r_variance(t, i, j)))
  s_c <- stationary_gth(walk_kernel(function(i, j) b * (ll[j] - ll[i]) + lk[j] - lk[i] + u[j] - u[i]))
  list(b = b, s_r = s_r, s_c = s_c,
       tv_r_vs_powered = tv(s_r, normalize(b * (ll + lk + u))), tv_r_vs_likelihood_path = tv(s_r, normalize(b * ll + lk + u)),
       tv_corrected_vs_likelihood_path = tv(s_c, normalize(b * ll + lk + u)))
})
alpha_linear_gap <- max(abs(vapply(1:60, function(i) r_variance(2, i, i + 1) - itl[2] * r_variance(1, i, i + 1), 0)))
r_swap <- function(t1, t2, i, j) {                   # rung t1 holds state i, rung t2 holds state j
  curr <- vector("list", 3); curr[[t1]] <- state(i); curr[[t2]] <- state(j)
  c(eval(swap_alpha, list2env(list(itl = itl, sw = c(t1, t2), curr = curr, ldig.kern = ldig, a = a1, b = b1))))
}
# Detailed balance for a swap between rung laws s1 (rung t1, state i) and s2 (rung t2, state j) needs the log ratio
# log s1(j) + log s2(i) - log s1(i) - log s2(j).
swap_gap <- function(s1, s2, ratio) {
  n <- length(u); g <- 0
  for (i in seq_len(n)) for (j in seq_len(n))
    g <- max(g, abs(ratio(i, j) - (log(s1[j]) + log(s2[i]) - log(s1[i]) - log(s2[j]))))
  g
}
swaps <- lapply(1:2, function(t1) {
  t2 <- t1 + 1
  lp <- function(t) normalize(itl[t] * ll + lk + u)
  list(b1 = itl[t1], b2 = itl[t2],
       r_vs_r_kernels = swap_gap(rungs[[t1]]$s_r, rungs[[t2]]$s_r, function(i, j) r_swap(t1, t2, i, j)),
       r_vs_likelihood_path = swap_gap(lp(t1), lp(t2), function(i, j) r_swap(t1, t2, i, j)),
       corrected_vs_corrected_kernels = swap_gap(rungs[[t1]]$s_c, rungs[[t2]]$s_c,
                                                 function(i, j) (itl[t1] - itl[t2]) * (ll[j] - ll[i])))
})
Sig <- Sigma_at(0.01)
G <- t(D_narrow) %*% solve(Sig, D_narrow)
S_d <- solve(diag(2) + G)
discrepancy <- lapply(2:3, function(t) {
  b <- itl[t]
  r_cov <- eval(discrep_cov, list2env(list(discrep.S = S_d, itl = itl, t = t)))
  list(b = b, r_vs_powered = max(abs(r_cov - S_d / b)) / max(abs(S_d / b)),
       r_vs_likelihood_path = max(abs(r_cov - solve(diag(2) + b * G))) / max(abs(S_d / b)))
})
# R's swap between states that differ only in v (sigma^2 = 0.01 in both), rungs 1 and 2.
v_pairs <- list(c(0, 0), c(1, -1), c(0.5, 2))
log_powered <- function(t, v) itl[t] * (llik_at(0.01, v) + ldig(0.01, a1, b1) + log(0.01) - 0.5 * sum(v^2))
v_gap <- max(vapply(v_pairs[-1], function(v2) {
  v1 <- v_pairs[[1]]
  curr <- vector("list", 3)
  curr[[1]] <- list(llik = llik_at(0.01, v1), s2 = 0.01); curr[[2]] <- list(llik = llik_at(0.01, v2), s2 = 0.01)
  r <- c(eval(swap_alpha, list2env(list(itl = itl, sw = c(1, 2), curr = curr, ldig.kern = ldig, a = a1, b = b1))))
  abs(r - (log_powered(1, v2) + log_powered(2, v1) - log_powered(1, v1) - log_powered(2, v2)))
}, 0))
tempering <- list(
  ladder = 1 / itl, log_sigma2_range = diff(range(u)), alpha_linear_gap = alpha_linear_gap,
  rungs = lapply(rungs, function(x) x[c("b", "tv_r_vs_powered", "tv_r_vs_likelihood_path", "tv_corrected_vs_likelihood_path")]),
  swaps = swaps, discrepancy_covariance = discrepancy, swap_vs_powered_v_prior_gap = v_gap)

tempering$grid <- list(u = u, llik = ll, log_prior_kernel = lk)
inputs <- list(y = y, trunc_error = emulator$dat$trunc.error, emu_s2 = sapply(emulator$mod.list, function(k) k$s2), coef = coef,
               bounds = bounds, groups = groups, D_narrow = D_narrow, sd_est = c(0.1, 0.1), s2_df = c(10, 10))
run_state <- lapply(seq_len(ntemps), function(t) {
  cur <- run$curr[[t]]
  list(theta_start = run$theta[1, t, ], theta_end = run$theta[2, t, ], ii = run$ii[2], s2 = cur$s2, emu_vars = cur$emu.vars,
       discrep_vars = c(cur$discrep.vars), pred = cur$pred, llik_cached = c(cur$llik), ldet = c(cur$cov$ldet))
})
jsonlite::write_json(list(m = m, theta_star = theta_star, stale_likelihood = stale_likelihood, branch_by_shape = branch_by_shape,
                          coordinates = coordinates, unit_mismatch = unit_mismatch, s2_df_switch = s2_df_switch,
                          tempering = tempering, inputs = inputs, run = run_state), args[2],
                     digits = I(17), auto_unbox = TRUE, pretty = TRUE)
