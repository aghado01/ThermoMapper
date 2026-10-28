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
#
# Emits { m, theta_star, stale_likelihood, branch_by_shape, coordinates, unit_mismatch, s2_df_switch }.

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

emit_json(list(m = m, theta_star = theta_star, stale_likelihood = stale_likelihood, branch_by_shape = branch_by_shape,
               coordinates = coordinates, unit_mismatch = unit_mismatch, s2_df_switch = s2_df_switch), args[2])
