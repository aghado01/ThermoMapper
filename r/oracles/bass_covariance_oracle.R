# bass_covariance_oracle.R — reference values for BASS calibration's covariance formulations (port brief, calibration
# segment, D09; C-O2).
#
# Usage (cwd = r/):  Rscript oracles/bass_covariance_oracle.R <basis.csv> <out.json>
#   basis: headerless CSV, one row per output location: y_m, y_s, then the basis columns K_1, K_2 (m x 4, m even),
#     with y_s not all 1.
#
# Inputs built here (seeded): truncation residuals (m x 10) with variances v_trunc and covariance C_trunc; two variance
# groups (first and second halves of the outputs); emulator component variances e; a vector v to solve against.
# U = y_s * K is the basis in physical units (D07). Sections:
#   functional: Sigma = diag(d) + U diag(e) U', d = y_s^2 v_trunc + s2_g, in three cases: well (e = 0.3, 0.7),
#     moderate (every fourth d at 1e-4, e = 100, 300) and extreme (every fourth d at 1e-8, e = 1e4, 3e4). The
#     reference is the dense Cholesky solve and log-determinant, checked against a dense symmetric eigendecomposition
#     (a second algorithm); where the two disagree the case is past what a dense reference can certify. Against it:
#     R's Woodbury helpers swm and swm.ldet (calib.R), evaluated as written, which form m x m inverses; and D09's
#     low-rank route, the Cholesky factor of I + E^{1/2} U' D^-1 U E^{1/2}, with log|Sigma| = sum log d + log|that|.
#     backward_error: each route's solve x judged by the structured residual, which needs no reference:
#     max_i |d_i x_i + (U E U' x)_i - v_i| / (|d_i x_i| + (|U| E |U'| |x|)_i + |v_i|).
#   multivariate: Sigma = diag(y_s) C_trunc diag(y_s) + diag(sigma_g) R diag(sigma_g) + Sigma_delta + U diag(e) U', with
#     R an AR(1) correlation (0.6) within groups and Sigma_delta a squared-exponential covariance; dense Cholesky
#     against eigendecomposition.
#   units: the posterior over theta_1 on a 201-point grid (uniform prior), emulator w(theta) = (sin(theta/3),
#     cos(theta/4)), y = y_m + y_s * (K w(14)) + noise. physical: log N(y; y_m + y_s * K w, Sigma) with the
#     functional Sigma above; standardized: log N((y - y_m) / y_s; K w, diag(1/y_s) Sigma diag(1/y_s)). They differ
#     by the constant sum log y_s, so the normalized posteriors agree. For contrast, R's mismatch (C-O1): the
#     physical residual with the covariance in R's units, diag(v_trunc + s2_g) + K diag(e) K'.
#
#   inputs: what the C# parity tests rebuild the covariances from: the residuals te (m x 10), the group variances, the
#     vector v, each functional case's d and e, the multivariate case's R and Sigma_delta, and the units section's y and
#     log-likelihoods on the grid. Written at 17 significant digits, so every double round-trips.
#
# Emits { m, functional: [ { case, condition, reference: { ldet, solve }, eigen_vs_reference, r_woodbury, low_rank,
#         backward_error } ],
#         multivariate: { ldet, solve, eigen_vs_reference }, units: { ... }, inputs: { ... } }; errors are relative (ldet
#         absolute).

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_covariance_oracle.R <basis.csv> <out.json>")
Bm <- read_matrix(args[1])
y_m <- Bm[, 1]; y_s <- Bm[, 2]; K <- Bm[, 3:4]
m <- nrow(Bm)
if (m %% 2 != 0) stop("m must be even")
bass <- bass_source(c("calib.R"))

set.seed(41)
te <- matrix(rnorm(m * 10, sd = 0.05), m, 10)
v_trunc <- apply(te, 1, var)
C_trunc <- cov(t(te))
groups <- rep(1:2, each = m / 2)
s2_g <- c(0.02, 0.08)
U <- y_s * K
v <- rnorm(m)

# R's helpers, evaluated as written in one environment (swm calls my.solve).
r_env <- new.env()
for (name in c("my.solve", "swm", "swm.ldet"))
  assign(name, eval(source_assignment(bass$calibrate.bassBasis, name), r_env), envir = r_env)

dense <- function(S) {
  L <- chol(S)
  list(ldet = 2 * sum(log(diag(L))), solve = backsolve(L, forwardsolve(t(L), v)))
}
eigen_route <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  list(ldet = sum(log(e$values)), solve = c(e$vectors %*% (crossprod(e$vectors, v) / e$values)))
}
rel_err <- function(a, ref) max(abs(a - ref)) / max(abs(ref))
compare <- function(route, ref) list(ldet_abs = abs(route$ldet - ref$ldet), solve_rel = rel_err(route$solve, ref$solve))

functional_case <- function(name, d, e) {
  S <- diag(d) + U %*% diag(e) %*% t(U)
  ref <- dense(S)
  inv <- r_env$swm(diag(1 / d), U, diag(1 / e), t(U))
  r_w <- list(ldet = c(r_env$swm.ldet(diag(1 / d), U, diag(1 / e), t(U), sum(log(d)), sum(log(e)))), solve = c(inv %*% v))
  Ue <- sweep(U, 2, sqrt(e), "*")
  L <- chol(diag(length(e)) + crossprod(Ue / d, Ue))
  Dv <- v / d
  low <- list(ldet = sum(log(d)) + 2 * sum(log(diag(L))),
              solve = Dv - (Ue / d) %*% backsolve(L, forwardsolve(t(L), crossprod(Ue, Dv))))
  low$solve <- c(low$solve)
  eig <- eigen_route(S)
  backward <- function(x) {
    r <- d * x + c(U %*% (e * crossprod(U, x))) - v
    max(abs(r) / (abs(d * x) + c(abs(U) %*% (e * crossprod(abs(U), abs(x)))) + abs(v)))
  }
  list(case = name, condition = kappa(S, exact = TRUE), reference = ref, eigen_vs_reference = compare(eig, ref),
       r_woodbury = compare(r_w, ref), low_rank = compare(low, ref),
       backward_error = list(cholesky = backward(ref$solve), eigen = backward(eig$solve), r_woodbury = backward(r_w$solve),
                             low_rank = backward(low$solve)))
}
d_well <- y_s^2 * v_trunc + s2_g[groups]
with_small <- function(level) { d <- d_well; d[seq(1, m, by = 4)] <- level; d }
case_inputs <- list(well = list(d = d_well, e = c(0.3, 0.7)), moderate = list(d = with_small(1e-4), e = c(100, 300)),
                    extreme = list(d = with_small(1e-8), e = c(1e4, 3e4)))
functional <- lapply(names(case_inputs), function(name) functional_case(name, case_inputs[[name]]$d, case_inputs[[name]]$e))

ar1 <- function(n, rho) rho^abs(outer(seq_len(n), seq_len(n), "-"))
R_groups <- matrix(0, m, m)
for (g in 1:2) R_groups[groups == g, groups == g] <- ar1(sum(groups == g), 0.6)
tt <- seq(0, 1, length.out = m)
Sigma_delta <- 0.01 * exp(-outer(tt, tt, "-")^2 / (2 * 0.2^2)) + diag(1e-10, m)
S_mult <- diag(y_s) %*% C_trunc %*% diag(y_s) + diag(sqrt(s2_g[groups])) %*% R_groups %*% diag(sqrt(s2_g[groups])) +
  Sigma_delta + U %*% diag(c(0.3, 0.7)) %*% t(U)
ref_mult <- dense(S_mult)
multivariate <- c(ref_mult, list(condition = kappa(S_mult, exact = TRUE), eigen_vs_reference = compare(eigen_route(S_mult), ref_mult)))

# ---- units ----------------------------------------------------------------------------------------------------------
w_of <- function(theta) c(sin(theta / 3), cos(theta / 4))
e_units <- c(0.3, 0.7)
S_phys <- diag(d_well) + U %*% diag(e_units) %*% t(U)
S_std <- diag(1 / y_s) %*% S_phys %*% diag(1 / y_s)
S_r <- diag(v_trunc + s2_g[groups]) + K %*% diag(e_units) %*% t(K)
y <- y_m + y_s * c(K %*% w_of(14)) + c(t(chol(S_phys)) %*% rnorm(m))
loglik <- function(r, S) { L <- chol(S); -sum(log(diag(L))) - 0.5 * sum(forwardsolve(t(L), r)^2) }
grid <- seq(10, 20, length.out = 201)
ll_phys <- vapply(grid, function(th) loglik(y - y_m - y_s * c(K %*% w_of(th)), S_phys), 0)
ll_std <- vapply(grid, function(th) loglik((y - y_m) / y_s - c(K %*% w_of(th)), S_std), 0)
ll_r <- vapply(grid, function(th) loglik(y - y_m - y_s * c(K %*% w_of(th)), S_r), 0)
post <- function(ll) { w <- exp(ll - max(ll)); w / sum(w) }
units <- list(theta = grid, posterior_physical = post(ll_phys), posterior_standardized = post(ll_std),
              posterior_r_units = post(ll_r),
              loglik_difference_spread = diff(range(ll_std - ll_phys)),
              loglik_difference_mean = mean(ll_std - ll_phys), sum_log_y_s = sum(log(y_s)),
              max_abs_posterior_physical_vs_standardized = max(abs(post(ll_phys) - post(ll_std))),
              tv_physical_vs_r_units = 0.5 * sum(abs(post(ll_phys) - post(ll_r))),
              posterior_mean = c(physical = sum(grid * post(ll_phys)), standardized = sum(grid * post(ll_std)),
                                 r_units = sum(grid * post(ll_r))))

units$y <- y
units$loglik_physical <- ll_phys
units$loglik_standardized <- ll_std
units$loglik_r_units <- ll_r
inputs <- list(te = te, s2_g = s2_g, groups = groups, v = v, cases = case_inputs, R_groups = R_groups, Sigma_delta = Sigma_delta,
               e_units = e_units)
jsonlite::write_json(list(m = m, functional = functional, multivariate = multivariate, units = units, inputs = inputs), args[2],
                     digits = I(17), auto_unbox = TRUE, pretty = TRUE)
