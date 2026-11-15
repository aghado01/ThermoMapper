# bppr_oracle.R — BayesPPR 0.2.0's ridge basis and birth ratio on fixed draws, for E36's parity (the BPPR brief, slice 3; decision 8).
#
# Usage (cwd = r/):  Rscript oracles/bppr_oracle.R <fixture.json> <out.json>
#   fixture: { x: [[raw values of input 0], ...], y: [...], K, var_coefs, n_act_max, n_ridge_mean,
#              ridges: [{ feat: [0-based inputs], dir: [weights], knot0 }], born: { feat, dir, knot0 } }
#   Every input has more than K distinct values, so BayesPPR reads each as continuous ("cont").
#
# Emits, from BayesPPR's own functions where it has them (get_mns_basis, get_qf_info, dwallenius, from its namespace) and its bppr()'s
# arithmetic transcribed where it inlines it (the standardization, the knots, the birth's ratio terms):
#   x_st: the standardized columns; ridges and born: each one's knots and basis columns;
#   qf, qf_prop: yᵀB(BᵀB)⁻¹Bᵀy on the intercept and the ridges, without and with the born one (NA where get_qf_info refuses);
#   lik_ratio: the birth's collapsed likelihood ratio, −n/2·(ln sse′ − ln sse) − K/2·ln(var_coefs + 1);
#   log_act_feat: the birth's active-count and variable-set term with adapt_act_feat = TRUE, its weights 1 plus the current counts;
#   log_count: ln(n_ridge_mean/(M + 1)), the Poisson count ratio.

source("oracles/_common.R")
args <- commandArgs(trailingOnly = TRUE)
fx <- jsonlite::read_json(args[1], simplifyVector = FALSE)
ns <- asNamespace("BayesPPR")
if (as.character(utils::packageVersion("BayesPPR")) != "0.2.0") stop("BayesPPR is not the pinned 0.2.0")

X <- do.call(cbind, lapply(fx$x, function(v) as.numeric(unlist(v))))
y <- as.numeric(unlist(fx$y))
n <- length(y)
p <- ncol(X)
K <- fx$K
g <- fx$var_coefs

# bppr()'s standardization for a continuous input: X_st[, j] <- (X[, j] - mn_X[j])/sd_X[j].
X_st <- X
for (j in 1:p) X_st[, j] <- (X[, j] - mean(X[, j])) / sd(X[, j])

knot_quants <- seq(0, 1, length.out = K + 1)
ridge <- function(r) {
  feat <- as.integer(unlist(r$feat)) + 1
  dir <- matrix(as.numeric(unlist(r$dir)))
  proj <- X_st[, feat, drop = FALSE] %*% dir
  knots <- c(r$knot0, quantile(proj[proj > r$knot0], knot_quants))
  list(feat = feat, knots = as.numeric(knots), basis = ns$get_mns_basis(proj, knots))
}
current <- lapply(fx$ridges, ridge)
born <- ridge(fx$born)

design <- function(rs) do.call(cbind, c(list(rep(1, n)), lapply(rs, function(r) r$basis)))
quadratic <- function(B) {
  info <- ns$get_qf_info(t(B) %*% B, t(B) %*% y)
  if (is.null(info)) NA_real_ else as.numeric(info$qf)
}
qf <- quadratic(design(current))
qf_prop <- quadratic(design(c(current, list(born))))
ssy <- sum(y^2)
c_var <- g / (g + 1)
sse <- ssy - c_var * qf
sse_prop <- ssy - c_var * qf_prop
lik_ratio <- -n / 2 * (log(sse_prop) - log(sse)) - K / 2 * log(g + 1)

# The birth's adaptive terms (bppr(), adapt_act_feat = TRUE): weights 1 plus each accepted ridge's count.
n_act_max <- fx$n_act_max
w_n_act <- rep(1, n_act_max)
w_feat <- rep(1, p)
for (r in current) {
  w_n_act[length(r$feat)] <- w_n_act[length(r$feat)] + 1
  w_feat[r$feat] <- w_feat[r$feat] + 1
}
a <- length(born$feat)
log_act_feat <- -(log(n_act_max) + log(w_n_act[a] / sum(w_n_act)))
if (a > 1) log_act_feat <- log_act_feat - (lchoose(p, a) + log(ns$dwallenius(w_feat, born$feat)))
log_count <- log(fx$n_ridge_mean / (length(current) + 1))

columns <- function(B) lapply(seq_len(ncol(B)), function(k) as.numeric(B[, k]))
emit_json(list(
  x_st = columns(X_st),
  ridges = lapply(current, function(r) list(knots = r$knots, basis = columns(r$basis))),
  born = list(knots = born$knots, basis = columns(born$basis)),
  qf = qf, qf_prop = qf_prop, sse = sse, sse_prop = sse_prop, lik_ratio = lik_ratio,
  log_act_feat = log_act_feat, log_count = log_count), args[2])
