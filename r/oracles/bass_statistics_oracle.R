# bass_statistics_oracle.R — R BASS's own sufficient statistics for a design (port brief, stage 4's S4 and S5): the
# parity material where R is right (the hub's parity rule), under R's uncentred convention with the intercept first.
#
# Usage (cwd = r/):  Rscript oracles/bass_statistics_oracle.R <design.csv> <out.json> [scalar] [tau=<t>]
#                    Rscript oracles/bass_statistics_oracle.R <A.csv> <out.json> functional <H.csv> <Y.csv> [tau=<t>]
#   scalar: headerless CSV, y in the first column, then the design B with the intercept's column of ones first.
#   functional: A the case-side factors (n x k) and H the output-side factors (m x k), each with the intercept's
#     column of ones first, and Y the n x m response; column j of the design is A_j (x) H_j, rows (case, location)
#     with the case outer, and y = vec(Y') in that order (bass_admission_oracle.R's hadamard convention).
#   tau: the fixed tau of s2_rate (default 1/N).
#
# R's statements are evaluated as written (source_calls), never transcribed: data$n and data$ssy from bass(); the
# intercept's start, curr.list[[i]]$Xty[1] and $XtX[1,1], and the case-side and output-side bases' columns of ones,
# from bass(); then one birth per column after the intercept, in the design's order, through birth_des's (scalar) or
# birth_des_func's (functional) ata, Xta and aty and their three assignments into curr$Xty and curr$XtX, which fill
# the upper triangle only, as in R. Between births, addBasis's nbasis and nc and addBasisDes's and addBasisFunc's
# cbind. R's getQf is then read twice: as R calls it, for its verdict (its pivot-ratio rule, max_{j >= 2} r_jj /
# min_j r_jj > 1e3, a defect fixture under D06, emitted as a flag and not used), and through its own R, bhat and qf
# statements with the rule left out, so the statistics exist whatever the verdict. s2_rate is updateMCMC's
# curr$s2.rate expression at tau, (y'y - qf/(1 + tau))/2, which cancels when qf is close to y'y (D06);
# s2_rate_direct is the same from the residual, (|r|^2 + qf tau/(1 + tau))/2.
#
# Emits { mode, n, m, N, k, XtX (k x k, symmetric), Xty, yty, qf, bhat, R_diag (diag of R's chol factor),
#         log_det (log|X'X| = 2 sum log R_diag), positive_definite, getqf_refuses, getqf_ratio,
#         getqf_refuses_at (the first birth, 1-based over the columns after the intercept, at which R's getQf
#           refused; 0 when none did), tau, s2_rate, s2_rate_direct,
#         check { xtx_rel_error, xty_rel_error, qf_rel_error, rss_direct } } :
#   XtX, Xty: R's, against crossprod of the explicit design (n x k scalar, nm x k functional);
#   qf = Xty' (X'X)^-1 Xty, against a direct QR solve of the explicit design.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_statistics_oracle.R <design.csv> <out.json> [scalar | functional <H.csv> <Y.csv>] [tau=<t>]")
named <- length(args) >= 3 && !grepl("=", args[3], fixed = TRUE)
mode <- if (named) args[3] else "scalar"
if (!mode %in% c("scalar", "functional")) stop("mode must be scalar or functional")
if (mode == "functional" && length(args) < 5) stop("functional mode needs <H.csv> <Y.csv>")
rest <- args[-seq_len(if (mode == "functional") 5 else if (named) 3 else 2)]
keys <- list()
for (kv in rest) { parts <- strsplit(kv, "=", fixed = TRUE)[[1]]; if (parts[1] != "tau") stop("unknown key ", parts[1]); keys$tau <- as.numeric(parts[2]) }

bass <- bass_source(c("bass.R", "mcmc_funcs.R", "mcmc_update.R", "rjmcmc_all.R", "rjmcmc_des.R", "rjmcmc_des_func.R"))
first_assignment <- function(fn, lhs) {
  found <- source_calls(fn, function(e) as.character(e[[1]])[1] %in% c("<-", "=") && identical(paste(deparse(e[[2]]), collapse = ""), lhs))
  if (!length(found)) stop("no assignment to ", lhs)
  found[[1]]
}

if (mode == "scalar") {
  F <- read_matrix(args[1])
  y <- F[, 1]; A <- F[, -1, drop = FALSE]; H <- NULL
  n <- length(y); m <- 1
  if (!all(A[, 1] == 1)) stop("the design's first column must be the intercept's ones")
  X <- A; yv <- y
  data <- list(y = y, ndes = n, nfunc = 1)
  birth <- bass$birth_des
} else {
  A <- read_matrix(args[1]); H <- read_matrix(args[4]); Y <- read_matrix(args[5])
  n <- nrow(A); m <- nrow(H)
  if (ncol(A) != ncol(H)) stop("A and H need the same number of columns")
  if (!all(dim(Y) == c(n, m))) stop("Y must be n x m")
  if (!all(A[, 1] == 1) || !all(H[, 1] == 1)) stop("A's and H's first columns must be the intercept's ones")
  X <- sapply(seq_len(ncol(A)), function(j) kronecker(A[, j], H[, j]))
  yv <- as.vector(t(Y))
  data <- list(y = Y, ndes = n, nfunc = m)
  birth <- bass$birth_des_func
}
k <- ncol(A)
env <- new.env(parent = bass)
env$data <- data
for (lhs in c("data$n", "data$ssy")) eval(first_assignment(bass$bass, lhs), env)
N <- env$data$n
tau <- if (is.null(keys$tau)) 1 / N else keys$tau

# The intercept's start, as bass() writes it into each chain.
env$i <- 1
env$curr.list <- list(list(nbasis = 0, nc = 1, Xty = rep(NA, k + 1), XtX = matrix(NA, k + 1, k + 1)))
for (lhs in c("curr.list[[i]]$des.basis", "curr.list[[i]]$Xty[1]", "curr.list[[i]]$XtX[1, 1]",
              if (mode == "functional") "curr.list[[i]]$func.basis"))
  eval(first_assignment(bass$bass, lhs), env)
env$curr <- env$curr.list[[1]]

# One birth per further column: R's statements in order, then the bookkeeping addBasis* does on acceptance.
birth_steps <- lapply(c("ata", "Xta", "aty", "curr$Xty[curr$nc + 1]", "curr$XtX[1:curr$nc, curr$nc + 1]",
                        "curr$XtX[curr$nc + 1, curr$nc + 1]", "qf.cand.list"), function(lhs) first_assignment(birth, lhs))
accept_steps <- c(lapply(c("curr$nbasis", "curr$nc"), function(lhs) first_assignment(bass$addBasis, lhs)),
                  list(first_assignment(bass$addBasisDes, "curr$des.basis")),
                  if (mode == "functional") list(first_assignment(bass$addBasisFunc, "curr$func.basis")))
refuses_at <- 0
for (j in seq_len(k)[-1]) {
  env$cand.des <- list(basis = A[, j])
  if (mode == "functional") env$cand.func <- list(basis = H[, j])
  for (e in birth_steps) eval(e, env)
  if (is.null(env$qf.cand.list) && refuses_at == 0) refuses_at <- j - 1
  env$cand <- env$cand.des
  for (e in accept_steps[1:3]) eval(e, env)
  if (mode == "functional") { env$cand <- env$cand.func; eval(accept_steps[[4]], env) }
}
XtX_upper <- env$curr$XtX[1:k, 1:k, drop = FALSE]
Xty <- env$curr$Xty[1:k]
XtX <- XtX_upper; XtX[lower.tri(XtX)] <- t(XtX_upper)[lower.tri(XtX)]

# getQf as R calls it, for its verdict; then its own statements without the rule, for the statistics.
verdict <- bass$getQf(XtX_upper, Xty)
qf_env <- new.env(parent = bass); qf_env$XtX <- XtX_upper; qf_env$Xty <- Xty
eval(call("<-", as.name("R"), source_assignment(bass$getQf, "R")), qf_env)
positive_definite <- isTRUE(as.logical(qf_env$R[1, 1])) && nrow(qf_env$R) == k
out <- list(mode = mode, n = n, m = m, N = N, k = k, XtX = XtX, Xty = Xty, yty = env$data$ssy,
            positive_definite = positive_definite, getqf_refuses = is.null(verdict), getqf_refuses_at = refuses_at, tau = tau)
if (positive_definite) {
  for (name in c("bhat", "qf")) eval(call("<-", as.name(name), source_assignment(bass$getQf, name)), qf_env)
  dr <- diag(qf_env$R)
  s2_env <- new.env(parent = bass)
  s2_env$data <- env$data
  s2_env$curr <- list(bhat = qf_env$bhat, Xty = Xty, nc = k, beta.prec = tau)
  s2_rate <- c(eval(first_assignment(bass$updateMCMC, "curr$s2.rate")[[3]], s2_env))
  fit <- qr(X, tol = 0)
  fitted <- qr.fitted(fit, yv)
  qf_direct <- sum(fitted^2); rss_direct <- sum((yv - fitted)^2)
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  out <- c(out, list(qf = c(qf_env$qf), bhat = c(qf_env$bhat), R_diag = dr, log_det = 2 * sum(log(dr)),
                     getqf_ratio = if (k > 1) max(dr[-1]) / min(dr) else 1,
                     s2_rate = s2_rate, s2_rate_direct = 0.5 * (rss_direct + qf_direct * tau / (1 + tau)),
                     check = list(xtx_rel_error = rel(XtX, crossprod(X)), xty_rel_error = rel(Xty, c(crossprod(X, yv))),
                                  qf_rel_error = abs(c(qf_env$qf) - qf_direct) / qf_direct, rss_direct = rss_direct)))
}
emit_json(out, args[2])
