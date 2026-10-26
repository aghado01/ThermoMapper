# bass_stage1_checks.R — reproduces the BASS port brief's stage-1 checks (D03 propriety, D06 admission under both
# intercept conventions, D07 truncation residuals) by running the BASS oracles on fixed fixtures, exactly as the C#
# bridge would.
#
# Usage (cwd = r/):  Rscript scripts/bass_stage1_checks.R

rscript <- file.path(R.home("bin"), "Rscript")
write_fixture <- function(fixture) {
  csv <- tempfile(fileext = ".csv")
  write.table(fixture, csv, sep = ",", row.names = FALSE, col.names = FALSE)
  csv
}
run <- function(oracle, fixture, ...) {
  csv <- write_fixture(fixture); out <- tempfile(fileext = ".json")
  status <- system2(rscript, c(file.path("oracles", oracle), csv, out, ...))
  if (status != 0) stop(oracle, " failed")
  jsonlite::read_json(out, simplifyVector = TRUE)
}
hinge <- function(x, t, s) pmax(s * (x - t), 0)
show_marginal <- function(label, res) {
  cat(sprintf("%-40s k=%2d exact=%-5s log integral: %s | slope %6.2f (exact-fit rate %5.1f)\n", label, res$k,
              res$exact_fit, paste(sprintf("%8.2f", res$log_integral), collapse = " "), res$slope,
              res$predicted_slope))
}

# D03: a 40-point design with two opposite hinges at one knot, which with the intercept span 1 and x exactly.
set.seed(1)
N <- 40
x <- runif(N)
B3 <- cbind(1, hinge(x, 0.5, 1), hinge(x, 0.5, -1))
y_noisy <- 1 + 2 * x + rnorm(N, sd = 0.1)
y_exact <- 1 + 2 * x
BN <- cbind(1, sapply(sort(x)[2:(N - 1)], function(t) hinge(x, t, 1)), x)    # k = N: interpolates any y

cat("D03 propriety: lower limit of the tau integral from 1e-4 to 1e-64\n")
show_marginal("noisy y, improper sigma^2", run("bass_marginal_oracle.R", cbind(y_noisy, B3)))
show_marginal("exact fit, improper sigma^2", run("bass_marginal_oracle.R", cbind(y_exact, B3)))
show_marginal("exact fit, sigma^2 ~ IG(1e-3, 1e-3)", run("bass_marginal_oracle.R", cbind(y_exact, B3), "r", "1e-3", "1e-3"))
show_marginal("k = N, improper sigma^2", run("bass_marginal_oracle.R", cbind(y_noisy, BN)))
show_marginal("Liang, noisy y", run("bass_marginal_oracle.R", cbind(y_noisy, B3[, -1]), "liang"))
show_marginal("Liang, exact fit", run("bass_marginal_oracle.R", cbind(y_exact, B3[, -1]), "liang"))

# D06: three columns whose newest-column angle depends on insertion order.
e <- diag(3)
B_order <- cbind(e[, 1], e[, 1] + 1e-3 * e[, 2], e[, 2] + 1e-3 * e[, 3])
cat("\nD06 order dependence (s_min = 1e-10)\n")
for (ord in list(c(1, 2, 3), c(2, 3, 1))) {
  res <- run("bass_admission_oracle.R", B_order[, ord])
  cat(sprintf("  order %s: newest column's s = %.3g, leave-one-out min = %.3g\n", paste(ord, collapse = ""),
              tail(res$in_order, 1), res$min_loo))
}

# D06: R's getQf rule refuses an orthogonal local term once n makes the other columns' norms large.
n <- 40000
xl <- seq(0, 1, length.out = n)
local <- numeric(n); local[seq(n / 2 - 9, n / 2 + 10)] <- 0.01                 # 20 nonzero points, small values
B_local <- cbind(1, xl, local)
res <- run("bass_admission_oracle.R", B_local)
cat(sprintf("\nD06 R's getQf rule at n = %d: ratio %.0f, refuses = %s; leave-one-out min = %.3f\n", n,
            res$getqf_ratio, res$getqf_refuses, res$min_loo))

# D06 under the Liang convention (F-O7): the gate reads the centred Gram.
verdict <- function(ok) if (isTRUE(ok)) "pass" else "FAIL"
cat("\nD06 under the Liang convention (F-O7)\n")
set.seed(4)
N <- 500
x1 <- runif(N); x2 <- runif(N)
B_hinge <- cbind(1, x1, hinge(x1, 0.2, 1), hinge(x1, 0.6, -1), hinge(x2, 0.5, 1) * hinge(x1, 0.3, 1),
                 as.numeric(x2 > 0.05))
cen <- run("bass_admission_oracle.R", B_hinge, "centred")
unc <- run("bass_admission_oracle.R", B_hinge)
gap <- max(abs(cen$loo_uncentred / (cen$loo * cen$rho) - 1))
cross <- max(abs(unc$loo[-1] / cen$loo_uncentred - 1))
cat(sprintf("  hinge design, N = %d: rho %s\n", N, paste(sprintf("%.3f", cen$rho), collapse = " ")))
cat(sprintf("    uncentred s = centred s * rho: %.0e %s; uncentred mode against the QR angles: %.0e %s\n",
            gap, verdict(gap <= 1e-12), cross, verdict(cross <= 1e-10)))

# A categorical indicator true for all but k cases, nearly collinear with a second column once both are centred.
near_constant <- function(N, k, eta) {
  set.seed(5)
  b1 <- rep(1, N); b1[sample(N, k)] <- 0
  x <- runif(N)
  cbind(1, b1, (1 - b1) + eta * sqrt(k) * rnorm(N), x, hinge(x, 0.4, 1))
}
s_min <- 1e-10
cat(sprintf("  near-constant column, N = 1e5, s_min = %g\n", s_min))
for (k in c(1, 1000)) {
  res <- run("bass_admission_oracle.R", near_constant(1e5, k, 1e-7), "centred")
  cat(sprintf("    k = %4d (rho %.0e): centred min s %.2e, uncentred %.2e; rel. error downdate %.1e, explicit %.1e (left to right %.1e)\n",
              k, res$rho[1], res$min_loo, min(res$loo_uncentred), res$downdate_rel_error, res$explicit_rel_error,
              res$explicit_recursive_rel_error))
  if (k == 1) {
    worse <- res$downdate_rel_error >= 100 * res$explicit_rel_error && res$downdate_rel_error < 0.5
    cat(sprintf("      Liang admits (%s), R's convention refuses (%s): %s; downdate >= 100x explicit, inside the band: %s\n",
                res$min_loo >= s_min, min(res$loo_uncentred) < s_min,
                verdict(res$min_loo >= s_min && min(res$loo_uncentred) < s_min), verdict(worse)))
  }
}

# Direct functional: the centred Hadamard form, on factors whose means are large against their spread. One term
# is design-only (H_j = 1) and one functional-only (A_j = 1).
set.seed(6)
n_cases <- 60; m_locs <- 80
A <- cbind(1, matrix(runif(n_cases * 3, 0.9, 1), n_cases), runif(n_cases, 0.9, 1), 1)
H <- cbind(1, matrix(runif(m_locs * 3, 0.95, 1), m_locs), 1, runif(m_locs, 0.95, 1))
had <- run("bass_admission_oracle.R", A, "hadamard", write_fixture(H))
cat(sprintf("  centred Hadamard form, n = %d, m = %d, min s %.2e: Gram %.0e, s %.0e %s (downdate: Gram %.0e, s %.0e)\n",
            had$n, had$m, had$min_loo, had$gram_three_term_error, had$three_term_rel_error,
            verdict(had$gram_three_term_error <= 1e-13 && had$three_term_rel_error <= 1e-10),
            had$gram_downdate_error, had$downdate_rel_error))

# D07 (O-O4): truncation residuals, on curves whose phase and amplitude move with two inputs.
curves <- function(n, m = 200) {
  t <- seq(0, 1, length.out = m)
  z <- matrix(runif(2 * n), n, 2)
  y <- t(sapply(seq_len(n), function(i) sin(2 * pi * (t - 0.3 * z[i, 1])) * (1 + z[i, 2]) + 0.5 * z[i, 2] * t^2))
  y + matrix(rnorm(n * m, 0, 0.01), n)
}
show_truncation <- function(label, Y, ...) {
  res <- run("bass_truncation_oracle.R", Y, ...)
  scale <- max(abs(sweep(Y, 2, colMeans(Y))))
  coords <- max(res$coords_vs_output_in_sample, res$coords_vs_output_loo) / scale
  ratio <- if (is.null(res$ratio)) "undefined" else sprintf("%.2f", res$ratio)
  cat(sprintf("  %-26s r=%2d of %2d: RMS in-sample %.4f, leave-one-out %.4f, ratio %-9s coords vs output %.0e %s\n",
              label, res$r, res$centred_rank, res$in_sample_rms, res$loo_rms, ratio, coords, verdict(coords <= 1e-10)))
  invisible(res)
}

cat("\nD07 truncation residuals (O-O4)\n")
set.seed(2)
Y5 <- curves(5)
show_truncation("n = 5, 99%", Y5)
for (n_cases in c(10, 40)) show_truncation(sprintf("n = %d, 99%%", n_cases), curves(n_cases))
full <- show_truncation("n = 5, full retention", Y5, "100")
cat(sprintf("    in-sample max |e| %.1e within tolerance %.1e: %s; leave-one-out ranks %s (capped at n - 2)\n",
            full$in_sample_max_abs, full$zero_tolerance, verdict(full$in_sample_max_abs <= full$zero_tolerance),
            paste(full$loo_rank, collapse = " ")))
Y2 <- curves(2)
two <- show_truncation("n = 2", Y2)
gap <- max(abs(two$loo[1, ] - (Y2[1, ] - Y2[2, ])), abs(two$loo[2, ] - (Y2[2, ] - Y2[1, ])))
cat(sprintf("    leave-one-out equals the difference of the two cases: %.0e %s\n", gap,
            verdict(gap <= 1e-12 * max(abs(Y2)))))
