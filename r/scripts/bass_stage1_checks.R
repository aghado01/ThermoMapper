# bass_stage1_checks.R — reproduces the BASS port brief's stage-1 checks (D03 propriety, D06 admission) by
# running the two BASS oracles on fixed fixtures, exactly as the C# bridge would.
#
# Usage (cwd = r/):  Rscript scripts/bass_stage1_checks.R

rscript <- file.path(R.home("bin"), "Rscript")
run <- function(oracle, fixture, ...) {
  csv <- tempfile(fileext = ".csv"); out <- tempfile(fileext = ".json")
  write.table(fixture, csv, sep = ",", row.names = FALSE, col.names = FALSE)
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
