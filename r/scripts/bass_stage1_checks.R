# bass_stage1_checks.R — reproduces the BASS port brief's stage-1 checks (D03 propriety, D06 admission under both
# intercept conventions, D07 truncation residuals, and the reference fixtures F-O5, F-O6, O-O1 to O-O3, C-O1 and C-O2)
# by running the BASS oracles on fixed fixtures, exactly as the C# bridge would.
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

# ---- reference fixtures that do not wait on the design pass (issues/regression/bass/reference-fixtures-brief.md) ------
# These oracles source the pinned R BASS checkout (BASS_R_SOURCE, default D:/aghado01/packages/BASS) at 2671014.

# F-O5 (D01): Wallenius probabilities against brute force over draw orders.
cat("\nD01 Wallenius variable-set probabilities (F-O5)\n")
show_wallenius <- function(label, z, orders, max_subsets = 2000, draws = 0) {
  res <- run("bass_wallenius_oracle.R", matrix(z, 1), paste(orders, collapse = ","), max_subsets, draws)$orders
  for (i in seq_len(nrow(res))) {
    o <- res[i, ]
    extra <- c(if (!is.na(o$sum_brute)) sprintf("sum %.15f", o$sum_brute),
               if (!is.na(o$max_rel_equal)) sprintf("vs 1/C(p,K) %.0e", o$max_rel_equal),
               if (!is.na(o$max_z_score)) sprintf("R's sampler max |z| %.2f", o$max_z_score))
    cat(sprintf("  %-24s K=%d (%4d subsets): R's dmwnchBass %.1e, quadrature %.1e  %s\n", if (i == 1) label else "",
                o$K, o$n_subsets, o$max_rel_r, o$max_rel_quadrature, paste(extra, collapse = ", ")))
  }
  invisible(res)
}
w_eq <- show_wallenius("equal weights, p = 6", rep(5, 6), 1:5, draws = 20000)
w_r <- show_wallenius("R-like counts, p = 8", c(5, 8, 12, 5, 6, 17, 5, 9), 1:5, draws = 20000)
w_un <- show_wallenius("strongly unequal, p = 6", 10^(-3:2), 1:5)
w_50 <- show_wallenius("equal weights, p = 50", rep(5, 50), 2:5, max_subsets = 200)
sums <- c(w_eq$sum_brute, w_r$sum_brute, w_un$sum_brute)
cat(sprintf("    brute force sums to 1 (%.0e) and equals 1/C(p,K) under equal weights (%.0e): %s\n", max(abs(sums - 1)),
            max(c(w_eq$max_rel_equal, w_50$max_rel_equal)),
            verdict(max(abs(sums - 1)) <= 1e-14 && max(c(w_eq$max_rel_equal, w_50$max_rel_equal)) <= 1e-14)))
z_max <- max(c(w_eq$max_z_score, w_r$max_z_score))
cat(sprintf("    R's sample() draws the brute-force law (max |z| %.2f over 20000 draws per order): %s\n", z_max,
            verdict(z_max < 4.5)))
r_range <- max(c(w_eq$max_rel_r[2:3], w_r$max_rel_r[2:3]))
cat(sprintf("    dmwnchBass in R's range (K = 2, 3; equal and R-like weights) to %.0e: %s; at K = 1, which R never calls, it returns %.0fx the probability: %s\n",
            r_range, verdict(r_range <= 1e-12), 1 + w_eq$max_rel_r[1],
            verdict(abs(w_eq$max_rel_r[1] - 1) <= 1e-12 && abs(w_r$max_rel_r[1] - 1) <= 1e-12)))
q_max <- max(c(w_eq$max_rel_quadrature, w_r$max_rel_quadrature, w_un$max_rel_quadrature, w_50$max_rel_quadrature))
cat(sprintf("    quadrature of the positive integrand within %.0e of brute force everywhere: %s\n", q_max,
            verdict(q_max <= 1e-11)))

# F-O6 (a) (D05): R's tempering and swap density on an enumerable structure space.
cat("\nD05 R's tempering and swap density (F-O6 a)\n")
set.seed(7)
grid12 <- expand.grid(x1 = c(0, 1 / 3, 2 / 3, 1), x2 = c(0, 0.5, 1))
grid12 <- rbind(grid12, grid12)
y_swap <- 1 + 1.5 * pmax(grid12$x1 - 1 / 3, 0) - pmax(0.5 - grid12$x2, 0) + rnorm(nrow(grid12), sd = 0.3)
show_swap <- function(label, sw) {
  s <- sw$swap
  cat(sprintf("  %s, N = %d, K_max = 2, M_max = 2, sigma^2 = 0.1, tau = 1/N, lambda = 2: %d terms (%d with >= %d nonzero values), %d states, %d moving transitions\n",
              label, sw$n, sw$n_terms, sw$n_supported, sw$options$npart, sw$n_states, s$n_transitions))
  cat(sprintf("    getQf agrees with the strict rank in both column orders (%d mismatches): %s; R's alpha is b times its b = 1 value (gap %.0e): %s\n",
              s$getqf_mismatches, verdict(s$getqf_mismatches == 0), s$alpha_linear_in_b_gap,
              verdict(s$alpha_linear_in_b_gap <= 1e-12)))
  cat(sprintf("    at b = 1 R's kernel with the state's move probabilities is the posterior (TV %.0e): %s; with R's fixed thirds, TV %.3f (F-O4's defect)\n",
              s$tv_r_isolated_b1_posterior, verdict(s$tv_r_isolated_b1_posterior <= 1e-12), s$tv_r_literal_b1_posterior))
  for (i in seq_len(nrow(s$rung))) {
    tv <- s$rung$tv_r_to[i, ]
    cat(sprintf("    b = %.1f: D05's kernel vs its target TV %.0e; R's kernel vs like_int %.3f, like_lin %.3f, pow_coll %.3f, pow_joint %.3f; reversibility gap R %.0e, D05 %.0e\n",
                s$rung$b[i], s$rung$tv_corrected_target[i], tv$like_int, tv$like_lin, tv$pow_coll, tv$pow_joint,
                s$rung$reversibility_r[i], s$rung$reversibility_corrected[i]))
  }
  for (i in seq_len(nrow(s$pairs))) {
    pr <- s$pairs[i, ]; fam <- pr$r_swap_vs_family
    cat(sprintf("    swap %.1f <-> %.1f, spread of log s1 - log s2 - (b1 - b2) lp: R's kernels %.2f; like_int %.2f, like_lin %.2f, pow_coll %.2f, pow_joint %.2f; D05's swap on its kernels %.0e\n",
                pr$b1, pr$b2, pr$r_swap_vs_r_kernels, fam$like_int, fam$like_lin, fam$pow_coll, fam$pow_joint,
                pr$corrected_swap_vs_corrected_kernels))
  }
  cat(sprintf("    R's swap fails balance against its own kernels and every candidate target (spread >= 0.1): %s; D05's swap balances (<= 1e-10): %s; D05's kernels reach their target and are reversible (<= 1e-12): %s\n",
              verdict(min(s$pairs$r_swap_vs_r_kernels, unlist(s$pairs$r_swap_vs_family)) >= 0.1),
              verdict(max(s$pairs$corrected_swap_vs_corrected_kernels) <= 1e-10),
              verdict(max(s$rung$tv_corrected_target, s$rung$reversibility_corrected) <= 1e-12)))
  invisible(s)
}
s2p <- show_swap("p = 2", run("bass_enumeration_oracle.R", cbind(y_swap, grid12$x1, grid12$x2), "swap"))
cat(sprintf("    with p = 2 every interaction draws the one pair, and R's tempered kernel is reversible (gap <= 1e-12): %s\n",
            verdict(max(s2p$rung$reversibility_r) <= 1e-12)))
set.seed(8)
grid27 <- expand.grid(x1 = c(0, 0.5, 1), x2 = c(0, 0.5, 1), x3 = c(0, 0.5, 1))
y_swap3 <- 1 + 1.5 * pmax(grid27$x1 - 0.5, 0) - pmax(0.5 - grid27$x2, 0) + 0.8 * grid27$x3 * grid27$x1 +
  rnorm(nrow(grid27), sd = 0.3)
s3p <- show_swap("p = 3", run("bass_enumeration_oracle.R", cbind(y_swap3, grid27$x1, grid27$x2, grid27$x3), "swap", "npart=4"))
cat(sprintf("    with p = 3 the Wallenius factor enters R's tempered acceptance as (q_rev / q_fwd)^b, and its kernel loses reversibility below b = 1 (gaps %s): %s\n",
            paste(sprintf("%.0e", s3p$rung$reversibility_r[s3p$rung$b < 1]), collapse = ", "),
            verdict(min(s3p$rung$reversibility_r[s3p$rung$b < 1]) >= 1e-9 && s3p$rung$reversibility_r[s3p$rung$b == 1] <= 1e-12)))

# F-O6 (b) (D03): R's Jeffreys mode against the flat and Jeffreys priors' coherent factors.
cat("\nD03 R's Jeffreys mode (F-O6 b)\n")
set.seed(3)
xj <- runif(40)
Bj <- cbind(1, hinge(xj, 0.3, 1), hinge(xj, 0.6, -1))
yj <- 1 + 2 * Bj[, 2] - Bj[, 3] + rnorm(40, sd = 0.2)
jf <- run("bass_jeffreys_oracle.R", cbind(yj, Bj), "0.04")
quad_gap <- max(abs(unlist(jf$log_I$quadrature) - unlist(jf$log_I$closed_form)))
cat(sprintf("  N = %d, sigma^2 = %.2f: the integral over a by quadrature against its closed form, %.0e: %s\n", jf$n, jf$sigma2,
            quad_gap, verdict(quad_gap <= 1e-10)))
for (i in seq_len(nrow(jf$birth))) {
  b <- jf$birth[i, ]
  cat(sprintf("    birth %d -> %d columns, log factor: R %.4f (its death %.4f), flat %.4f, Jeffreys %.4f, prior sigma^-2k |G|^3/4 %.4f\n",
              b$k, b$k + 1, b$r, b$r_death, b$flat, b$jeffreys, b$r_moves_prior))
}
ch <- jf$change
cat(sprintf("    change, log factor: R %.4f, flat %.4f, Jeffreys %.4f, prior sigma^-2k |G|^3/4 %.4f\n", ch$r, ch$flat,
            ch$jeffreys, ch$r_moves_prior))
sh <- jf$sigma_shape
cat(sprintf("    sigma^2 shape (k = 2): R %.4f, flat %.4f (closed %.1f), Jeffreys %.4f (closed %.1f), prior sigma^-(k-1) %.4f\n",
            sh$r, sh$flat, sh$closed_flat, sh$jeffreys, sh$closed_jeffreys, sh$r_sigma_prior))
moves_match <- max(abs(c(jf$birth$r - jf$birth$r_moves_prior, ch$r - ch$r_moves_prior, jf$birth$r + jf$birth$r_death)))
coherent_match <- max(abs(c(jf$birth$flat - jf$birth$closed_flat, jf$birth$jeffreys - jf$birth$closed_jeffreys,
                            sh$flat - sh$closed_flat, sh$jeffreys - sh$closed_jeffreys)))
cat(sprintf("    coherent factors match D03's closed forms (%.0e): %s; R's moves imply sigma^-2k |G|^3/4 (%.0e): %s; R's sigma^2 draw implies sigma^-(k-1) (%.0e): %s\n",
            coherent_match, verdict(coherent_match <= 1e-8), moves_match, verdict(moves_match <= 1e-8),
            abs(sh$r - sh$r_sigma_prior), verdict(abs(sh$r - sh$r_sigma_prior) <= 1e-8)))
gaps <- c(abs(jf$birth$r - jf$birth$flat), abs(jf$birth$r - jf$birth$jeffreys), abs(sh$r - sh$flat), abs(sh$r - sh$jeffreys))
cat(sprintf("    R matches neither coherent prior (smallest gap %.2f): %s\n", min(gaps), verdict(min(gaps) >= 0.1)))

# F-O6 (c) (D01): the count law each birth ratio implies under a constant likelihood.
cat("\nD01 the count law of each birth ratio (F-O6 c)\n")
grid4 <- expand.grid(x1 = c(0, 1), x2 = c(0, 1))
show_count <- function(label, cl) {
  cn <- cl$count
  cat(sprintf("  %s, two knots per variable, K_max = 2, M_max = %d, lambda = 2: %d terms, %d states\n", label,
              max(cn$M), cl$n_terms, cl$n_states))
  cat(sprintf("    count law M = 0..%d: stated %s | D01's ratio %s | KHAOS's %s | stated x 1/M! %s\n", max(cn$M),
              paste(sprintf("%.4f", cn$stated), collapse = " "), paste(sprintf("%.4f", cn$d01), collapse = " "),
              paste(sprintf("%.4f", cn$khaos), collapse = " "), paste(sprintf("%.4f", cn$stated_over_factorial), collapse = " ")))
  cat(sprintf("    D01's ratio gives the stated law (TV %.0e): %s; KHAOS's gives it times 1/M! (TV %.0e): %s, not the stated law (TV %.3f): %s\n",
              cn$tv_d01_stated, verdict(cn$tv_d01_stated <= 1e-12), cn$tv_khaos_over_factorial,
              verdict(cn$tv_khaos_over_factorial <= 1e-12), cn$tv_khaos_stated, verdict(cn$tv_khaos_stated >= 0.1)))
}
show_count("p = 2", run("bass_enumeration_oracle.R", cbind(0, grid4$x1, grid4$x2), "count", "M_max=3"))
grid8 <- expand.grid(x1 = c(0, 1), x2 = c(0, 1), x3 = c(0, 1))
show_count("p = 3", run("bass_enumeration_oracle.R", cbind(0, grid8$x1, grid8$x2, grid8$x3), "count", "M_max=2"))

# O-O1 (D07): R's output-path defects.
cat("\nD07 R's output-path defects (O-O1)\n")
set.seed(21)
t_out <- seq(0, 1, length.out = 40)
Y_out <- t(sapply(1:8, function(i) sin(2 * pi * (t_out + 0.1 * i)) +
                    (i %% 3 == 0) * 2 * exp(-((t_out - 0.2 - 0.07 * i) / 0.04)^2) + rnorm(40, sd = 0.02)))
od <- run("bass_output_defects_oracle.R", Y_out)
ts <- od$truncation_sampler
cat(sprintf("  truncation sampler: %d locations x %d cases, R stops (\"%s\"); 3 x %d, R draws cases %s, corrected %s: %s\n",
            ts$wide$locations, ts$wide$cases, ts$wide$r_error, ts$narrow$cases, paste(ts$narrow$r_cases_drawn, collapse = ","),
            paste(ts$narrow$corrected_cases_drawn, collapse = ","),
            verdict(!is.null(ts$wide$r_error) && length(ts$narrow$r_cases_drawn) == 3 &&
                    length(ts$narrow$corrected_cases_drawn) == ts$narrow$cases)))
nr <- od$noise_recycling
cat(sprintf("  component noise: R's spread across cases within a draw %g, correlation %.2f; corrected %.2f, %.2f: %s\n",
            nr$r_max_spread_across_cases, nr$r_correlation_cases_1_2, nr$corrected_max_spread_across_cases,
            nr$corrected_correlation_cases_1_2,
            verdict(nr$r_max_spread_across_cases == 0 && abs(nr$corrected_correlation_cases_1_2) < 0.2)))
rs <- od$residual
cat(sprintf("  residual (r = %d, y_s in [%.2f, %.2f]): R stores minus D07's e (gap %g); training reconstruction error R %.2f, corrected %.0e (scale %.2f); skewness D07 %.2f, R %.2f: %s\n",
            rs$r, rs$y_s_range[1], rs$y_s_range[2], rs$r_stored_is_negative_of_corrected, rs$r_reconstruction_max_error,
            rs$corrected_reconstruction_max_error, rs$scale, rs$skewness_corrected, rs$skewness_r,
            verdict(rs$r_stored_is_negative_of_corrected == 0 && rs$corrected_reconstruction_max_error <= 1e-12 * rs$scale &&
                    rs$r_reconstruction_max_error > 0.1)))
fr <- od$full_retention
cat(sprintf("  perc.var = 100: R stops (\"%s\"), its n.pc expression gives %s (last cumulative share - 100 = %g); corrected rank %d: %s\n",
            fr$r_error, if (is.null(fr$r_n_pc)) "NA" else fr$r_n_pc, fr$cumulative_share_last_minus_100, fr$corrected_rank,
            verdict(!is.null(fr$r_error) && is.null(fr$r_n_pc) && fr$corrected_rank == nrow(Y_out) - 1)))
cc <- od$constant_columns
cat(sprintf("  constant input: training %g, prediction %s, corrected: refused at entry; constant output: R stops (\"%s\", y_s = %g), corrected y_s = 1: %s\n",
            cc$input$r_training, paste(cc$input$r_prediction, collapse = " and "), cc$output$r_error, cc$output$r_y_s,
            verdict(cc$input$r_training == 0 && all(cc$input$r_prediction %in% c("NaN", "Inf")) && !is.null(cc$output$r_error))))
do <- od$draw_order
cat(sprintf("  draw order: requested %s, R returns %s; R's dims %s, %s, %s; nugget noise R on draw %s, corrected on draw %s: %s\n",
            paste(do$requested, collapse = ","), paste(do$r_rows, collapse = ","), paste(do$r_dim_4x3, collapse = "x"),
            paste(do$r_dim_4x1, collapse = "x"), paste(do$r_dim_1x1, collapse = "x"),
            paste(do$r_rows_with_noise, collapse = ","), paste(do$corrected_rows_with_noise, collapse = ","),
            verdict(!identical(do$r_rows, do$requested) && length(do$r_dim_4x1) == 1 &&
                    !identical(do$r_rows_with_noise, do$corrected_rows_with_noise))))
iq <- od$interval_quantiles
cat(sprintf("  95%% intervals: means spread, R (%.2f, %.2f) vs mixture (%.2f, %.2f); sds spread, R (%.2f, %.2f) vs mixture (%.2f, %.2f): %s\n",
            iq$r[[1]][1], iq$r[[1]][2], iq$mixture[[1]][1], iq$mixture[[1]][2], iq$r[[2]][1], iq$r[[2]][2],
            iq$mixture[[2]][1], iq$mixture[[2]][2], verdict(min(abs(unlist(iq$r) - unlist(iq$mixture))) > 0.5)))

# O-O2 (D08): R's sensitivity defects.
cat("\nD08 R's sensitivity defects (O-O2)\n")
t_sens <- seq(0, 1, length.out = 12)
sd8 <- run("bass_sensitivity_defects_oracle.R", cbind(10 + t_sens, 2 + 3 * t_sens, sin(pi * t_sens), cos(pi * t_sens)))
cw <- sd8$category_weights
cat(sprintf("  categories (a, b, c), D1 = {a}, D2 = {a, b}: R's P(D1) %.3f, P(D2) %.3f, P(D1 & D2) %.3f on both paths, whatever the weights\n",
            cw$r_direct$sub_cnt[1], cw$r_direct$sub_cnt[2], cw$r_direct$C2Cat_D1_D2))
for (i in seq_len(nrow(cw$cases))) {
  e <- cw$cases$enumerated[i, ]
  cat(sprintf("    enumerated under (%s): %.3f, %.3f, %.3f\n", paste(cw$cases$weights[[i]], collapse = ", "), e$P_D1, e$P_D2,
              e$P_D1_and_D2))
}
cat(sprintf("    two levels, P(a) = 0.9: R's E[1(a)^2] %.2f, enumerated %.2f; no moment function reads prior.cat: %s\n",
            cw$two_level$r_second_moment, cw$two_level$enumerated,
            verdict(!any(unlist(cw$static_reads_prior_cat)) && abs(cw$two_level$r_second_moment - 0.5) < 1e-15 &&
                    abs(cw$r_basis$C2Basis_D1_D2 - 1 / 3) < 1e-15)))
tb <- sd8$truncation_bounds$cases
cat("  hinge products under U(0.25, 0.75), on R's scale: basis path (C22Basis) | direct path (C2) | quadrature\n")
for (i in seq_len(nrow(tb))) cat(sprintf("    knots (%.1f, %.1f), signs (%+d, %+d): %.7f | %.7f | %.7f\n", tb$t[[i]][1], tb$t[[i]][2],
                                         tb$s[[i]][1], tb$s[[i]][2], tb$r_basis[i], tb$r_direct[i], tb$quadrature[i]))
cat(sprintf("    the direct path matches quadrature (%.0e): %s; the basis path gives 0.45 for the inventory's 0.1105903: %s\n",
            max(abs(tb$r_direct - tb$quadrature)), verdict(max(abs(tb$r_direct - tb$quadrature)) <= 1e-12),
            verdict(abs(tb$r_basis[1] - 0.45) < 1e-12 && abs(tb$quadrature[1] - 0.1105903) < 1e-7)))
du <- sd8$decomposition_units
cat(sprintf("  sobolBasis, y_s in [%.0f, %.0f]: Var.tot vs decomposition units %.0e, vs physical %.2f; effect variances %.0e vs %.2f; indices vs reference %.0e: %s\n",
            min(du$y_s), max(du$y_s), du$max_rel_r_vs_decomposition, du$max_rel_r_vs_physical,
            du$max_rel_effects_r_vs_decomposition, du$max_rel_effects_r_vs_physical, du$max_abs_indices_r_vs_reference,
            verdict(du$max_rel_r_vs_decomposition <= 1e-12 && du$max_rel_effects_r_vs_decomposition <= 1e-12 &&
                    du$max_abs_indices_r_vs_reference <= 1e-12 && du$max_rel_r_vs_physical > 0.5)))
ef <- sd8$effects
fm <- ef$first_moment
cat("  effects' first moment of a scaled hinge: C1.temp (effects, f0) | C1 (variances) | quadrature\n")
for (i in seq_len(nrow(fm))) cat(sprintf("    %-7s on [%.2f, %.2f], (%+d, %.1f): %.6f | %.6f | %.6f\n", fm$measure[i], fm$trunc[[i]][1],
                                         fm$trunc[[i]][2], fm$s[i], fm$t[i], fm$r_c1_temp[i], fm$r_c1[i], fm$quadrature[i]))
cat(sprintf("    C1 matches quadrature (%.0e) and C1.temp only under U(0, 1): %s\n", max(abs(fm$r_c1 - fm$quadrature)),
            verdict(max(abs(fm$r_c1 - fm$quadrature)) <= 1e-12 && abs(fm$r_c1_temp[5] - fm$quadrature[5]) <= 1e-12 &&
                    min(abs(fm$r_c1_temp[1:4] - fm$quadrature[1:4])) > 0.01)))
cat(sprintf("    sobol_des's makeBasis degree arguments: %s; a degree-2 effect curve as written is the degree-1 hinge (gap %g), %.2f from the degree-2 one: %s\n",
            paste(ef$makebasis_degree_args, collapse = ", "), ef$degree2_curve$r_max_abs_vs_degree1,
            ef$degree2_curve$r_max_abs_vs_degree2,
            verdict(all(ef$makebasis_degree_args == "1") && ef$degree2_curve$r_max_abs_vs_degree1 == 0)))

# O-O3 (D08): moments against quadrature and category enumeration.
cat("\nD08 moments against quadrature and enumeration (O-O3)\n")
first <- expand.grid(t = c(-0.1, 0.1, 0.3, 0.5, 0.8, 1.1), s = c(-1, 1), q = 1:3)
kn <- c(0.1, 0.3, 0.5, 0.8)
pairs_kn <- which(upper.tri(diag(4), diag = TRUE), arr.ind = TRUE)
cross <- do.call(rbind, lapply(1:3, function(q) do.call(rbind, lapply(seq_len(nrow(pairs_kn)), function(i)
  do.call(rbind, lapply(list(c(1, 1), c(1, -1), c(-1, 1), c(-1, -1)), function(sg)
    c(kn[pairs_kn[i, 1]], sg[1], kn[pairs_kn[i, 2]], sg[2], q)))))))
mo <- run("bass_moments_oracle.R", rbind(cbind(first$t, first$s, 0, 0, first$q), cross))
names_m <- mo$measures$name
cat(sprintf("  %d cases (first moments with knots inside and outside [0, 1]; cross moments with equal and distinct knots, all sign pairs; degrees 1-3) x %d measures\n",
            nrow(mo$cases), length(names_m)))
cat(sprintf("    adaptive quadrature vs composite Gauss-Legendre: %.0e: %s\n", mo$max_rel_quadrature_vs_gl,
            verdict(mo$max_rel_quadrature_vs_gl <= 1e-12)))
groups_m <- mo$groups
show_table <- function(title, covered, maxrel, rows = seq_along(names_m)) {
  cat(sprintf("    %s: covered cases (largest relative error) for %s\n", title, paste(groups_m, collapse = ", ")))
  for (k in seq_along(rows)) {
    j <- rows[k]
    cells <- vapply(seq_along(groups_m), function(g)
      if (covered[k, g] == 0) "    -          " else sprintf("%3d (%-8s)", covered[k, g], sprintf("%.0e", maxrel[k, g])), "")
    cat(sprintf("      %-23s %s%s\n", names_m[j], paste(cells, collapse = " "),
                if (!is.na(mo$r_error[j]) && title == "R as written") sprintf("  [%s]", mo$r_error[j]) else ""))
  }
}
show_table("R as written", mo$r_covered, mo$r_max_rel)
student <- grep("^student", names_m)
show_table("R's Student forms with a working 2F1", mo$repaired_covered[student, , drop = FALSE],
           mo$repaired_max_rel[student, , drop = FALSE], student)
cat(sprintf("    R's uniform moments (every degree) and normal moments (degree 1) match quadrature (%.0e): %s\n",
            mo$max_rel_r_vs_quadrature, verdict(mo$max_rel_r_vs_quadrature <= 1e-12)))
trivial_only <- all(vapply(student, function(j) all(mo$quadrature[!is.na(mo$r[, j]), j] == 0), TRUE))
cat(sprintf("    R's Student moments as written give a value only where a hinge misses the interval (the rest fail on hypergeo's complex or NA results): %s\n",
            verdict(trivial_only)))
g1 <- match(c("first_q1", "cross_q1"), groups_m); g23 <- match(c("cross_q2", "cross_q3"), groups_m)
rep_nu23 <- mo$repaired_max_rel[match(c("student_nu2", "student_nu3"), names_m), , drop = FALSE]
nu1 <- match("student_nu1", names_m)
cat(sprintf("    with a working 2F1: degree 1 matches at nu = 2, 3 (%.0e): %s; degrees 2-3 of the cross moment come out as the degree-1 value (errors >= %.2f): %s; nu = 1 divides by nu - 1: %s\n",
            max(rep_nu23[, g1]), verdict(max(rep_nu23[, g1]) <= 1e-11), min(rep_nu23[, g23]),
            verdict(min(rep_nu23[, g23]) >= 0.5),
            verdict(all(mo$quadrature[!is.na(mo$r_repaired[, nu1]), nu1] == 0))))
mc <- run("bass_moments_oracle.R", rbind(c(0, 0.2, 0.3, 0.5), c(1, 0, 0, 0), c(1, 1, 0, 0), c(0, 1, 1, 0), c(0, 0, 1, 1)),
          "categorical")
cat(sprintf("  categories with probabilities (%s): P(D) %s by enumeration, R's ratios %s; P(D1 & D3) %.2f (R %.2f): %s\n",
            paste(mc$probabilities, collapse = ", "), paste(mc$P, collapse = " "), paste(mc$r_P, collapse = " "),
            mc$P_joint[1, 3], mc$r_P_joint[1, 3], verdict(mc$P[1] == 0 && mc$P_joint[1, 3] == 0 && mc$r_P[1] == 0.25)))

# C-O1 (D09): R's calibration defects.
cat("\nD09 R's calibration defects (C-O1)\n")
t_cal <- seq(0, 1, length.out = 20)
basis_cal <- cbind(5 + t_cal, 0.5 + 2 * t_cal, sin(pi * t_cal), cos(pi * t_cal))
cd <- run("bass_calibration_defects_oracle.R", basis_cal)
sl <- cd$stale_likelihood
cat(sprintf("  stale likelihood, one iteration over %d rungs: theta rejected on rung %s; cached minus recomputed llik %s: %s\n",
            sl$rungs, paste(which(sl$theta_rejected), collapse = ","),
            paste(sprintf("%.3f", sl$cached_minus_recomputed), collapse = " "),
            verdict(any(sl$theta_rejected) && all(abs(sl$cached_minus_recomputed[!sl$theta_rejected]) <= 1e-10) &&
                    all(abs(sl$cached_minus_recomputed[sl$theta_rejected]) > 1e-3))))
bs <- cd$branch_by_shape
cat(sprintf("  branch: type 1 with a square discrepancy runs %s, type 2 with a narrow one %s, type 1 narrow %s: %s\n",
            bs$type1_square$r_branch, bs$type2_narrow$r_branch, bs$type1_narrow$r_branch,
            verdict(bs$type1_square$r_branch == "multivariate" && bs$type2_narrow$r_branch == "functional")))
co <- cd$coordinates
cat(sprintf("  coordinates: R returns theta_1 in [%.2f, %.2f] while the emulator saw [%.2f, %.2f]; starting theta vs the emulator's first call, R %.2f, corrected %.0e: %s\n",
            co$r_returned_range[1, 1], co$r_returned_range[2, 1], co$emulator_input_range[1, 1], co$emulator_input_range[2, 1],
            co$r_start_vs_first_call, co$corrected_start_vs_first_call,
            verdict(all(co$r_returned_range >= 0 & co$r_returned_range <= 1) && co$r_start_vs_first_call > 1 &&
                    co$corrected_start_vs_first_call <= 1e-12)))
um <- cd$unit_mismatch
cat(sprintf("  units: R's log-determinant vs dense in R's units %.0e, vs dense physical at least %.3f: %s\n",
            max(abs(um$r_ldet - um$dense_r_units)), min(abs(um$r_ldet - um$dense_physical)),
            verdict(max(abs(um$r_ldet - um$dense_r_units)) <= 1e-10 && min(abs(um$r_ldet - um$dense_physical)) > 0.1)))
s2s <- cd$s2_df_switch
cat(sprintf("  s2.df = (0, 10): group 2's log prior kernel at sigma^2 = %s is R %s against its IG %s; with (10, 10) R uses the IG: %s\n",
            paste(s2s$sigma2, collapse = ","), paste(sprintf("%.3f", s2s$r_group2), collapse = " "),
            paste(sprintf("%.2f", s2s$corrected_group2), collapse = " "),
            verdict(max(abs(s2s$r_group2 + log1p(s2s$sigma2))) <= 1e-15 &&
                    max(abs(s2s$r_group2_when_all_positive - s2s$corrected_group2)) <= 1e-12)))
tp <- cd$tempering
cat(sprintf("  tempering (ladder %s): R's variance alpha is b times likelihood + prior + Jacobian (gap %.0e): %s\n",
            paste(tp$ladder, collapse = ", "), tp$alpha_linear_gap, verdict(tp$alpha_linear_gap <= 1e-12)))
for (i in seq_len(nrow(tp$rungs))) {
  rg <- tp$rungs[i, ]
  cat(sprintf("    b = %.1f, log sigma^2 grid: R's kernel vs its powered target TV %.0e, vs the likelihood path %.3f; D09's kernel vs the likelihood path %.0e\n",
              rg$b, rg$tv_r_vs_powered, rg$tv_r_vs_likelihood_path, rg$tv_corrected_vs_likelihood_path))
}
for (i in seq_len(nrow(tp$swaps))) {
  sp <- tp$swaps[i, ]
  cat(sprintf("    swap %.1f <-> %.1f, largest gap from detailed balance over all state pairs: R vs its own kernels %.2f (= (b1 - b2) x the grid's log sigma^2 range %.2f), vs the likelihood path %.0f; D09's swap vs its kernels %.0e\n",
              sp$b1, sp$b2, sp$r_vs_r_kernels, (sp$b1 - sp$b2) * tp$log_sigma2_range, sp$r_vs_likelihood_path,
              sp$corrected_vs_corrected_kernels))
}
dc <- tp$discrepancy_covariance
cat(sprintf("    discrepancy draw at b = %s: R's covariance vs the powered conditional %s, vs the likelihood path's %s; R's swap misses its own target's discrepancy prior by up to %.2f\n",
            paste(dc$b, collapse = ", "), paste(sprintf("%.0e", dc$r_vs_powered), collapse = ", "),
            paste(sprintf("%.3f", dc$r_vs_likelihood_path), collapse = ", "), tp$swap_vs_powered_v_prior_gap))
cat(sprintf("    R's moves are exact for its powered target and its swap leaves it (gaps >= 0.1): %s; D09's likelihood path balances (<= 1e-10): %s\n",
            verdict(max(tp$rungs$tv_r_vs_powered, dc$r_vs_powered) <= 1e-12 &&
                    min(tp$swaps$r_vs_r_kernels, tp$swap_vs_powered_v_prior_gap) >= 0.1),
            verdict(max(tp$rungs$tv_corrected_vs_likelihood_path) <= 1e-12 &&
                    max(tp$swaps$corrected_vs_corrected_kernels) <= 1e-10)))

# C-O2 (D09): dense against low-rank covariance; units.
cat("\nD09 covariance formulations and units (C-O2)\n")
cv <- run("bass_covariance_oracle.R", basis_cal)
fn <- cv$functional
for (i in seq_len(nrow(fn))) {
  cat(sprintf("  functional, %-8s (cond %.0e): eigen vs Cholesky ldet %.0e, solve %.0e; R's swm %.0e, %.0e; low rank %.0e, %.0e; backward errors <= %.0e\n",
              fn$case[i], fn$condition[i], fn$eigen_vs_reference$ldet_abs[i], fn$eigen_vs_reference$solve_rel[i],
              fn$r_woodbury$ldet_abs[i], fn$r_woodbury$solve_rel[i], fn$low_rank$ldet_abs[i], fn$low_rank$solve_rel[i],
              max(unlist(fn$backward_error[i, ]))))
}
cat(sprintf("    low rank matches dense where dense is a reference (well <= 1e-12, moderate <= 1e-8): %s; every route backward stable (<= 1e-14): %s\n",
            verdict(fn$low_rank$ldet_abs[1] <= 1e-12 && fn$low_rank$solve_rel[1] <= 1e-12 &&
                    fn$low_rank$ldet_abs[2] <= 1e-8 && fn$low_rank$solve_rel[2] <= 1e-8),
            verdict(max(unlist(fn$backward_error)) <= 1e-14)))
mv <- cv$multivariate
cat(sprintf("  multivariate (cond %.0e): eigen vs Cholesky ldet %.0e, solve %.0e: %s\n", mv$condition,
            mv$eigen_vs_reference$ldet_abs, mv$eigen_vs_reference$solve_rel,
            verdict(mv$eigen_vs_reference$ldet_abs <= 1e-10 && mv$eigen_vs_reference$solve_rel <= 1e-10)))
un <- cv$units
cat(sprintf("  units: standardized minus physical log-likelihood %.6f = sum log y_s %.6f (spread %.0e); posteriors agree to %.0e: %s; R's units: TV %.2f, posterior mean %.2f vs %.2f\n",
            un$loglik_difference_mean, un$sum_log_y_s, un$loglik_difference_spread,
            un$max_abs_posterior_physical_vs_standardized,
            verdict(un$loglik_difference_spread <= 1e-10 && abs(un$loglik_difference_mean - un$sum_log_y_s) <= 1e-10 &&
                    un$max_abs_posterior_physical_vs_standardized <= 1e-12),
            un$tv_physical_vs_r_units, un$posterior_mean[[3]], un$posterior_mean[[1]]))

# D03's tau prior: the unit-information Gamma family Gamma(a, rate a*N), Zellner-Siow (a = 1/2) against Francom et
# al. 2018's shape a = 1/N, which the correction replaces by a shape independent of N.
cat("\nD03 tau prior: Gamma(a, rate a*N), Zellner-Siow a = 1/2 against Francom 2018's a = 1/N\n")
tau_cases <- rbind(c(40, 1, 0.05), c(40, 1, 0.2), c(40, 1, 0.5), c(40, 2, 0.25),
                   c(100, 1, 0.05), c(100, 1, 0.2), c(100, 1, 0.5), c(100, 2, 0.25),
                   c(1000, 1, 0.05), c(1000, 1, 0.2), c(1000, 1, 0.5), c(1000, 2, 0.25))
tp <- run("bass_tau_prior_oracle.R", tau_cases)
zs <- tp$zs; fr <- tp$francom
moments_ok <- all(abs(c(zs$mass, fr$mass) - 1) <= 1e-10) &&
  all(abs(c(zs$log_mean - zs$log_mean_closed, fr$log_mean - fr$log_mean_closed)) <= 1e-8 * pmax(1, abs(c(zs$log_mean, fr$log_mean)))) &&
  all(abs(c(zs$log_sd - zs$log_sd_closed, fr$log_sd - fr$log_sd_closed)) <= 1e-8 * c(zs$log_sd, fr$log_sd))
cat(sprintf("  log-moments by quadrature against psi and psi': %s\n", verdict(moments_ok)))
for (N in c(40, 100, 1000)) {
  i <- which(zs$N == N)[1]
  cat(sprintf("  N = %4d: E[ln tau] - ln(1/N): ZS %6.2f, Francom %8.2f | sd(ln tau): ZS %.2f, Francom %7.2f | median: ZS %.1e, Francom e^%.0f\n",
              N, zs$log_mean[i] + log(N), fr$log_mean[i] + log(N), zs$log_sd[i], fr$log_sd[i], exp(zs$log_median[i]), fr$log_median[i]))
}
zs_offset <- zs$log_mean + log(zs$N); fr_offset <- fr$log_mean + log(fr$N)
cat(sprintf("  ZS's offset from unit information is the same at every N (%.2f, within 1.3 nats); Francom's exceeds 0.9 N nats: %s\n",
            zs_offset[1], verdict(diff(range(zs_offset)) <= 1e-10 && all(abs(zs_offset) <= 1.3) && all(abs(fr_offset) >= 0.9 * fr$N))))
cat(sprintf("  Bayes factors converged over twice the range (<= 1e-8): %s\n",
            verdict(all(abs(c(zs$ln_bf - zs$ln_bf_wide, fr$ln_bf - fr$ln_bf_wide)) <= 1e-8))))
gap <- zs$ln_bf - fr$ln_bf
null_ok <- TRUE; pair_ok <- TRUE
for (R2 in c(0.05, 0.2, 0.5)) {
  j <- which(zs$k == 1 & zs$R2 == R2)
  cat(sprintf("  leaving the null, one term at R2 = %.2f: Francom's extra penalty %s nats at N = 40, 100, 1000\n",
              R2, paste(sprintf("%.2f", gap[j]), collapse = ", ")))
  null_ok <- null_ok && all(gap[j] > 0) && all(diff(gap[j]) > 0)
}
for (N in c(40, 100, 1000)) {
  one <- which(zs$N == N & zs$k == 1 & zs$R2 == 0.2); two <- which(zs$N == N & zs$k == 2)
  pair <- gap[two] - gap[one]
  cat(sprintf("  between two structures (2 terms at R2 = 0.25 against 1 at 0.20), N = %4d: the priors differ by %.2f nats\n", N, pair))
  pair_ok <- pair_ok && abs(pair) < 0.3 && abs(pair) < gap[one]
}
cat(sprintf("  Francom's prior penalizes leaving the null, more as N grows: %s; between structures that both carry tau it nearly cancels (< 0.3 nats): %s\n",
            verdict(null_ok), verdict(pair_ok)))

# The shared family's hyper-g members (2026-10-29), a = 3: the law on u normalizes, g'/(1 + g') has its Beta(1, a/2 - 1)
# mean, and hyper-g's Bayes factors by quadrature meet Liang et al.'s closed form.
hg <- tp$hyperg; hgn <- tp$hypergn
cat("\nD03 tau prior, the shared family's hyper-g members (a = 3): hyper-g on g = 1/tau, hyper-g/n on g/N\n")
cat(sprintf("  the density of ln tau normalizes, and E[g'/(1 + g')] = 2/a, by quadrature (<= 1e-10): %s\n",
            verdict(all(abs(c(hg$mass, hgn$mass) - 1) <= 1e-10) && all(abs(c(hg$beta_mean, hgn$beta_mean) - 2 / 3) <= 1e-10))))
hg_gap <- max(abs(hg$ln_bf - hg$ln_bf_closed))
cat(sprintf("  Bayes factors converged over twice the range (<= 1e-8): %s; hyper-g's against Liang et al.'s eq. (17) (%.0e): %s\n",
            verdict(all(abs(c(hg$ln_bf - hg$ln_bf_wide, hgn$ln_bf - hgn$ln_bf_wide)) <= 1e-8)), hg_gap, verdict(hg_gap <= 1e-8)))
for (N in c(40, 1000)) {
  i <- which(hg$N == N & hg$k == 1 & hg$R2 == 0.2)
  cat(sprintf("  N = %4d, one term at R2 = 0.20: ln BF against the null, ZS %.2f, hyper-g %.2f, hyper-g/n %.2f\n",
              N, zs$ln_bf[i], hg$ln_bf[i], hgn$ln_bf[i]))
}

cat("\nD03 intercept: Liang's flat intercept against level centring (F-O9, the D03 alternative)\n")
int_cases <- rbind(c(40, 1, 0.05), c(40, 1, 0.2), c(40, 1, 0.5), c(40, 2, 0.25),
                   c(100, 1, 0.05), c(100, 1, 0.2), c(100, 1, 0.5), c(100, 2, 0.25),
                   c(1000, 1, 0.05), c(1000, 1, 0.2), c(1000, 1, 0.5), c(1000, 2, 0.25))
ic <- run("bass_intercept_oracle.R", int_cases, "grid")
cat(sprintf("  Bayes factors converged over twice the range (<= 1e-8): %s\n",
            verdict(all(abs(c(ic$liang - ic$liang_wide, ic$level - ic$level_wide)) <= 1e-8))))
credit <- ic$level - ic$liang
for (N in c(40, 100, 1000)) {
  j <- which(ic$N == N)
  cat(sprintf("  N = %4d: level centring's extra credit against the null, one term at R2 = .05/.2/.5: %s; two terms at .25: %.2f\n",
              N, paste(sprintf("%.2f", credit[j][ic$k[j] == 1]), collapse = "/"), credit[j][ic$k[j] == 2]))
}
pair_credit <- sapply(c(40, 100, 1000), function(N) {
  credit[ic$N == N & ic$k == 2] - credit[ic$N == N & ic$k == 1 & ic$R2 == 0.2]
})
cat(sprintf("  between two structures one term apart (2 at R2 = 0.25 against 1 at 0.20): %s nats more, at N = 40, 100, 1000\n",
            paste(sprintf("%.3f", pair_credit), collapse = ", ")))
cat(sprintf("  positive everywhere, 0.4-0.5 for one term, 0.2-0.3 between structures, flat in N (spread <= 0.02): %s\n",
            verdict(all(credit > 0) && all(credit[ic$k == 1] > 0.4 & credit[ic$k == 1] < 0.5) &&
                    all(pair_credit > 0.2 & pair_credit < 0.3) && diff(range(pair_credit)) <= 0.02)))

set.seed(9)
xd <- sort(runif(40)); yd <- 2 + 0.8 * sin(2 * pi * xd) + rnorm(40, sd = 0.3)
idd <- run("bass_intercept_oracle.R", cbind(xd, yd), "design")
for (i in seq_len(nrow(idd))) {
  cat(sprintf("  design, offset %6.0f: level %s | liang %s | zero %s (one hinge, two hinges)\n", idd$offset[i],
              paste(sprintf("%7.3f", idd$level[[i]]), collapse = " "), paste(sprintf("%7.3f", idd$liang[[i]]), collapse = " "),
              paste(sprintf("%8.3f", idd$zero[[i]]), collapse = " ")))
}
lev <- do.call(rbind, idd$level); lg <- do.call(rbind, idd$level_g); li <- do.call(rbind, idd$liang); ze <- do.call(rbind, idd$zero)
cat(sprintf("  the explicit design's level marginal equals the g-form identity (<= 1e-8): %s\n", verdict(max(abs(lev - lg)) <= 1e-8)))
cat(sprintf("  level centring and Liang do not move with y's offset (<= 1e-8); R's zero centre does (> 1 nat at -250): %s\n",
            verdict(max(abs(sweep(lev, 2, lev[1, ]))) <= 1e-8 && max(abs(sweep(li, 2, li[1, ]))) <= 1e-8 &&
                    max(abs(ze[3, ] - ze[1, ])) > 1)))

# F-O3 and F-O4 (D01, D02 as decided): count, order and lambda laws, and the boundary ratios, under a constant
# likelihood with the support rule and D06's centred angle gate.
cat("\nD01/D02 laws under a constant likelihood (F-O3) and the boundary ratios (F-O4)\n")
show_laws <- function(label, lw_run) {
  lw <- lw_run$laws
  cat(sprintf("  %s: %d terms, %d supported, %d states in F (%d refused by the angle gate)\n", label, lw_run$n_terms,
              lw_run$n_supported, lw_run$n_states, lw$n_gate_refused))
  cat(sprintf("    count law M = 0..%d: target %s | D01 %s | prior-cancelling %s | R's thirds %s\n", max(lw$M),
              paste(sprintf("%.4f", lw$count$target), collapse = " "), paste(sprintf("%.4f", lw$count$d01), collapse = " "),
              paste(sprintf("%.4f", lw$count$cancel), collapse = " "), paste(sprintf("%.4f", lw$count$thirds), collapse = " ")))
  laws_ok <- max(lw$tv$d01, lw$tv$cancel, lw$tv_order$d01, lw$tv_order$cancel) <= 1e-12
  cat(sprintf("    D01 and the prior-cancelling schedule give the target's state, count and order laws (TV <= 1e-12): %s; R's thirds do not (TV %.3f): %s\n",
              verdict(laws_ok), lw$tv$thirds, verdict(lw$tv$thirds >= 0.05)))
  b <- lw$boundary
  fmt <- function(v) paste(sprintf("%.6f", v[!is.na(v)]), collapse = " ")
  cat(sprintf("    boundary ratios over the target's (low, high%s): D01 %s | cancelling %s | R's thirds %s\n",
              if (max(lw$M) >= 3) ", middle" else "", fmt(b$d01), fmt(b$cancel), fmt(b$thirds)))
  exact <- c(b$d01, b$cancel); exact <- exact[!is.na(exact)]
  th <- b$thirds   # low, high, middle (NA below M_max = 3); JSON keeps the order, not the names
  thirds_ok <- abs(th[1] - 3) <= 1e-9 && abs(th[2] - 2 / 3) <= 1e-9 && (is.na(th[3]) || abs(th[3] - 1) <= 1e-9)
  cat(sprintf("    D01's and the cancelling schedule's are exact (<= 1e-9): %s; R's thirds give 3 at the bottom and 2/3 at the top%s: %s\n",
              verdict(max(abs(exact - 1)) <= 1e-9), if (max(lw$M) >= 3) ", and are right in the middle" else "", verdict(thirds_ok)))
  probs_ok <- all(lw$cancel_probs >= 0 & lw$cancel_probs <= 1) && max(abs(rowSums(lw$cancel_probs) - 1)) <= 1e-12
  cat(sprintf("    the cancelling schedule's move probabilities are proper at every M: %s\n", verdict(probs_ok)))
  share <- lw$admitted_mass / sum(lw$admitted_mass)
  tilt_ok <- max(abs(lw$order_at_one$target - share), abs(lw$order_at_one$d01 - share)) <= 1e-12
  cat(sprintf("    the tilt at M = 1: admitted prior mass per order %s, order law %s, against the unconstrained %s: %s\n",
              paste(sprintf("%.4f", lw$admitted_mass), collapse = "/"), paste(sprintf("%.4f", lw$order_at_one$d01), collapse = "/"),
              paste(sprintf("%.4f", rep(1 / length(share), length(share))), collapse = "/"),
              verdict(tilt_ok && max(abs(share - 1 / length(share))) >= 0.01)))
  lam_gap <- max(abs(lw$lambda$mixture / lw$lambda$tilted - 1))
  cat(sprintf("    lambda at retention, Gamma(h1 + M, h2 + 1) over the count law, is the prior tilted by F (%.0e): %s\n",
              lam_gap, verdict(lam_gap <= 1e-10)))
}
show_laws("p = 3, {0, 1}^3, K_max = 2, M_max = 3, npart = 2",
          run("bass_enumeration_oracle.R", cbind(0, grid8$x1, grid8$x2, grid8$x3), "laws", "M_max=3", "npart=2"))
grid16 <- expand.grid(x1 = c(0, 1, 2, 3) / 3, x2 = c(0, 1, 2, 3) / 3)
show_laws("p = 2, four levels per variable, K_max = 2, M_max = 2, npart = 3",
          run("bass_enumeration_oracle.R", cbind(0, grid16$x1, grid16$x2), "laws", "M_max=2", "npart=3"))
