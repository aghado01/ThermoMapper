# bass_stage4_checks.R — reproduces the BASS port's stage-4 oracle checks (the fitting segment's stage-4 plan; chip
# brief issues/regression/bass/stage4-oracle-modes-brief.md): the exact posterior over an enumerable structure space
# with data (S6), its evidence and tempered laws (S8), and R BASS's own sufficient statistics (S4, S5), by running the
# oracles on fixed fixtures, exactly as the C# bridge would. Exits with status 1 if any check fails.
#
# Usage (cwd = r/):  Rscript scripts/bass_stage4_checks.R

rscript <- file.path(R.home("bin"), "Rscript")
write_fixture <- function(fixture) {
  csv <- tempfile(fileext = ".csv")
  write.table(fixture, csv, sep = ",", row.names = FALSE, col.names = FALSE)
  csv
}
run_files <- function(oracle, csv, ...) {
  out <- tempfile(fileext = ".json")
  status <- system2(rscript, c(file.path("oracles", oracle), csv, out, ...))
  if (status != 0) stop(oracle, " failed")
  # Numbers simplified; the term and state lists, whose shapes vary, as parsed.
  res <- jsonlite::read_json(out, simplifyVector = TRUE)
  raw <- jsonlite::read_json(out, simplifyVector = FALSE)
  if (!is.null(raw$posterior)) {
    res$terms <- lapply(raw$terms, function(t) list(cont = matrix(as.numeric(unlist(t$cont)), ncol = 3, byrow = TRUE),
                                                    cat = lapply(t$cat, function(f) list(var = f[[1]], levels = unlist(f[[2]])))))
    res$states <- lapply(raw$states, function(s) as.integer(unlist(s)))
    res$posterior$tempered <- raw$posterior$tempered
  }
  res
}
run <- function(oracle, fixture, ...) run_files(oracle, write_fixture(fixture), ...)
failures <- 0
report <- function(text, ok) {
  if (!isTRUE(ok)) failures <<- failures + 1
  cat(sprintf("    %s: %s\n", text, if (isTRUE(ok)) "PASS" else "FAIL"))
}
hinge <- function(x, t, s) pmax(s * (x - t), 0)
lse <- function(l) { m <- max(l); m + log(sum(exp(l - m))) }
tv <- function(a, b) 0.5 * sum(abs(a - b))
num <- function(x) if (is.null(x)) NA_real_ else x
rung_field <- function(res, form, field) vapply(res$posterior$tempered[[form]], function(r) num(r[[field]]), 0)
rung_law <- function(res, form, b) {
  r <- res$posterior$tempered[[form]][[match(b, unlist(res$posterior$tempered$rungs))]]
  unlist(r$law)
}

# A state's columns from its terms, read in the C# conventions the posterior mode emits (0-based variables, knots
# indexing each variable's ascending distinct values, levels indexing each factor's), as plain hinges times
# indicators: the marginal reads the column space only.
columns_of <- function(res, S) {
  X <- res$posterior$data$x; Z <- res$posterior$data$z
  knots <- lapply(seq_len(ncol(X)), function(v) sort(unique(X[, v])))
  B <- matrix(1, nrow(X), length(S))
  for (j in seq_along(S)) {
    cont <- res$terms[[S[j] + 1]]$cont
    for (f in seq_len(nrow(cont))) B[, j] <- B[, j] * hinge(X[, cont[f, 1] + 1], knots[[cont[f, 1] + 1]][cont[f, 3] + 1], cont[f, 2])
    for (f in res$terms[[S[j] + 1]]$cat) {
      levels <- sort(unique(Z[, f$var + 1]))
      B[, j] <- B[, j] * (Z[, f$var + 1] %in% levels[f$levels + 1])
    }
  }
  B
}
# D01's term measure, recomputed from the emitted terms: J uniform on {0..K_max} beside a categorical part (on
# {1..K_max} without one), L uniform on {0..L_max}, variable sets uniform, signs and knots uniform, and each subset
# by its measure; the state's prior Gamma(h1 + M)/(h2 + 1)^M prod pi(t), normalized over the emitted F.
log_prior_of <- function(res, K_max, h1 = 10, h2 = 10) {
  X <- res$posterior$data$x; Z <- res$posterior$data$z; cg <- res$posterior$categorical
  n_v <- vapply(seq_len(ncol(X)), function(v) length(unique(X[, v])), 0)
  mixed <- is.list(cg)
  lt <- vapply(res$terms, function(t) {
    J <- nrow(t$cont); L <- length(t$cat)
    v <- -(if (mixed) log(K_max + 1) else log(K_max)) - lchoose(ncol(X), J) - J * log(2) - sum(log(n_v[t$cont[, 1] + 1]))
    if (mixed) {
      n_l <- vapply(seq_len(ncol(Z)), function(u) length(unique(Z[, u])), 0)
      v <- v - log(cg$L_max + 1) - lchoose(ncol(Z), L)
      for (f in t$cat) {
        n <- n_l[f$var + 1]
        v <- v + if (cg$subsets == "uniform") -log(2^n - 2) else -log(n - 1) - lchoose(n, length(f$levels))
      }
    }
    v
  }, 0)
  ls <- vapply(res$states, function(S) lgamma(h1 + length(S)) - length(S) * log(h2 + 1) + sum(lt[S + 1]), 0)
  ls - lse(ls)
}
map_of_size <- function(res, M) { i <- which(res$posterior$M == M); i[which.max(res$posterior$posterior[i])] }

# Fixtures.
#   A: three values per input, the 9-point grid three times (N = 27), y = 1 + 2 x1 + 3 (x1 - 1/2)+ (1/2 - x2)+ plus
#      N(0, 0.3^2) noise at seed 11: S6's chain fixture at M_max = 3, and the options' fixture at M_max = 2.
#   B: four values per input, the 16-point grid twice (N = 32), y = 1 + 2 (x1 - 1/3)+ + 3 (x1 - 1/3)+ (2/3 - x2)+ plus
#      N(0, 0.25^2) at seed 7, at M_max = 2: the laws mode's comparison and the ladder.
#   E: fixture A's inputs with y = 1 + 2 (x1 - 1/2)+ exactly, which one hinge and the intercept fit.
g3 <- expand.grid(x1 = c(0, 0.5, 1), x2 = c(0, 0.5, 1))[rep(1:9, 3), ]
set.seed(11)
FA <- cbind(1 + 2 * g3$x1 + 3 * hinge(g3$x1, 0.5, 1) * hinge(g3$x2, 0.5, -1) + rnorm(27, sd = 0.3), g3$x1, g3$x2)
g4 <- expand.grid(x1 = (0:3) / 3, x2 = (0:3) / 3); g4 <- rbind(g4, g4)
set.seed(7)
FB <- cbind(1 + 2 * hinge(g4$x1, 1 / 3, 1) + 3 * hinge(g4$x1, 1 / 3, 1) * hinge(g4$x2, 2 / 3, -1) + rnorm(32, sd = 0.25),
            g4$x1, g4$x2)
FE <- cbind(1 + 2 * hinge(g3$x1, 0.5, 1), g3$x1, g3$x2)

# ---- Task A: the posterior (S6) ---------------------------------------------------------------------------------
cat("The exact posterior over F with data (S6)\n")
A <- run("bass_enumeration_oracle.R", FA, "posterior", "M_max=3"); P <- A$posterior
cat(sprintf("  fixture A, M_max = 3: %d terms, %d supported, %d states in F (%d refused by the gate)\n", A$n_terms,
            A$n_supported, A$n_states, P$n_gate_refused))
cat(sprintf("    count law M = 0..3: prior %s | posterior %s; order law (n1, n2): %s\n",
            paste(sprintf("%.4f", P$count$prior), collapse = " "), paste(sprintf("%.4f", P$count$posterior), collapse = " "),
            paste(sprintf("(%s) %.4f", apply(P$order$levels, 1, paste, collapse = ","), P$order$posterior), collapse = " ")))
cat(sprintf("    log Z = %.8f; tau: mean %.6f, E[log tau] %.4f, quantiles (2.5, 50, 97.5%%) %.6f %.6f %.6f\n", P$log_Z,
            P$tau$mean, P$tau$log_mean, P$tau$quantiles[1], P$tau$quantiles[4], P$tau$quantiles[7]))
top <- order(P$posterior, decreasing = TRUE)[1:3]
cat(sprintf("    top states: %s\n", paste(sprintf("{%s} %.4f", vapply(A$states[top], paste, "", collapse = ","),
                                                  P$posterior[top]), collapse = "  ")))
incl <- numeric(A$n_terms)
for (i in seq_len(A$n_states)) incl[A$states[[i]] + 1] <- incl[A$states[[i]] + 1] + P$posterior[i]
gap <- max(abs(sum(P$posterior) - 1), abs(sum(exp(P$log_prior)) - 1), abs(sum(P$count$posterior) - 1),
           abs(sum(P$order$posterior) - 1), abs(incl - P$inclusion))
report(sprintf("the posterior and prior normalize, and the count, order and inclusion laws are views of it (%.0e)", gap), gap <= 1e-12)
gap <- abs(lse(P$log_prior + P$log_marginal) - P$log_Z)
report(sprintf("log Z equals the direct sum over the states (%.0e)", gap), gap <= 1e-12)
gap <- max(abs(log_prior_of(A, 2) - P$log_prior))
report(sprintf("the prior equals D01's measure and D02's count prior recomputed from the emitted terms (%.0e)", gap), gap <= 1e-12)
report(sprintf("the closed forms against quadrature over sigma^2 (%.0e) and the tau rule against integrate() (%.0e)",
               P$quadrature$sigma2_gap, P$quadrature$tau_gap), max(P$quadrature$sigma2_gap, P$quadrature$tau_gap) <= 1e-9)

# Each option on fixture A at M_max = 2.
cat("  the options on fixture A, M_max = 2\n")
A2 <- list()
for (o in list(c(name = "mean (default)"), c(name = "flat", "centre=flat"), c(name = "zero", "centre=zero"),
               c(name = "jeffreys", "sigma2_prior=jeffreys"), c(name = "flat + jeffreys", "centre=flat", "sigma2_prior=jeffreys"),
               c(name = "fixed g", "coef=fixed"))) {
  r <- run("bass_enumeration_oracle.R", FA, "posterior", "M_max=2", unname(o[-1]))
  A2[[o[["name"]]]] <- r
  q <- r$posterior$quadrature
  report(sprintf("%-15s normalizes (%.0e); closed forms against sigma^2 quadrature %.0e, tau rule %.0e; log Z %s", o[["name"]],
                 abs(sum(r$posterior$posterior) - 1), q$sigma2_gap, q$tau_gap,
                 if (is.null(r$posterior$log_Z)) "null (improper)" else sprintf("%.6f", r$posterior$log_Z)),
         abs(sum(r$posterior$posterior) - 1) <= 1e-12 && max(q$sigma2_gap, q$tau_gap) <= 1e-9 &&
           is.null(r$posterior$log_Z) == (o[["name"]] %in% c("flat", "jeffreys", "flat + jeffreys")))
}

# One structure's marginal against bass_marginal_oracle.R's quadrature, with the constants it leaves out restored:
# log p = C - (n/2) log 2 pi [- log(N)/2 under flat] + lgamma(g1 + n/2) + its log integral, n = N (N - 1 under flat).
cat("  marginals against bass_marginal_oracle.R (tau by integrate() on log tau, the lower limit at 1e-64)\n")
for (name in c("mean (default)", "flat", "zero", "jeffreys")) {
  r <- A2[[name]]; P <- r$posterior; N <- r$n
  flat <- P$response$centre == "flat"
  g1 <- P$response$sigma2_shape; g2 <- P$response$sigma2_rate
  gaps <- vapply(1:2, function(M) {
    i <- map_of_size(r, M); B <- columns_of(r, r$states[[i]])
    y <- r$posterior$data$y
    fixture <- if (flat) cbind(y, B) else cbind(if (P$response$centre == "zero") y else y - mean(y), 1, B)
    m <- run("bass_marginal_oracle.R", fixture, if (flat) "liang" else "r", g1, g2)
    n_eff <- N - flat
    constant <- (if (g1 > 0) g1 * log(g2) - lgamma(g1) else 0) - n_eff / 2 * log(2 * pi) - flat * 0.5 * log(N) + lgamma(g1 + n_eff / 2)
    abs(P$log_marginal[i] - (constant + m$log_integral[5]))
  }, 0)
  report(sprintf("%-15s MAP states at M = 1, 2: %.0e %.0e", name, gaps[1], gaps[2]), max(gaps) <= 1e-7)
}

# Bayes factors against the intercept-only structure under the Jeffreys sigma^2, against bass_intercept_oracle.R's
# grid mode: Liang (flat) and level centring, each from (N, M, R^2) with R^2 = Q / (rss + Q).
cat("  Bayes factors against bass_intercept_oracle.R (Liang and level centring, Jeffreys sigma^2)\n")
lj <- A2[["flat + jeffreys"]]$posterior; mj <- A2[["jeffreys"]]$posterior
pick <- unlist(lapply(1:2, function(M) { i <- which(mj$M == M); i[order(mj$posterior[i], decreasing = TRUE)[1:3]] }))
r2 <- mj$Q[pick] / (mj$rss[pick] + mj$Q[pick])
ib <- run("bass_intercept_oracle.R", cbind(A2[["jeffreys"]]$n, mj$M[pick], r2), "grid")
null_l <- lj$log_marginal[lj$M == 0]; null_m <- mj$log_marginal[mj$M == 0]
gap_l <- max(abs(lj$log_marginal[pick] - null_l - ib$liang)); gap_m <- max(abs(mj$log_marginal[pick] - null_m - ib$level))
report(sprintf("six states at M = 1, 2: Liang %.0e, level %.0e (Bayes factors %s nats)", gap_l, gap_m,
               paste(sprintf("%.2f", ib$level), collapse = " ")), max(gap_l, gap_m) <= 1e-7)
same <- max(abs(lj$rss - mj$rss), abs(lj$Q - mj$Q))
report(sprintf("flat and level centring read the same residual and fit energy (%.0e)", same), same <= 1e-10)

# Shift and scale (D03's independent checks).
cat("  invariances on fixture A, M_max = 2\n")
base <- A2[["mean (default)"]]$posterior
shifted <- run("bass_enumeration_oracle.R", cbind(FA[, 1] + 10, FA[, -1]), "posterior", "M_max=2")$posterior
shifted_flat <- run("bass_enumeration_oracle.R", cbind(FA[, 1] + 10, FA[, -1]), "posterior", "M_max=2", "centre=flat")$posterior
shifted_zero <- run("bass_enumeration_oracle.R", cbind(FA[, 1] + 10, FA[, -1]), "posterior", "M_max=2", "centre=zero")$posterior
scaled <- run("bass_enumeration_oracle.R", cbind(FA[, 1] * 1e3, FA[, -1]), "posterior", "M_max=2")$posterior
report(sprintf("y + 10: level centring and Liang do not move (%.0e, %.0e), R's zero centre does (TV %.3f)",
               max(abs(shifted$posterior - base$posterior), abs(shifted$log_Z - base$log_Z)),
               max(abs(shifted_flat$posterior - A2[["flat"]]$posterior$posterior)),
               tv(shifted_zero$posterior, A2[["zero"]]$posterior$posterior)),
       max(abs(shifted$posterior - base$posterior), abs(shifted$log_Z - base$log_Z),
           abs(shifted_flat$posterior - A2[["flat"]]$posterior$posterior)) <= 1e-10 &&
         tv(shifted_zero$posterior, A2[["zero"]]$posterior$posterior) >= 0.01)
report(sprintf("y x 1000 under the scaled sigma^2: the posterior does not move (%.0e), log Z moves by -N log 1000 (%.0e)",
               max(abs(scaled$posterior - base$posterior)), abs(scaled$log_Z - base$log_Z + 27 * log(1e3))),
       max(abs(scaled$posterior - base$posterior)) <= 1e-10 && abs(scaled$log_Z - base$log_Z + 27 * log(1e3)) <= 1e-8)

# The exact-fit exclusion (D03): under the Jeffreys sigma^2 F drops the structures that fit y to working precision;
# under the default's proper sigma^2 they stay, with finite marginals.
cat("  the exact-fit exclusion on fixture E, M_max = 2\n")
Ed <- run("bass_enumeration_oracle.R", FE, "posterior", "M_max=2"); Ej <- run("bass_enumeration_oracle.R", FE, "posterior", "M_max=2", "sigma2_prior=jeffreys")
yt <- FE[, 1] - mean(FE[, 1])
exact <- sqrt(Ed$posterior$rss) <= 1e3 * .Machine$double.eps * sqrt(sum(yt^2))
cat(sprintf("    default: %d states, %d exact fits holding %.4f of the posterior; Jeffreys: %d states, %d refused\n", Ed$n_states,
            sum(exact), sum(Ed$posterior$posterior[exact]), Ej$n_states, Ej$posterior$n_exact_fit_refused))
report("the Jeffreys option refuses exactly the default's exact fits, and both posteriors normalize",
       sum(exact) > 0 && Ej$posterior$n_exact_fit_refused == sum(exact) && Ej$n_states == Ed$n_states - sum(exact) &&
         all(is.finite(Ed$posterior$log_marginal)) && abs(sum(Ed$posterior$posterior) - 1) <= 1e-12 &&
         abs(sum(Ej$posterior$posterior) - 1) <= 1e-12)

# Under a constant likelihood the posterior is the laws mode's target (F-O3), and tau's is its prior.
cat("  a constant likelihood on fixture B's inputs, M_max = 2, against the laws mode\n")
F0 <- cbind(0, FB[, -1])
laws <- run("bass_enumeration_oracle.R", F0, "laws", "M_max=2")$laws
flat_like <- run("bass_enumeration_oracle.R", F0, "posterior", "M_max=2", "likelihood=constant"); C <- flat_like$posterior
gap <- max(abs(C$count$posterior - laws$count$target), abs(C$order$posterior - laws$order$target), abs(C$posterior - exp(C$log_prior)))
report(sprintf("count and order laws equal the laws mode's target (%.0e), on %d states", gap, flat_like$n_states), gap <= 1e-12)
q_tau <- qgamma(C$tau$levels, 0.5, rate = flat_like$n / 2)
gap <- max(abs(C$tau$quantiles / q_tau - 1), abs(C$tau$mean * flat_like$n - 1))
report(sprintf("tau's law is its prior: mean 1/N and qgamma's quantiles (%.0e relative)", gap), gap <= 1e-10)

# ---- Task E: a categorical predictor (S6's mixed chain) --------------------------------------------------------
cat("\nA categorical predictor beside the continuous ones (S6's mixed chain)\n")
gc <- expand.grid(x1 = c(0, 0.5, 1), x2 = c(0, 0.5, 1), z = 0:2)
set.seed(21)
FC <- cbind(1 + 1.5 * gc$x1 + 1.5 * hinge(gc$x2, 0.5, 1) * (gc$z == 0) + rnorm(nrow(gc), sd = 0.4), gc$x1, gc$x2, gc$z)
Cm <- run("bass_enumeration_oracle.R", FC, "posterior", "M_max=2", "K_max=1", "categorical=1"); P <- Cm$posterior
cat(sprintf("  fixture C (x1, x2 on three values, a three-level factor, N = %d, K_max = 1, L_max = 1, M_max = 2): %d terms, %d supported, %d states in F (%d refused by the gate)\n",
            Cm$n, Cm$n_terms, Cm$n_supported, Cm$n_states, P$n_gate_refused))
type_mass <- vapply(seq_len(nrow(P$order$types)), function(j) sum(P$order$posterior[P$order$levels[, j] > 0]), 0)
cat(sprintf("    count law M = 0..2: prior %s | posterior %s; posterior mass with a term of type (J, L): %s; log Z = %.6f\n",
            paste(sprintf("%.4f", P$count$prior), collapse = " "), paste(sprintf("%.4f", P$count$posterior), collapse = " "),
            paste(sprintf("(%s) %.4f", apply(P$order$types, 1, paste, collapse = ","), type_mass), collapse = " "), P$log_Z))
gap <- max(abs(sum(P$posterior) - 1), P$quadrature$sigma2_gap, P$quadrature$tau_gap)
report(sprintf("normalizes, and the closed forms hold against quadrature (%.0e)", gap), gap <= 1e-9)
gap <- max(abs(log_prior_of(Cm, 1) - P$log_prior))
report(sprintf("the prior equals D01's measure recomputed from the emitted terms, the subsets size-uniform (%.0e)", gap), gap <= 1e-12)
nonzero <- colSums(columns_of(Cm, seq_len(Cm$n_terms) - 1) != 0)
in_F <- sort(unique(unlist(Cm$states)))
report(sprintf("support counts the continuous part times the categorical part: terms in F have >= 3 nonzero values, the %d refused of %d fewer",
               sum(!P$supported), Cm$n_terms),
       all(nonzero[in_F + 1] >= 3) && all((nonzero >= 3) == P$supported))
mains <- which(vapply(Cm$terms, function(t) nrow(t$cont) == 0 && length(t$cat) == 1, TRUE)) - 1
levels_of <- function(t) Cm$terms[[t + 1]]$cat[[1]]$levels
single <- function(t) which(vapply(Cm$states, function(S) length(S) == 1 && S == t, TRUE))
pair_gap <- max(vapply(mains, function(t) {
  u <- mains[vapply(mains, function(v) setequal(levels_of(v), setdiff(0:2, levels_of(t))), TRUE)]
  abs(P$log_marginal[single(t)] - P$log_marginal[single(u)])
}, 0))
both <- any(vapply(Cm$states, function(S) length(S) == 2 && all(S %in% mains) &&
                     setequal(c(levels_of(S[1]), levels_of(S[2])), 0:2) && !length(intersect(levels_of(S[1]), levels_of(S[2]))), TRUE))
report(sprintf("a main effect on D and one on its complement define one model: equal marginals (%.0e), and F holds no state with both",
               pair_gap), pair_gap <= 1e-10 && !both)
has_cat <- vapply(Cm$states, function(S) any(vapply(Cm$terms[S + 1], function(t) length(t$cat) > 0, TRUE)), TRUE)
i <- which(has_cat)[which.max(P$posterior[has_cat])]
m <- run("bass_marginal_oracle.R", cbind(P$data$y - mean(P$data$y), 1, columns_of(Cm, Cm$states[[i]])), "r",
         P$response$sigma2_shape, P$response$sigma2_rate)
g1 <- P$response$sigma2_shape; g2 <- P$response$sigma2_rate
gap <- abs(P$log_marginal[i] - (g1 * log(g2) - lgamma(g1) - Cm$n / 2 * log(2 * pi) + lgamma(g1 + Cm$n / 2) + m$log_integral[5]))
report(sprintf("the most probable state with a categorical part, rebuilt from its terms, against bass_marginal_oracle.R (%.0e)", gap), gap <= 1e-7)

cat("  the subset measures under a constant likelihood, M_max = 1\n")
g4l <- expand.grid(x1 = c(0, 0.5, 1), z = 0:3)[rep(1:12, 2), ]
F4 <- cbind(0, g4l$x1, g4l$z)
s4 <- run("bass_enumeration_oracle.R", F4, "posterior", "M_max=1", "likelihood=constant", "categorical=1")
u4 <- run("bass_enumeration_oracle.R", F4, "posterior", "M_max=1", "likelihood=constant", "categorical=1", "subsets=uniform")
s3 <- run("bass_enumeration_oracle.R", cbind(0, FC[, -1]), "posterior", "M_max=1", "K_max=1", "likelihood=constant", "categorical=1")
u3 <- run("bass_enumeration_oracle.R", cbind(0, FC[, -1]), "posterior", "M_max=1", "K_max=1", "likelihood=constant", "categorical=1", "subsets=uniform")
gap <- max(abs(log_prior_of(s4, 1) - s4$posterior$log_prior), abs(log_prior_of(u4, 1) - u4$posterior$log_prior),
           abs(log_prior_of(s3, 1) - s3$posterior$log_prior), abs(log_prior_of(u3, 1) - u3$posterior$log_prior))
report(sprintf("both measures equal D01's recomputed from the terms, at three and four levels (%.0e)", gap), gap <= 1e-12)
report(sprintf("they agree at three levels (%.0e) and differ at four (largest log-prior gap %.3f)",
               max(abs(s3$posterior$log_prior - u3$posterior$log_prior)), max(abs(s4$posterior$log_prior - u4$posterior$log_prior))),
       max(abs(s3$posterior$log_prior - u3$posterior$log_prior)) <= 1e-12 && max(abs(s4$posterior$log_prior - u4$posterior$log_prior)) >= 0.01)

# ---- Task C: the evidence and tempered laws (S8) -----------------------------------------------------------------
cat("\nThe evidence and tempered laws on F (S8)\n")
h <- 1e-4
rungs <- c(0, 0.1 - h, 0.1, 0.1 + h, 0.5 - h, 0.5, 0.5 + h, 1)
B <- run("bass_enumeration_oracle.R", FB, "posterior", "M_max=2", paste0("rungs=", paste(rungs, collapse = ",")))
P <- B$posterior
cat(sprintf("  fixture B, M_max = 2: %d states; log Z = %.8f\n", B$n_states, P$log_Z))
for (form in c("integrated", "linear")) {
  lz <- rung_field(B, form, "log_Z"); ms <- rung_field(B, form, "mean_slope")
  cat(sprintf("    %-10s log Z_b at b = 0, 0.1, 0.5, 1: %s; E_b[slope] %s\n", form,
              paste(sprintf("%.6f", lz[c(1, 3, 6, 8)]), collapse = " "), paste(sprintf("%.4f", ms[c(1, 3, 6, 8)]), collapse = " ")))
  ends <- max(abs(rung_law(B, form, 1) - P$posterior), abs(rung_law(B, form, 0) - exp(P$log_prior)))
  report(sprintf("%s: b = 1 gives the posterior and b = 0 the prior on F (%.0e); log Z_1 = log Z, log Z_0 = 0 (%.0e)", form,
                 ends, max(abs(lz[8] - P$log_Z), abs(lz[1]))), ends <= 1e-12 && max(abs(lz[8] - P$log_Z), abs(lz[1])) <= 1e-12)
  fd <- c((lz[4] - lz[2]) / (2 * h), (lz[7] - lz[5]) / (2 * h))
  gap <- max(abs(fd - ms[c(3, 6)]) / abs(ms[c(3, 6)]))
  report(sprintf("%s: d log Z_b/db by central differences equals E_b[slope] at b = 0.1, 0.5 (%.0e relative)", form, gap), gap <= 1e-6)
}
report(sprintf("b = 0's law is the laws mode's target on the same F (count %.0e)", max(abs(laws$count$target - unlist(B$posterior$tempered$integrated[[1]]$count)))),
       max(abs(laws$count$target - unlist(B$posterior$tempered$integrated[[1]]$count))) <= 1e-12)
report(sprintf("the two forms differ between the ends (TV at b = 0.5: %.3f)", tv(rung_law(B, "integrated", 0.5), rung_law(B, "linear", 0.5))),
       tv(rung_law(B, "integrated", 0.5), rung_law(B, "linear", 0.5)) >= 0.01)
fixed0 <- A2[["fixed g"]]$posterior$tempered$integrated[[1]]$mean_slope
cat(sprintf("    the integrated form's slope at b = 0: null under Zellner-Siow (E[1/tau] is infinite: %s); %.1f at fixed g, about (N/2) psi(c) = %.1f\n",
            if (is.null(P$tempered$integrated[[1]]$mean_slope)) "null" else "a number", fixed0, 27 / 2 * digamma(1e-3)))

# ---- Task B: R's sufficient statistics (S4, S5) ------------------------------------------------------------------
cat("\nR BASS's sufficient statistics (S4, S5)\n")
z <- A2[["zero"]]; i <- map_of_size(z, 2); Bz <- cbind(1, columns_of(z, z$states[[i]]))
st <- run("bass_statistics_oracle.R", cbind(z$posterior$data$y, Bz))
report(sprintf("scalar, fixture A's MAP state at M = 2: R's X'X and X'y against the explicit design (%.0e, %.0e), qf against a direct solve (%.0e)",
               st$check$xtx_rel_error, st$check$xty_rel_error, st$check$qf_rel_error),
       max(st$check$xtx_rel_error, st$check$xty_rel_error, st$check$qf_rel_error) <= 1e-12 && !st$getqf_refuses)
gap <- max(abs(st$qf - z$posterior$Q[i]) / st$qf, abs(st$yty - st$qf - z$posterior$rss[i]) / st$yty)
report(sprintf("R's qf and y'y - qf equal the zero-centre posterior's Q and rss for that state (%.0e)", gap), gap <= 1e-12)
report(sprintf("log|X'X| = 2 sum log diag R (%.4f); R's s2.rate at tau = 1/N against the residual's (%.0e)", st$log_det,
               abs(st$s2_rate / st$s2_rate_direct - 1)),
       abs(st$log_det - c(determinant(crossprod(Bz))$modulus)) <= 1e-10 && abs(st$s2_rate / st$s2_rate_direct - 1) <= 1e-10)

set.seed(5)
nx <- 15; m <- 11
xa <- runif(nx); xb <- runif(nx); tt <- seq(0, 1, length.out = m)
Af <- cbind(1, hinge(xa, 0.4, 1), 1, hinge(xb, 0.6, -1) * hinge(xa, 0.2, 1), hinge(xb, 0.3, 1))   # a functional-only term third
Hf <- cbind(1, 1, hinge(tt, 0.5, 1), hinge(tt, 0.2, 1), 1)                                        # design-only terms second and last
Yf <- outer(xa, tt, function(a, b) 1 + a * b) + 0.5 * xb + matrix(rnorm(nx * m, sd = 0.05), nx)
sf <- run_files("bass_statistics_oracle.R", write_fixture(Af), "functional", write_fixture(Hf), write_fixture(Yf))
report(sprintf("functional (n = %d, m = %d, k = %d): the Hadamard X'X and a'Yh against the explicit nm x k design (%.0e, %.0e), qf (%.0e)",
               sf$n, sf$m, sf$k, sf$check$xtx_rel_error, sf$check$xty_rel_error, sf$check$qf_rel_error),
       max(sf$check$xtx_rel_error, sf$check$xty_rel_error, sf$check$qf_rel_error) <= 1e-12 && sf$N == nx * m)

n <- 2000
xl <- seq(0, 1, length.out = n)
local <- numeric(n); local[seq(n / 2 - 9, n / 2 + 10)] <- 1e-3            # 20 nonzero points, small values
sr <- run("bass_statistics_oracle.R", cbind(1 + xl + rnorm(n, sd = 0.1), 1, xl, local))
report(sprintf("R's getQf refuses a small local term at n = %d (ratio %.0f, at birth %d), flagged; the statistics stand (qf %.0e)",
               n, sr$getqf_ratio, sr$getqf_refuses_at, sr$check$qf_rel_error),
       sr$getqf_refuses && sr$getqf_refuses_at == 2 && sr$getqf_ratio > 1e3 && sr$check$qf_rel_error <= 1e-12)

cat(sprintf("\n%s\n", if (failures == 0) "stage-4 checks: all PASS" else sprintf("stage-4 checks: %d FAIL", failures)))
if (failures > 0) quit(status = 1)
