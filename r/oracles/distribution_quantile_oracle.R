# distribution_quantile_oracle.R — reference for the distributions' quantiles and truncated moments, in base R alone.
#
# Usage (cwd = r/):  Rscript oracles/distribution_quantile_oracle.R <cases.csv> <out.json> quantile
#                    Rscript oracles/distribution_quantile_oracle.R <cases.csv> <out.json> moments
#
# quantile: a row is family, a, b, log_tail, lower, x. family: 1 gamma(shape a, rate 1), 2 beta(a, b), 3 t(a), 4 normal.
#   lower: 1 for ln P(X <= x) = log_tail, 0 for ln P(X > x) = log_tail. x: the quantile under test.
#   Emits { q, r_log_tail, x_log_tail, warned }: R's quantile; R's own log tail back at it, whose distance from log_tail says
#   whether R itself solved the case (qt in the far tail stops short); R's log tail at the x under test, a round trip through
#   tails that share nothing with it; and whether R warned. An oracle is a target only where it is right.
#
# moments: a row is family, a, b, lower, upper. family: 1 normal(mean a, sd b), 2 gamma(shape a, rate b), 3 inverse gamma
#   (shape a > 2, scale b). lower and upper may be infinite.
#   Emits { log_mass, mean, variance }, each a closed form in R's log tails, which never reads a truncated draw:
#     normal: mean + sd (phi(alpha) - phi(beta)) / M and sd^2 (1 + (alpha phi(alpha) - beta phi(beta)) / M - ((phi(alpha) - phi(beta)) / M)^2);
#     gamma: E[X^k; l, u] = (Gamma(a + k) / Gamma(a)) b^-k M(a + k), with M(a) = P(l < Gamma(a, b) <= u);
#     inverse gamma: E[X^k; l, u] = (Gamma(a - k) / Gamma(a)) b^k M(a - k) on the gamma's interval b/u, b/l.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: distribution_quantile_oracle.R <cases.csv> <out.json> quantile|moments")
cases <- read_matrix(args[1])

# ln(e^hi - e^lo) for hi >= lo.
log_diff <- function(hi, lo) if (is.infinite(lo)) hi else hi + log(-expm1(lo - hi))

# ln P(l < X <= u) from the tail that holds it.
log_interval <- function(p, l, u) {
  fl <- p(l, TRUE); su <- p(u, FALSE); sl <- p(l, FALSE); fu <- p(u, TRUE)
  if (fl >= log(0.5)) return(log_diff(sl, su))
  if (fu <= log(0.5)) return(log_diff(fu, fl))
  log1p(-(exp(fl) + exp(su)))
}

quantiles <- function() {
  n <- nrow(cases)
  q <- back <- at_x <- numeric(n)
  warned <- logical(n)
  for (i in seq_len(n)) {
    family <- cases[i, 1]; a <- cases[i, 2]; b <- cases[i, 3]; tail <- cases[i, 4]; lower <- cases[i, 5] == 1
    solve <- function() switch(family,
      qgamma(tail, a, 1, lower.tail = lower, log.p = TRUE),
      qbeta(tail, a, b, lower.tail = lower, log.p = TRUE),
      qt(tail, a, lower.tail = lower, log.p = TRUE),
      qnorm(tail, lower.tail = lower, log.p = TRUE))
    log_tail_at <- function(x) suppressWarnings(switch(family,
      pgamma(x, a, 1, lower.tail = lower, log.p = TRUE),
      pbeta(x, a, b, lower.tail = lower, log.p = TRUE),
      pt(x, a, lower.tail = lower, log.p = TRUE),
      pnorm(x, lower.tail = lower, log.p = TRUE)))
    q[i] <- withCallingHandlers(solve(), warning = function(w) { warned[i] <<- TRUE; invokeRestart("muffleWarning") })
    back[i] <- log_tail_at(q[i])
    at_x[i] <- log_tail_at(cases[i, 6])
  }
  list(q = q, r_log_tail = back, x_log_tail = at_x, warned = warned)
}

moments <- function() {
  n <- nrow(cases)
  log_mass <- mean <- variance <- numeric(n)
  for (i in seq_len(n)) {
    family <- cases[i, 1]; a <- cases[i, 2]; b <- cases[i, 3]; l <- cases[i, 4]; u <- cases[i, 5]
    if (family == 1) {
      alpha <- (l - a) / b; beta <- (u - a) / b
      log_m <- log_interval(function(z, low) pnorm(z, lower.tail = low, log.p = TRUE), alpha, beta)
      ratio <- function(z) if (is.finite(z)) exp(dnorm(z, log = TRUE) - log_m) else 0
      weighted <- function(z) if (is.finite(z)) z * exp(dnorm(z, log = TRUE) - log_m) else 0
      shift <- ratio(alpha) - ratio(beta)
      log_mass[i] <- log_m
      mean[i] <- a + b * shift
      variance[i] <- b^2 * (1 + weighted(alpha) - weighted(beta) - shift^2)
    } else if (family == 2) {
      m <- function(s) log_interval(function(x, low) pgamma(x, s, b, lower.tail = low, log.p = TRUE), l, u)
      log_m <- m(a)
      m1 <- exp(lgamma(a + 1) - lgamma(a) - log(b) + m(a + 1) - log_m)
      m2 <- exp(lgamma(a + 2) - lgamma(a) - 2 * log(b) + m(a + 2) - log_m)
      log_mass[i] <- log_m; mean[i] <- m1; variance[i] <- m2 - m1^2
    } else {
      if (a <= 2) stop("the inverse gamma's variance needs shape above 2")
      m <- function(s) log_interval(function(x, low) pgamma(x, s, b, lower.tail = low, log.p = TRUE), 1 / u, 1 / l)
      log_m <- m(a)
      m1 <- exp(lgamma(a - 1) - lgamma(a) + log(b) + m(a - 1) - log_m)
      m2 <- exp(lgamma(a - 2) - lgamma(a) + 2 * log(b) + m(a - 2) - log_m)
      log_mass[i] <- log_m; mean[i] <- m1; variance[i] <- m2 - m1^2
    }
  }
  list(log_mass = log_mass, mean = mean, variance = variance)
}

emit_json(if (args[3] == "quantile") quantiles() else moments(), args[2], na = "string")
