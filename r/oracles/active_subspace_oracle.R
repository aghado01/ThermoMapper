# active_subspace_oracle.R — coactivity's C_bass (knrumsey/coactivity@35e4573, the closed form of Rumsey et al., arXiv 2307.11241) on
# given hinge draws, beside a reference that shares nothing with its code, and its defects, each confirmed or struck (E8's build brief,
# issues/regression/bass/active-subspace-brief.md, slice 0).
#
# Usage (cwd = r/):  Rscript oracles/active_subspace_oracle.R <draws.csv> <out.json>
#   draws: headerless numeric CSV, six columns, draws and terms counted from 1, inputs from 1, term 0 the intercept:
#     1,draw,term,coefficient,0,0      one per coefficient, the intercept's included;
#     2,draw,term,input,sign,knot      one per hinge factor, the knot in [0, 1] (BASS's scaled units, the box [0, 1] per input);
#     3,draw,0,p,degree,0              one per draw: its inputs and its hinges' degree q;
#     4,draw,input,1,lower,upper       the input's law: uniform on [lower, upper] within [0, 1];
#     4,draw,input,2,mean,sd           or one normal component; repeated rows make an equal-weight mixture, truncated to [0, 1].
#   A draw is f(x) = b0 + sum_m b_m prod_v [s_v (x_v - t_v)]_+^q / ((s_v + 1)/2 - s_v t_v)^q, BASS's normalized hinges (R BASS's
#   makeBasis and const). Each draw is packed alone as a fitted bass object, its knots given by value (wasfunc), with its degree
#   and maxInt at least two, and C_bass reads it with the prior its law gives.
#
# Emits, per draw: C (C_bass as written, p x p by row); second_moment_fixed and corrected (C_bass with the patches below: the
#   second moment alone, then the weights too); reference (the stated law: the mixture renormalized on [0, 1] as one law, R's
#   convention); as_read (the law C_bass's prior handling integrates against, set out below); and the gaps. The reference and
#   as_read integrate the draw's own gradient by tensor Gauss-Legendre on panels split at the input's knots and the law's ends,
#   24 nodes a panel: exact for a uniform law, to rounding for a normal one.
#
# Readings of the code, each checked by running it:
#   degree_read: C_bass reads no degree. Its integrals I1, I2, I3 are degree one's, and its hinge scale 1/((s + 1)/2 - s t) is
#     degree one's const. Shown on the draws: the parity test sends a degree-two draw and the same knots at degree one.
#   normal_second_moment: XI_FUNC_TNORM_SVW's Z2 = Delta1 - AA phi(AA) - BB phi(BB), the paper's xi_TNorm, flips the lower end's term:
#     the integral of z^2 phi(z) over [AA, BB] is Delta1 + AA phi(AA) - BB phi(BB). Emitted at four intervals beside integrate and
#     the patched code; on the draws, entries whose term pairs share a hinged normal input carry it.
#   normal_weights: C_bass sets each component's weight to w/cc, cc = sum_l w_l z_l (z_l the component's mass on the truncation),
#     then multiplies by the component's own truncated moment, which divides by z_l again. The law it integrates is
#     sum_l (w_l/cc) phi_l(x)/z_l on the truncation, of mass sum_l w_l/(cc z_l) and not 1; the comment beside it says the intent is
#     to keep the mixture's shape. as_read is that law, and the draws show C_bass meeting it once the second moment is patched.
#   additive: maxInt = 1 packs one factor slot a term; the slices drop a dimension and C_bass stops (the paper's supplement says
#     concordance cannot do maxInt = 1). Recorded as the error it raises.
#   gamma_prior: build_prior's gamma takes shape and scale; C_bass's gamma reads prior$scale as the shape and prior$rate as the rate,
#     multiplies the rate by the scaled truncation's upper end (1 by default, so the rate stays in native units), and does not
#     renormalize on the truncation. Recorded: build_prior's prior through C_bass, and a hand-set prior's E0 on [0, 1] against the
#     law's mass there.
#   default_draws: C_bass without mcmc.use reads the last draw only.
#   native_scale_default: use_native_scale's default in the signature (the roxygen says TRUE).

source("oracles/_common.R")
suppressPackageStartupMessages(library(coactivity))
cns <- asNamespace("coactivity")
sha <- utils::packageDescription("coactivity")$RemoteSha
if (is.null(sha) || !startsWith(sha, "35e4573")) stop("coactivity is not at the pinned 35e4573")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: active_subspace_oracle.R <draws.csv> <out.json>")
D <- read_matrix(args[1])

# ---- the draws -----------------------------------------------------------------------------------------------------------------
read_draws <- function(D) {
  ids <- sort(unique(D[, 2]))
  lapply(ids, function(d) {
    R <- D[D[, 2] == d, , drop = FALSE]
    head <- R[R[, 1] == 3, , drop = FALSE]
    p <- as.integer(head[1, 4]); q <- as.integer(head[1, 5])
    coef <- R[R[, 1] == 1, , drop = FALSE]
    M <- if (nrow(coef)) as.integer(max(coef[, 3])) else 0L
    beta <- numeric(M + 1); beta[coef[, 3] + 1] <- coef[, 4]
    fac <- R[R[, 1] == 2, , drop = FALSE]
    terms <- lapply(seq_len(M), function(m) {
      f <- fac[fac[, 3] == m, , drop = FALSE]
      data.frame(input = as.integer(f[, 4]), sign = f[, 5], knot = f[, 6])
    })
    law <- R[R[, 1] == 4, , drop = FALSE]
    laws <- lapply(seq_len(p), function(v) {
      l <- law[law[, 3] == v, , drop = FALSE]
      if (!nrow(l) || l[1, 4] == 1) list(kind = "uniform", lower = if (nrow(l)) l[1, 5] else 0, upper = if (nrow(l)) l[1, 6] else 1)
      else list(kind = "normal", mean = l[, 5], sd = l[, 6])
    })
    list(draw = d, p = p, q = q, beta = beta, terms = terms, laws = laws)
  })
}

# The draw as a fitted bass object with one model and one MCMC draw: knots by value (wasfunc), the box [0, 1] per input.
pack <- function(d, maxInt = max(2L, max(c(1L, vapply(d$terms, nrow, 1L))))) {
  M <- length(d$terms)
  signs <- vars <- knots <- array(NA_real_, c(1, max(M, 1L), maxInt))
  for (m in seq_len(M)) {
    k <- nrow(d$terms[[m]])
    vars[1, m, seq_len(k)] <- d$terms[[m]]$input
    signs[1, m, seq_len(k)] <- d$terms[[m]]$sign
    knots[1, m, seq_len(k)] <- d$terms[[m]]$knot
  }
  structure(list(func = FALSE, pfunc = 0, wasfunc = TRUE, pdes = d$p, degree = d$q, nbasis = M, model.lookup = 1L,
                 beta = matrix(d$beta, nrow = 1), signs.des = signs, vars.des = vars, knots.des = knots,
                 range.des = matrix(rep(c(0, 1), d$p), nrow = 2)),
            class = "bass")
}

prior_of <- function(d) lapply(d$laws, function(l)
  if (l$kind == "uniform") list(dist = "uniform", trunc = c(l$lower, l$upper))
  else list(dist = "normal", mean = l$mean, sd = l$sd))

# ---- the reference -------------------------------------------------------------------------------------------------------------
# One input's rule on [0, 1] for a law given as a density on [lo, hi]: Gauss-Legendre panels split at the breakpoints.
rule <- function(density, lo, hi, breaks, n = 24) {
  cuts <- sort(unique(c(lo, hi, breaks[breaks > lo & breaks < hi])))
  x <- w <- numeric(0)
  for (k in seq_len(length(cuts) - 1)) {
    g <- gauss_legendre(n, cuts[k], cuts[k + 1])
    x <- c(x, g$x); w <- c(w, g$w * density(g$x))
  }
  list(x = x, w = w)
}

# The stated law: uniform on [lower, upper]; the equal-weight normal mixture renormalized on [0, 1].
stated_rule <- function(l, breaks) {
  if (l$kind == "uniform") return(rule(function(x) rep(1 / (l$upper - l$lower), length(x)), l$lower, l$upper, breaks))
  w <- rep(1 / length(l$mean), length(l$mean))
  z <- pnorm((1 - l$mean) / l$sd) - pnorm((0 - l$mean) / l$sd)
  rule(function(x) Reduce(`+`, lapply(seq_along(w), function(k) w[k] * dnorm(x, l$mean[k], l$sd[k]))) / sum(w * z),
       0, 1, breaks)
}

# The law C_bass integrates against: the stated one for a uniform law; for a normal mixture sum_l (w_l/cc) phi_l/z_l on [0, 1].
as_read_rule <- function(l, breaks) {
  if (l$kind == "uniform") return(stated_rule(l, breaks))
  w <- rep(1 / length(l$mean), length(l$mean))
  z <- pnorm((1 - l$mean) / l$sd) - pnorm((0 - l$mean) / l$sd)
  cc <- sum(w * z)
  rule(function(x) Reduce(`+`, lapply(seq_along(w), function(k) (w[k] / cc) * dnorm(x, l$mean[k], l$sd[k]) / z[k])),
       0, 1, breaks)
}

# The draw's gradient at the rows of X (N x p), from the hinge definition.
gradient <- function(d, X) {
  G <- matrix(0, nrow(X), d$p)
  q <- d$q
  for (m in seq_along(d$terms)) {
    tm <- d$terms[[m]]
    scale <- ((tm$sign + 1) / 2 - tm$sign * tm$knot)^q
    value <- sapply(seq_len(nrow(tm)), function(k) pmax(tm$sign[k] * (X[, tm$input[k]] - tm$knot[k]), 0)^q / scale[k])
    value <- matrix(value, nrow = nrow(X))
    for (k in seq_len(nrow(tm))) {
      u <- tm$sign[k] * (X[, tm$input[k]] - tm$knot[k])
      slope <- q * tm$sign[k] * ifelse(u > 0, u^(q - 1), 0) / scale[k]
      others <- if (nrow(tm) > 1) apply(value[, -k, drop = FALSE], 1, prod) else 1
      G[, tm$input[k]] <- G[, tm$input[k]] + d$beta[m + 1] * slope * others
    }
  }
  G
}

tensor_C <- function(d, rules) {
  grid <- as.matrix(expand.grid(lapply(rules, function(r) seq_along(r$x))))
  X <- sapply(seq_len(d$p), function(v) rules[[v]]$x[grid[, v]])
  X <- matrix(X, ncol = d$p)
  W <- Reduce(`*`, lapply(seq_len(d$p), function(v) rules[[v]]$w[grid[, v]]))
  G <- gradient(d, X)
  crossprod(G * W, G)
}

breaks_of <- function(d, v) {
  k <- unlist(lapply(d$terms, function(tm) tm$knot[tm$input == v]))
  l <- d$laws[[v]]
  c(k, if (l$kind == "uniform") c(l$lower, l$upper) else numeric(0))
}

# ---- coactivity with corrections patched in ------------------------------------------------------------------------------------
# C_bass_univariate, XI_FUNC and XI_FUNC_TNORM_SVW copied into an environment under the namespace, each patch matched exactly once in
# the deparsed body:
#   second_moment: Z2 = Delta1 + AA phi(AA) - BB phi(BB), the truncated normal's second moment (the code and the paper's equation for
#     xi_TNorm have - AA phi(AA));
#   weights: each component's weight w z/cc, its share of the truncated mixture, so the component's truncated moment is weighted once.
patched <- function(fn, from, to) {
  text <- deparse(fn, width.cutoff = 500L)
  hits <- sum(lengths(regmatches(text, gregexpr(from, text, fixed = TRUE))))
  if (hits != 1) stop("expected one match of ", from, ", found ", hits)
  eval(parse(text = sub(from, to, text, fixed = TRUE)))
}
coactivity_corrected <- function(second_moment = TRUE, weights = TRUE) {
  env <- new.env(parent = cns)
  tnorm <- cns$XI_FUNC_TNORM_SVW
  if (second_moment) tnorm <- patched(tnorm, "Z2 <- Delta1 - AA * dnorm(AA) - BB * dnorm(BB)", "Z2 <- Delta1 + AA * dnorm(AA) - BB * dnorm(BB)")
  univariate <- cns$C_bass_univariate
  if (weights) univariate <- patched(univariate, "prior[[i]]$weights <- prior[[i]]$weights/cc", "prior[[i]]$weights <- prior[[i]]$weights * prior[[i]]$z/cc")
  xi <- cns$XI_FUNC
  for (f in c("tnorm", "univariate", "xi")) { g <- get(f); environment(g) <- env; assign(f, g) }
  assign("XI_FUNC_TNORM_SVW", tnorm, env); assign("XI_FUNC", xi, env); assign("C_bass_univariate", univariate, env)
  function(mod, prior) univariate(mod, prior = prior, mcmc.use = NULL, use_native_scale = FALSE, verbose = FALSE)
}
second_fixed <- coactivity_corrected(second_moment = TRUE, weights = FALSE)
both_fixed <- coactivity_corrected(second_moment = TRUE, weights = TRUE)

draws <- read_draws(D)
results <- lapply(draws, function(d) {
  C <- cns$C_bass(pack(d), prior = prior_of(d))
  C2 <- second_fixed(pack(d), prior_of(d))
  Cc <- both_fixed(pack(d), prior_of(d))
  ref <- tensor_C(d, lapply(seq_len(d$p), function(v) stated_rule(d$laws[[v]], breaks_of(d, v))))
  as_read <- tensor_C(d, lapply(seq_len(d$p), function(v) as_read_rule(d$laws[[v]], breaks_of(d, v))))
  gap <- function(a, b) max(abs(a - b)) / max(abs(b))
  list(draw = d$draw, p = d$p, degree = d$q, C = unname(C), second_moment_fixed = unname(C2), corrected = unname(Cc),
       reference = unname(ref), as_read = unname(as_read),
       gap_reference = gap(C, ref), gap_second_fixed_as_read = gap(C2, as_read), gap_corrected_reference = gap(Cc, ref))
})

# ---- normal_second_moment ------------------------------------------------------------------------------------------------------
# XI_FUNC_TNORM_SVW(2, a, b) as written against integrate of x^2 times the component's density on [a, b] over its mass on [0, 1].
normal_second_moment <- lapply(list(c(0.5, 0.5, 0.3, 1), c(0.2, 0.3, 0, 0.6), c(0.7, 0.1, 0.65, 0.9), c(0.5, 0.2, -Inf, 0.4)), function(s) {
  pr <- list(mean = s[1], sd = s[2], trunc = c(0, 1))
  lo <- max(s[3], 0); hi <- min(s[4], 1)
  z <- pnorm((1 - s[1]) / s[2]) - pnorm((0 - s[1]) / s[2])
  list(mean = s[1], sd = s[2], lower = lo, upper = hi, code = cns$XI_FUNC_TNORM_SVW(2, lo, hi, pr, 1),
       corrected = environment(second_fixed)$tnorm(2, lo, hi, pr, 1),
       integral = integrate(function(x) x^2 * dnorm(x, s[1], s[2]), lo, hi, rel.tol = 1e-13)$value / z)
})

# ---- additive: maxInt = 1 ------------------------------------------------------------------------------------------------------
additive_draw <- list(draw = 0, p = 2, q = 1, beta = c(0.3, 1.2, -0.8),
                      terms = list(data.frame(input = 1L, sign = 1, knot = 0.3), data.frame(input = 2L, sign = -1, knot = 0.6)),
                      laws = list(list(kind = "uniform", lower = 0, upper = 1), list(kind = "uniform", lower = 0, upper = 1)))
additive <- tryCatch(list(ran = TRUE, C = unname(cns$C_bass(pack(additive_draw, maxInt = 1L), prior = prior_of(additive_draw)))),
                     error = function(e) list(ran = FALSE, error = conditionMessage(e)))
additive$reference <- unname(tensor_C(additive_draw, lapply(1:2, function(v) stated_rule(additive_draw$laws[[v]], breaks_of(additive_draw, v)))))

# ---- gamma_prior ---------------------------------------------------------------------------------------------------------------
# One input, f = x (a single hinge s = +1, t = 0, so C = E[1] = the law's mass), the box [0, 1].
linear_draw <- list(draw = 0, p = 1, q = 1, beta = c(0, 1), terms = list(data.frame(input = 1L, sign = 1, knot = 0)),
                    laws = list(list(kind = "uniform", lower = 0, upper = 1)))
built <- coactivity::build_prior(dist = "gamma", shape = 2, scale = 0.5)
gamma_built <- tryCatch(list(ran = TRUE, C = unname(cns$C_bass(pack(linear_draw), prior = built))),
                        error = function(e) list(ran = FALSE, error = conditionMessage(e)))
gamma_hand <- list(shape = 2, rate = 3,
                   C = unname(cns$C_bass(pack(linear_draw), prior = list(list(dist = "gamma", scale = 2, rate = 3)))),
                   mass = pgamma(1, shape = 2, rate = 3))

# ---- default_draws -------------------------------------------------------------------------------------------------------------
two <- pack(draws[[1]])
two$beta <- rbind(two$beta, two$beta * 2)
two$nbasis <- rep(two$nbasis, 2)
two$model.lookup <- c(1L, 1L)
default_draws <- list(default = unname(cns$C_bass(two, prior = prior_of(draws[[1]]))),
                      first = unname(cns$C_bass(two, prior = prior_of(draws[[1]]), mcmc.use = 1)),
                      last = unname(cns$C_bass(two, prior = prior_of(draws[[1]]), mcmc.use = 2)))

emit_json(list(coactivity = sha, draws = results, normal_second_moment = normal_second_moment, additive = additive,
               gamma_prior = list(built = gamma_built, hand = gamma_hand),
               default_draws = default_draws,
               native_scale_default = formals(cns$C_bass)$use_native_scale),
          args[2], null = "null", na = "null", matrix = "rowmajor")
