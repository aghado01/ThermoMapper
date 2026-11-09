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
#     4,draw,input,2,mean,sd           or one normal component; repeated rows make an equal-weight mixture, truncated to [0, 1];
#     5,pair,drawF,drawG,0,0           a co-active pair (slice 3): Cfg_bass on the two draws, under drawF's laws;
#     6,draw,term,0,sign,knot          a functional factor of a term, the knot in the output's [0, 1] (a functional draw, slice 3);
#     7,0,0,index,t,0                  an output location in [0, 1] at which every functional draw is read.
#   Row 3's last column says what the draw is for: 0 read alone (C_bass), 1 a functional draw, 2 only within a pair.
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
#
# Slice 3: co-active pairs (Rumsey, Hardy, Ahrens and Vander Wiel, arXiv 2311.18146) and functional outputs.
#   pairs: Cfg_bass(mod1, mod2) on two packed draws as written (Cfg_bass_univariate), with slice 0's two corrections patched into a copy
#     of Cfg_bass_univariate (its weight line is C_bass_univariate's, matched once), and the reference: sum of w grad f grad g^T by the
#     same tensor rule, split at both draws' knots. Cfg_bass_univariate shares C_bass_univariate's integrals and so its readings.
#   cfg_default_draws: Cfg_bass without mcmc.use reads mcmc.use = min(length(s2), length(s2')), one pair: the shorter chain's last draw
#     with the same index of the other.
#   cfg_list_mod2: Cfg_bass's list dispatch tests class(mod1) in its third branch where mod2 is meant, so a list in mod2 alone falls
#     through to the model check and stops. (Its recursive calls also pass func.use[i], pca_method, verbose positionally into func.use,
#     func.true and pca_method, read and not run: a scalar model ignores them.)
#   functional: a functional draw read at an output location t the way coactivity reads a functional BASS model: bassfunc2bass_fixed_t
#     folds each term's functional hinges at t, degree one's normalized hinge, into its coefficient and merges the terms that coincide,
#     then C_bass_univariate (corrected) reads the scalar model. The reference folds the same hinges at t and integrates the draw's
#     gradient by the tensor rule. Two of the converter's stops are worked around in the pack and recorded, each run unpadded:
#     functional_one_draw: one MCMC row stops it, since beta[, 1:(1 + maxBasis)] drops a one-row matrix to a vector; the pack has two
#       identical rows.
#     functional_no_drop: it trims beta to max(nbasis) + 1 columns, then tests column max(nbasis after merging) + 2, which exists only
#       if the largest model lost a term at t (merged, or left with no design factor); otherwise it stops, subscript out of bounds. The
#       pack appends a term with one functional hinge and coefficient 0, which every t drops, so C is unchanged.
#     (Read, not run: with one basis function, knotInds[i, 1:1, ] drops to a vector and the duplicate search reads it as a column.)
#   functional_readings, counted in the pinned bodies: bassPCA2bass_fixed_t reads no func.true, its documented argument, and places the
#     basis rows at seq(0, 1, length.out = nt), equispaced on [0, 1] whatever the output locations are; C_bass_pca_v2 defaults func.use to
#     seq(0, 1, length = nt + 2) without its ends, nt points strictly between the grid's, so its default C is read between the training
#     locations; bassfunc2bass_fixed_t finds the functional variable as pdes + pfunc, the last one only; C_bass_pca_v1 and
#     Cfg_bass_pca_v1 return TRUE (run).

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
  ids <- sort(unique(D[D[, 1] %in% c(1, 2, 3, 4, 6), 2]))
  lapply(ids, function(d) {
    R <- D[D[, 2] == d & D[, 1] %in% c(1, 2, 3, 4, 6), , drop = FALSE]
    head <- R[R[, 1] == 3, , drop = FALSE]
    p <- as.integer(head[1, 4]); q <- as.integer(head[1, 5]); functional <- head[1, 6] == 1; alone <- head[1, 6] == 0
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
    ff <- R[R[, 1] == 6, , drop = FALSE]
    func <- lapply(seq_len(M), function(m) {
      f <- ff[ff[, 3] == m, , drop = FALSE]
      data.frame(sign = f[, 5], knot = f[, 6])
    })
    list(draw = d, p = p, q = q, beta = beta, terms = terms, laws = laws, functional = functional, alone = alone, func = func)
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
      others <- Reduce(`*`, lapply(setdiff(seq_len(nrow(tm)), k), function(j) value[, j]), 1)   # the other factors' product, row by row
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

# ---- coactivity with slice 0's corrections -------------------------------------------------------------------------------------
# Each patch replaces one expression of the pinned body, matched by identical() the number of times it states (patch_body in
# _common.R), in copies of XI_FUNC_TNORM_SVW and C_bass_univariate under an environment whose parent is the namespace; XI_FUNC is
# copied unpatched so that it calls the patched second moment.
# XI_FUNC_TNORM_SVW: the truncated normal's second moment, the integral of z^2 phi(z) over [AA, BB], is Delta1 + AA phi(AA) - BB phi(BB);
# the code and the paper's xi_TNorm write - AA phi(AA).
second_moment_patches <- list(
  patch(quote(Z2 <- Delta1 - AA * dnorm(AA) - BB * dnorm(BB)), quote(Z2 <- Delta1 + AA * dnorm(AA) - BB * dnorm(BB))))
# C_bass_univariate: a normal mixture's component takes its share of the truncated mixture, w z/cc, so its truncated moment, which
# divides by z, is weighted once; the code writes w/cc.
weight_patches <- list(
  patch(quote(prior[[i]]$weights <- prior[[i]]$weights/cc), quote(prior[[i]]$weights <- prior[[i]]$weights * prior[[i]]$z/cc)))

coactivity_corrected <- function(second_moment = TRUE, weights = TRUE) {
  env <- new.env(parent = cns)
  assign("XI_FUNC_TNORM_SVW", patch_body(cns$XI_FUNC_TNORM_SVW, if (second_moment) second_moment_patches else list(), env), env)
  assign("XI_FUNC", patch_body(cns$XI_FUNC, list(), env), env)
  univariate <- patch_body(cns$C_bass_univariate, if (weights) weight_patches else list(), env)
  assign("C_bass_univariate", univariate, env)
  function(mod, prior) univariate(mod, prior = prior, mcmc.use = NULL, use_native_scale = FALSE, verbose = FALSE)
}
second_fixed <- coactivity_corrected(second_moment = TRUE, weights = FALSE)
both_fixed <- coactivity_corrected(second_moment = TRUE, weights = TRUE)
tnorm_fixed <- patch_body(cns$XI_FUNC_TNORM_SVW, second_moment_patches, cns)

draws <- read_draws(D)
by_id <- setNames(draws, vapply(draws, function(d) as.character(d$draw), ""))
results <- lapply(Filter(function(d) d$alone, draws), function(d) {
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
       corrected = tnorm_fixed(2, lo, hi, pr, 1),
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

# ---- pairs: Cfg_bass ------------------------------------------------------------------------------------------------------------
# A copy of Cfg_bass_univariate with slice 0's corrections, under the same patched XI_FUNC_TNORM_SVW and XI_FUNC.
cfg_corrected <- function(second_moment = TRUE, weights = TRUE) {
  env <- new.env(parent = cns)
  assign("XI_FUNC_TNORM_SVW", patch_body(cns$XI_FUNC_TNORM_SVW, if (second_moment) second_moment_patches else list(), env), env)
  assign("XI_FUNC", patch_body(cns$XI_FUNC, list(), env), env)
  univariate <- patch_body(cns$Cfg_bass_univariate, if (weights) weight_patches else list(), env)
  function(mod1, mod2, prior) univariate(mod1, mod2, prior = prior, mcmc.use = 1, use_native_scale = FALSE, verbose = FALSE)
}
cfg_fixed <- cfg_corrected()

# The two draws packed with one shared xx.des, which Cfg_bass_univariate compares.
pack_pair <- function(df, dg) {
  pf <- pack(df); pg <- pack(dg)
  pf$xx.des <- pg$xx.des <- matrix(0.5, 1, df$p)
  list(f = pf, g = pg)
}

tensor_Cfg <- function(df, dg, rules) {
  grid <- as.matrix(expand.grid(lapply(rules, function(r) seq_along(r$x))))
  X <- matrix(sapply(seq_len(df$p), function(v) rules[[v]]$x[grid[, v]]), ncol = df$p)
  W <- Reduce(`*`, lapply(seq_len(df$p), function(v) rules[[v]]$w[grid[, v]]))
  crossprod(gradient(df, X) * W, gradient(dg, X))
}

P <- D[D[, 1] == 5, , drop = FALSE]
pairs <- lapply(seq_len(nrow(P)), function(k) {
  df <- by_id[[as.character(P[k, 3])]]; dg <- by_id[[as.character(P[k, 4])]]
  pk <- pack_pair(df, dg)
  Cfg <- suppressWarnings(cns$Cfg_bass(pk$f, pk$g, prior = prior_of(df), mcmc.use = 1))
  Cc <- cfg_fixed(pk$f, pk$g, prior_of(df))
  ref <- tensor_Cfg(df, dg, lapply(seq_len(df$p), function(v) stated_rule(df$laws[[v]], c(breaks_of(df, v), breaks_of(dg, v)))))
  gap <- function(a, b) max(abs(a - b)) / max(abs(b))
  list(pair = P[k, 2], f = df$draw, g = dg$draw, Cfg = unname(Cfg), corrected = unname(Cc), reference = unname(ref),
       gap_reference = gap(Cfg, ref), gap_corrected_reference = gap(Cc, ref))
})

# ---- cfg_default_draws and cfg_list_mod2 ---------------------------------------------------------------------------------------
cfg_default_draws <- NULL
cfg_list_mod2 <- NULL
if (length(pairs)) {
  df <- by_id[[as.character(P[1, 3])]]; dg <- by_id[[as.character(P[1, 4])]]
  pk <- pack_pair(df, dg)
  two <- function(m) {
    m$beta <- rbind(m$beta, m$beta * 2); m$nbasis <- rep(m$nbasis, 2); m$model.lookup <- c(1L, 1L); m$s2 <- c(1, 1); m
  }
  tf <- two(pk$f); tg <- two(pk$g)
  cfg_default_draws <- list(default = unname(suppressWarnings(cns$Cfg_bass(tf, tg, prior = prior_of(df)))),
                            first = unname(suppressWarnings(cns$Cfg_bass(tf, tg, prior = prior_of(df), mcmc.use = 1))),
                            last = unname(suppressWarnings(cns$Cfg_bass(tf, tg, prior = prior_of(df), mcmc.use = 2))))
  cfg_list_mod2 <- tryCatch(list(ran = TRUE, value = length(cns$Cfg_bass(pk$f, list(pk$g), prior = prior_of(df), mcmc.use = 1))),
                            error = function(e) list(ran = FALSE, error = conditionMessage(e)))
}

# ---- functional: bassfunc2bass_fixed_t then C_bass_univariate ------------------------------------------------------------------
locations <- D[D[, 1] == 7, 5]
# A functional draw as a fitted functional bass object: design knots by index into xx.des (one row a factor), functional knots by index
# into xx.func (the locations and the knots, sorted), rows identical MCMC draws of one model.
pack_func <- function(d, rows = 2L, pad = TRUE) {
  if (pad) {
    d$terms <- c(d$terms, list(data.frame(input = integer(0), sign = numeric(0), knot = numeric(0))))
    d$func <- c(d$func, list(data.frame(sign = 1, knot = locations[1])))
    d$beta <- c(d$beta, 0)
  }
  M <- length(d$terms)
  nd <- max(1L, max(c(0L, vapply(d$terms, nrow, 1L))))
  nf <- max(1L, max(c(0L, vapply(d$func, nrow, 1L))))
  grid <- sort(unique(c(locations, unlist(lapply(d$func, function(f) f$knot)))))
  factors <- sum(vapply(d$terms, nrow, 1L))
  xx <- matrix(0.5, max(1L, factors), d$p)
  signs.des <- vars.des <- knotInd.des <- array(NA_real_, c(1, M, nd))
  signs.func <- vars.func <- knotInd.func <- array(NA_real_, c(1, M, nf))
  n.int.des <- n.int.func <- matrix(0, 1, M)
  at <- 0L
  for (m in seq_len(M)) {
    tm <- d$terms[[m]]
    for (k in seq_len(nrow(tm))) {
      at <- at + 1L
      xx[at, tm$input[k]] <- tm$knot[k]
      vars.des[1, m, k] <- tm$input[k]; signs.des[1, m, k] <- tm$sign[k]; knotInd.des[1, m, k] <- at
    }
    n.int.des[1, m] <- nrow(tm)
    fm <- d$func[[m]]
    for (k in seq_len(nrow(fm))) {
      vars.func[1, m, k] <- 1; signs.func[1, m, k] <- fm$sign[k]; knotInd.func[1, m, k] <- match(fm$knot[k], grid)
    }
    n.int.func[1, m] <- nrow(fm)
  }
  structure(list(func = TRUE, pfunc = 1, wasfunc = FALSE, cat = FALSE, pdes = d$p, degree = d$q, n.models = 1L,
                 nbasis = rep(M, rows), model.lookup = rep(1L, rows), s2 = rep(1, rows),
                 beta = matrix(rep(d$beta, rows), nrow = rows, byrow = TRUE),
                 signs.des = signs.des, vars.des = vars.des, knotInd.des = knotInd.des, n.int.des = n.int.des, xx.des = xx,
                 signs.func = signs.func, vars.func = vars.func, knotInd.func = knotInd.func, n.int.func = n.int.func,
                 xx.func = matrix(grid, ncol = 1), y = matrix(0, 1, length(grid)),
                 range.des = matrix(rep(c(0, 1), d$p), nrow = 2), maxInt.des = nd, maxInt.func = nf, cx = rep("numeric", d$p)),
            class = "bass")
}

# The draw at t: each term's coefficient times its functional hinges at t (BASS's normalized hinge, degree q), as a scalar draw.
at_location <- function(d, t) {
  for (m in seq_along(d$terms)) {
    fm <- d$func[[m]]
    if (!nrow(fm)) next
    h <- pmax(fm$sign * (t - fm$knot), 0)^d$q / ((fm$sign + 1) / 2 - fm$sign * fm$knot)^d$q
    d$beta[m + 1] <- d$beta[m + 1] * prod(h)
  }
  d
}

functional <- lapply(Filter(function(d) d$functional, draws), function(d) {
  pk <- pack_func(d)
  per <- lapply(seq_along(locations), function(l) {
    t <- locations[l]
    mod <- suppressWarnings(cns$bassfunc2bass_fixed_t(pk, t))
    Cc <- both_fixed(mod, prior_of(d))
    dt <- at_location(d, t)
    ref <- tensor_C(dt, lapply(seq_len(d$p), function(v) stated_rule(d$laws[[v]], breaks_of(dt, v))))
    list(location = l, t = t, corrected = unname(Cc), reference = unname(ref))
  })
  list(draw = d$draw, locations = per)
})
functional_one_draw <- functional_no_drop <- NULL
fd <- Filter(function(d) d$functional, draws)
if (length(fd) && length(locations)) {
  functional_one_draw <- tryCatch({ cns$bassfunc2bass_fixed_t(pack_func(fd[[1]], rows = 1L), locations[1]); list(ran = TRUE) },
                                  error = function(e) list(ran = FALSE, error = conditionMessage(e)))
  functional_no_drop <- tryCatch({ suppressWarnings(cns$bassfunc2bass_fixed_t(pack_func(fd[[1]], pad = FALSE), locations[1])); list(ran = TRUE) },
                                 error = function(e) list(ran = FALSE, error = conditionMessage(e)))
}

# ---- functional_readings -------------------------------------------------------------------------------------------------------
occurs <- function(fn, e) replace_lang(body(fn), e, e)$count
functional_readings <- list(
  func_true_in_body = "func.true" %in% all.names(body(cns$bassPCA2bass_fixed_t)),
  equispaced_grid = occurs(cns$bassPCA2bass_fixed_t, quote(seq(0, 1, length.out = nt))),
  get_u_grid = occurs(cns$get_u, quote(seq(0, 1, length.out = nrow(U)))),
  pca_v2_default_func_use = occurs(cns$C_bass_pca_v2, quote(func.use <- seq(0, 1, length = nt + 2)[-c(1, nt + 2)])),
  default_points_15 = seq(0, 1, length = 17)[-c(1, 17)],
  last_functional_variable = occurs(cns$bassfunc2bass_fixed_t, quote(func_idx <- which(vars_curr == (bassfunc$pdes + bassfunc$pfunc)))),
  pca_v1 = isTRUE(suppressWarnings(cns$C_bass_pca_v1(NULL, NULL, NULL, NULL, NULL, NULL, NULL))),
  cfg_pca_v1 = isTRUE(cns$Cfg_bass_pca_v1(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL)))

emit_json(list(coactivity = sha, draws = results, normal_second_moment = normal_second_moment, additive = additive,
               gamma_prior = list(built = gamma_built, hand = gamma_hand),
               default_draws = default_draws,
               native_scale_default = formals(cns$C_bass)$use_native_scale,
               pairs = pairs, cfg_default_draws = cfg_default_draws, cfg_list_mod2 = cfg_list_mod2,
               functional = functional, functional_one_draw = functional_one_draw, functional_no_drop = functional_no_drop,
               functional_readings = functional_readings),
          args[2], null = "null", na = "null", matrix = "rowmajor")
