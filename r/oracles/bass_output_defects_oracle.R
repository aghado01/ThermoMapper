# bass_output_defects_oracle.R — defect fixtures for R BASS's output path (port brief, outputs segment, D07; O-O1).
#
# Usage (cwd = r/):  Rscript oracles/bass_output_defects_oracle.R <Y.csv> <out.json>
#   Y: headerless CSV, one case per row and one output location per column (n x m, m > n > 3), not constant in any
#     column, with residuals asymmetric enough that their sign shows.
#
# Each section runs the pinned R code (bassPCA.R, predict.R, misc.R, plot.R) as written, on Y or on small inputs built
# here, and records R's result beside the corrected expectation from D07:
#   truncation_sampler: truncErrSampN indexes the residual matrix's columns (cases) with its row count (locations).
#     wide: Y's own residuals (m > n), where R's draw fails; narrow: Y[, 1:3] (3 < n), where R draws only the first
#     three cases. Corrected: whole case vectors, uniformly over all n cases.
#   noise_recycling: predict1mod adds one rnorm per posterior draw, recycled across the prediction cases. On a
#     zero-mean mock component (200 draws, 5 cases, s2 = 1): R's largest spread across cases within a draw and the
#     correlation across cases; corrected noise is independent per case and draw.
#   residual: bassPCAsetup(scale = TRUE, perc.var = 90) stores reconstructed minus observed, in standardized units,
#     and predict.bassBasis adds it after restoring physical units. R's training reconstruction (predict1mcmc's
#     physical mean plus R's residual) against D07's y = y_m + y_s * (K w + e), e = observed minus reconstructed;
#     the residuals' pooled skewness under both orientations.
#   full_retention: bassPCAsetup(perc.var = 100), whose strict '>' comparison finds no rank; R's n.pc expression
#     evaluated as written; the corrected rank is the numerical rank of the centred data at max(n, m) eps d_1.
#   constant_columns: a constant input through scale_range (training) and scale_range_mat (prediction), at the
#     training value and at a new one (R's prediction values as strings, since NaN and Inf are not JSON); a constant
#     output column through bassPCAsetup(scale = TRUE). Corrected (D07, decided 2026-10-29): the input is refused at
#     entry, since it has no knot, as a one-level factor is; the output column gets y_s = 1.
#   draw_order: predict.bass on a mock fit whose mult function writes each draw's index, mcmc.use = (4, 1, 5, 2) over
#     model groups (1, 1, 2, 2, 3): the rows R returns, the dimensions it keeps for 4 x 3, 4 x 1 and 1 x 1 requests,
#     and with nugget = TRUE (s2 = 0 except draw 2) which rows carry noise. Corrected: requested order, explicit
#     draw x case x output axes, each row's noise from its own draw.
#   interval_quantiles: plot.bass's interval (the quantile over draws of each draw's Gaussian interval endpoint),
#     evaluated as written, against the posterior predictive mixture's quantile found by root finding on its CDF.
#
# Emits { n, m, truncation_sampler, noise_recycling, residual, full_retention, constant_columns, draw_order,
#         interval_quantiles }; null where R fails.

source("oracles/_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: bass_output_defects_oracle.R <Y.csv> <out.json>")
Y <- read_matrix(args[1])
n <- nrow(Y); m <- ncol(Y)
if (m <= n || n <= 3) stop("need m > n > 3")
bass <- bass_source(c("misc.R", "basis_funcs.R", "predict.R", "plot.R", "bassPCA.R"))
xx <- matrix(seq_len(n) / n)
error_of <- function(expr) tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))

# ---- truncation sampler ---------------------------------------------------------------------------------------------
wide <- bass$bassPCAsetup(xx, Y, perc.var = 90)
narrow <- bass$bassPCAsetup(xx, Y[, 1:3], perc.var = 90)
set.seed(11)
draws_narrow <- bass$truncErrSampN(2000, narrow$trunc.error)
case_of <- function(rows, te) apply(rows, 1, function(r) which(colSums(abs(te - r)) == 0)[1])
set.seed(11)
corrected_idx <- sample.int(n, 2000, replace = TRUE)
truncation_sampler <- list(
  wide = list(locations = nrow(wide$trunc.error), cases = ncol(wide$trunc.error),
              r_error = error_of(bass$truncErrSampN(100, wide$trunc.error))),
  narrow = list(locations = nrow(narrow$trunc.error), cases = ncol(narrow$trunc.error),
                r_cases_drawn = sort(unique(case_of(draws_narrow, narrow$trunc.error))),
                corrected_cases_drawn = sort(unique(corrected_idx))))

# ---- noise recycling ------------------------------------------------------------------------------------------------
predict.oracle_zero <- function(object, newdata, mcmc.use, ...) matrix(0, length(mcmc.use), nrow(newdata))
zero_mod <- structure(list(s2 = rep(1, 200)), class = "oracle_zero")
set.seed(12)
pm <- bass$predict1mod(zero_mod, matrix(0, 5, 1), 1:200, TRUE)
set.seed(12)
corrected_noise <- matrix(rnorm(200 * 5), 200, 5)
noise_recycling <- list(draws = 200, cases = 5,
                        r_max_spread_across_cases = max(apply(pm, 1, function(r) max(r) - min(r))),
                        r_correlation_cases_1_2 = cor(pm[, 1], pm[, 2]),
                        corrected_max_spread_across_cases = max(apply(corrected_noise, 1, function(r) max(r) - min(r))),
                        corrected_correlation_cases_1_2 = cor(corrected_noise[, 1], corrected_noise[, 2]))

# ---- residual sign and units ----------------------------------------------------------------------------------------
dat <- bass$bassPCAsetup(xx, Y, perc.var = 90, scale = TRUE)
yc <- t(sweep(sweep(Y, 2, dat$y.m), 2, dat$y.s, "/"))                 # m x n, standardized
e <- yc - dat$basis %*% dat$newy                                      # observed minus reconstructed
physical_mean <- vapply(seq_len(n), function(i) c(bass$predict1mcmc(matrix(dat$newy[, i], 1), dat)), numeric(m))
r_recon <- physical_mean + dat$trunc.error                            # R's predictive path, residual as stored
corrected_recon <- physical_mean + dat$y.s * e
skew <- function(v) mean((v - mean(v))^3) / mean((v - mean(v))^2)^1.5
residual <- list(r = dat$n.pc, y_s_range = range(dat$y.s),
                 r_stored_is_negative_of_corrected = max(abs(dat$trunc.error + e)),
                 r_reconstruction_max_error = max(abs(r_recon - t(Y))),
                 corrected_reconstruction_max_error = max(abs(corrected_recon - t(Y))),
                 scale = max(abs(t(Y))),
                 skewness_corrected = skew(c(e)), skewness_r = skew(c(dat$trunc.error)))

# ---- full retention -------------------------------------------------------------------------------------------------
Yc <- sweep(Y, 2, colMeans(Y))
d <- svd(Yc)$d
ev <- svd(t(Yc))$d^2
n_pc_expr <- source_assignment(bass$bassPCAsetup, "n.pc")
r_n_pc <- eval(n_pc_expr, list2env(list(ev = ev, perc.var = 100)))
full_retention <- list(r_error = error_of(bass$bassPCAsetup(xx, Y, perc.var = 100)), r_n_pc = r_n_pc,
                       cumulative_share_last_minus_100 = tail(cumsum(ev / sum(ev)) * 100, 1) - 100,
                       corrected_rank = sum(d > max(n, m) * .Machine$double.eps * d[1]))

# ---- constant columns -----------------------------------------------------------------------------------------------
x_const <- rep(0.3, n)
r_train <- bass$scale_range(x_const)
range_mat <- matrix(range(x_const), 2, 1)
r_predict <- c(bass$scale_range_mat(matrix(c(0.3, 0.5), 2, 1), range_mat))
Y_const <- Y; Y_const[, 2] <- 7
constant_columns <- list(
  input = list(training_value = 0.3, new_value = 0.5, r_training = r_train[1], r_prediction = as.character(r_predict),
               corrected = "refused"),
  output = list(r_error = error_of(bass$bassPCAsetup(xx, Y_const, perc.var = 90, scale = TRUE)),
                r_y_s = apply(Y_const, 2, sd)[2], corrected_y_s = 1))

# ---- draw order and axes --------------------------------------------------------------------------------------------
assign("mult_oracle_mock", function(model, mcmc.use.mod, object, tnewdata.des, newdata.cat, tnewdata.func)
  array(rep(mcmc.use.mod, ncol(tnewdata.des)), dim = c(length(mcmc.use.mod), ncol(tnewdata.des), 1)), envir = bass)
mock <- list(nmcmc = 5, nburn = 0, thin = 1, func = FALSE, pdes = 1, pcat = 0, cx = "numeric",
             range.des = matrix(c(0, 1), 2, 1), model.lookup = c(1, 1, 2, 2, 3), type = "_oracle_mock",
             s2 = c(0, 1, 0, 0, 0))
use <- c(4, 1, 5, 2)
out43 <- bass$predict.bass(mock, matrix(c(0.2, 0.5, 0.7)), mcmc.use = use)
out41 <- bass$predict.bass(mock, matrix(0.2), mcmc.use = use)
out11 <- bass$predict.bass(mock, matrix(0.2), mcmc.use = 4)
set.seed(13)
noisy <- bass$predict.bass(mock, matrix(c(0.2, 0.5, 0.7)), mcmc.use = use, nugget = TRUE)
draw_order <- list(requested = use, model_lookup = mock$model.lookup, r_rows = out43[, 1],
                   r_dim_4x3 = dim(out43), r_dim_4x1 = if (is.null(dim(out41))) length(out41) else dim(out41),
                   r_dim_1x1 = if (is.null(dim(out11))) length(out11) else dim(out11),
                   r_rows_with_noise = out43[, 1][apply(noisy != out43, 1, any)],
                   corrected_rows = use, corrected_dims = list(c(4, 3, 1), c(4, 1, 1), c(1, 1, 1)),
                   corrected_rows_with_noise = use[mock$s2[use] > 0])

# ---- interval quantiles ---------------------------------------------------------------------------------------------
# Case 1: means spread over the draws, unit sd; case 2: equal means, sds spread over the draws. plot.bass takes one sd
# per draw, so each case gets its own evaluation.
D <- 1000
u <- qnorm((seq_len(D) - 0.5) / D)
quants <- c(0.025, 0.975)
r_interval <- function(mu, sd) {
  e <- list2env(list(x = list(yhat = matrix(mu), s2 = sd^2), s = sd, quants = quants, margin = 2))
  c(eval(source_assignment(bass$plot.bass, "qq2"), e), eval(source_assignment(bass$plot.bass, "qq1"), e))
}
mixture_quantile <- function(mu, sd, prob)
  uniroot(function(q) mean(pnorm((q - mu) / sd)) - prob, range(mu) + c(-10, 10) * max(sd), tol = 1e-13)$root
interval_quantiles <- lapply(list(list(mu = u, sd = rep(1, D)), list(mu = rep(0, D), sd = exp(u / 2))), function(cs)
  list(r = r_interval(cs$mu, cs$sd), mixture = vapply(quants, function(q) mixture_quantile(cs$mu, cs$sd, q), 0)))

emit_json(list(n = n, m = m, truncation_sampler = truncation_sampler, noise_recycling = noise_recycling,
               residual = residual, full_retention = full_retention, constant_columns = constant_columns,
               draw_order = draw_order, interval_quantiles = interval_quantiles), args[2], na = "null")
