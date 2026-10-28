# bass_basis_oracle.R — R BASS's basis evaluation, for the C# port's parity (fitting segment, stage 4's S2).
#
# Usage (cwd = r/):  Rscript oracles/bass_basis_oracle.R <fixture.json> <out.json>
#   fixture: { degree,
#              x:     [[raw values of continuous predictor 0], ...],
#              cat:   [{ levels: [names], values: [one name per case] }, ...],
#              t:     [[raw values of output variable 0], ...]   (empty for a scalar response),
#              terms: [{ cont: [{ var, sign, case }], cat: [{ var, levels: [names] }], func: [{ var, sign, loc }] }] }
#   Indices are 0-based. A continuous factor's knot is the scaled value of training case `case`, a functional
#   factor's that of location `loc`, as R draws knots (knotInd).
#
# Emits { xxt_des, xxt_func, case, output }: the scaled columns (R's scale_range), and per term the case-side factor
#   (makeBasis on the continuous part times makeBasisCat on the categorical part, 1 where a part is empty) and the
#   output-side factor (makeBasis on the functional part). Factors are passed in variable order, as R's births sort
#   them. Everything is R BASS's own code at the pinned commit (basis_funcs.R, misc.R), sourced, not transcribed.

source("oracles/_common.R")
args <- commandArgs(trailingOnly = TRUE)
fx <- jsonlite::read_json(args[1], simplifyVector = FALSE)
env <- bass_source(c("basis_funcs.R", "misc.R"))
q <- fx$degree

scaled_rows <- function(columns) do.call(rbind, lapply(columns, function(v) env$scale_range(unlist(v))))
xxt.des <- if (length(fx$x) > 0) scaled_rows(fx$x) else NULL
xxt.func <- if (length(fx$t) > 0) scaled_rows(fx$t) else NULL
xx.cat <- if (length(fx$cat) > 0)
  as.data.frame(lapply(fx$cat, function(f) factor(unlist(f$values), levels = unlist(f$levels)))) else NULL
n <- if (!is.null(xxt.des)) ncol(xxt.des) else nrow(xx.cat)
m <- if (!is.null(xxt.func)) ncol(xxt.func) else 0

hinge_part <- function(factors, xxt, at_name) {
  vars <- vapply(factors, function(f) f$var, numeric(1)) + 1
  ord <- order(vars)
  vars <- vars[ord]
  signs <- vapply(factors, function(f) f$sign, numeric(1))[ord]
  at <- vapply(factors, function(f) f[[at_name]], numeric(1))[ord] + 1
  knots <- xxt[cbind(vars, at)]
  env$makeBasis(signs, vars, knots, xxt, q)
}

case <- list()
output <- list()
for (term in fx$terms) {
  a <- rep(1, n)
  if (length(term$cont) > 0) a <- hinge_part(term$cont, xxt.des, "case")
  if (length(term$cat) > 0) {
    vars <- vapply(term$cat, function(f) f$var, numeric(1)) + 1
    sub <- lapply(term$cat, function(f) unlist(f$levels))
    ord <- order(vars)
    a <- a * env$makeBasisCat(vars[ord], sub[ord], xx.cat)
  }
  case[[length(case) + 1]] <- a
  if (m > 0) {
    h <- rep(1, m)
    if (length(term$func) > 0) h <- hinge_part(term$func, xxt.func, "loc")
    output[[length(output) + 1]] <- h
  }
}

rows <- function(mat) if (is.null(mat)) list() else lapply(seq_len(nrow(mat)), function(j) mat[j, ])
emit_json(list(xxt_des = rows(xxt.des), xxt_func = rows(xxt.func), case = case, output = output), args[2])
