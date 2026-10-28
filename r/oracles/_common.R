# Shared oracle helpers.
# Fixtures are headerless CSV (a plain numeric matrix); results are JSON, read
# C#-side via the jso-jackson tools. Run oracles with the r/ dir as cwd.

read_matrix <- function(path) {
  as.matrix(read.csv(path, header = FALSE))
}

emit_json <- function(obj, path, ...) {
  jsonlite::write_json(obj, path, digits = 16, auto_unbox = TRUE, pretty = TRUE, ...)
}

# The pinned R BASS source (BASS 1.3.1, github.com/lanl/BASS), sourced rather than installed: R BASS is not in the
# project's renv, and its imports (truncdist, hypergeo) sit on paths the oracles do not call. BASS_R_SOURCE locates
# the checkout (default D:/aghado01/packages/BASS), which must be at the pinned commit.
bass_pinned_commit <- "26710146c623a258854162367a939041efedd091"

bass_source <- function(files, root = Sys.getenv("BASS_R_SOURCE", "D:/aghado01/packages/BASS")) {
  head <- suppressWarnings(tryCatch(system2("git", c("-C", root, "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
                                    error = function(e) NA_character_))
  if (!identical(as.vector(head), bass_pinned_commit))
    stop("the R BASS checkout at ", root, " is not at the pinned commit ", bass_pinned_commit)
  env <- new.env(parent = globalenv())
  for (f in files) sys.source(file.path(root, "R", f), envir = env)
  env
}

# Calls in a function's body (searched recursively) for which keep(call) is TRUE, in source order. Lets an oracle
# evaluate a pinned R BASS expression as written instead of a transcription of it.
source_calls <- function(fn, keep) {
  found <- list()
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    if (isTRUE(keep(e))) found[[length(found) + 1]] <<- e
    args <- as.list(e)[-1]
    for (i in seq_along(args)) if (!identical(args[[i]], quote(expr = ))) walk(args[[i]])
  }
  walk(body(fn))
  found
}

# The right-hand side of the n-th assignment `name <- ...` in a function's body.
source_assignment <- function(fn, name, n = 1) {
  found <- source_calls(fn, function(e) as.character(e[[1]])[1] %in% c("<-", "=") && identical(e[[2]], as.name(name)))
  if (length(found) < n) stop("no assignment number ", n, " to ", name)
  found[[n]][[3]]
}

# Gauss-Legendre nodes and weights on [a, b] (Golub-Welsch), exact for polynomials of degree 2n - 1.
gauss_legendre <- function(n, a = -1, b = 1) {
  k <- seq_len(n - 1)
  J <- matrix(0, n, n)
  J[cbind(k, k + 1)] <- J[cbind(k + 1, k)] <- k / sqrt(4 * k^2 - 1)
  e <- eigen(J, symmetric = TRUE)
  o <- order(e$values)
  list(x = (a + b) / 2 + (b - a) / 2 * e$values[o], w = (b - a) / 2 * 2 * e$vectors[1, o]^2)
}
