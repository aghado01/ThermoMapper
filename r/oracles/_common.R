# Shared oracle helpers.
# Fixtures are headerless CSV (a plain numeric matrix); results are JSON, read
# C#-side via the jso-jackson tools. Run oracles with the r/ dir as cwd.

read_matrix <- function(path) {
  as.matrix(read.csv(path, header = FALSE))
}

emit_json <- function(obj, path) {
  jsonlite::write_json(obj, path, digits = 16, auto_unbox = TRUE, pretty = TRUE)
}
