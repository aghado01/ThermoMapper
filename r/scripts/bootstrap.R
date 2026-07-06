# bootstrap.R — one-time provisioning of the R oracle's package library.
#
# Run with R turned on (. env-Rlang.ps1), from the r/ project dir:
#   Rscript scripts/bootstrap.R
# Re-runnable; network required. R is opt-in, so this is a deliberate manual step.

repos <- "https://cloud.r-project.org"

if (!requireNamespace("renv", quietly = TRUE))
  install.packages("renv", repos = repos)

# Initialise renv for THIS project (creates renv/, .Rprofile) on first run.
if (!file.exists("renv.lock"))
  renv::init(bare = TRUE, restart = FALSE)

# Oracle dependencies. jsonlite = fixture/result I/O; the rest are reference
# implementations the C# estimators are validated against.
pkgs <- c(
  "jsonlite",
  "Rdimtools",
  "maotai",
  "T4cluster",
  "TDAstats",
  "TDAkit",
  "Riemann",
  "T4transport",
  "SHT"
)
renv::install(pkgs)

# Pin exact versions from DESCRIPTION + transitive imports into renv.lock.
renv::snapshot(prompt = FALSE)

cat("oracle bootstrap complete:", paste(pkgs, collapse = ", "), "\n")
