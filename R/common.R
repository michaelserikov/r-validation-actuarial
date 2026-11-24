# Shared helpers: paths, Gini, deviance
`%||%` <- function(a, b) if (is.null(a)) b else a
data_dir <- function(...) file.path(Sys.getenv("R1_ROOT", unset = getwd()), "data", ...)
res_dir <- function(...) { d <- file.path(Sys.getenv("R1_ROOT", unset = getwd()), "results"); dir.create(d, showWarnings = FALSE); file.path(d, ...) }

lorenz_gini <- function(y, pred, w) {
  o <- order(pred)
  cw <- c(0, cumsum(w[o]) / sum(w)); cy <- c(0, cumsum((y * w)[o]) / sum(y * w))
  1 - 2 * sum(diff(cw) * (head(cy, -1) + tail(cy, -1)) / 2)
}
poisson_dev <- function(y, mu, w) {
  t <- ifelse(y > 0, y * log(y / mu), 0) - (y - mu)
  2 * sum(w * t) / sum(w)
}
