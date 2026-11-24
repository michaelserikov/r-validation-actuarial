# Module C: tails. GPD (peaks over threshold) on the freMTPL2 training claims: MLE (extRemes, evd), probability weighted moments, profile CI.
suppressPackageStartupMessages({library(extRemes); library(evd)})
source("R/common.R")
set.seed(1)
cl <- read.csv(data_dir("fremtpl2_claims.csv")); pol <- read.csv(data_dir("fremtpl2_policies.csv"))
x <- cl$ClaimAmount[cl$IDpol %in% pol$IDpol[pol$split == "train"]]
cat("training claims:", length(x), "\n")
# Hosking-Wallis probability weighted moments for the GPD (own implementation), unbiased plotting positions
pwm_gpd <- function(exc) {
  y <- sort(exc); n <- length(y); i <- seq_len(n)
  a0 <- mean(y); a1 <- sum((n - i) / (n - 1) * y) / n
  xi <- 2 - a0 / (a0 - 2 * a1); sg <- 2 * a0 * a1 / (a0 - 2 * a1)
  c(shape = xi, scale = sg)
}
# Own maximum likelihood (BFGS on log scale, started at the PWM estimate) and profile likelihood interval for the shape
mle_own <- function(exc) {
  st <- pwm_gpd(exc)
  nll <- function(par) -sum(evd::dgpd(exc, loc = 0, scale = exp(par[1]), shape = par[2], log = TRUE))
  o <- optim(c(log(st[["scale"]]), st[["shape"]]), nll, method = "Nelder-Mead", control = list(reltol = 1e-14, maxit = 5000))
  o <- optim(o$par, nll, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  c(shape = o$par[2], scale = exp(o$par[1]), ll = -o$value)
}
profile_ci <- function(exc, fit) {
  grid <- seq(0.2, 2.0, by = 0.005)
  pl <- sapply(grid, function(xi) max(sapply(c(0.5, 1, 2), function(f) -optimize(function(ls) -sum(evd::dgpd(exc, 0, exp(ls), xi, log = TRUE)), log(fit[["scale"]]) + c(-3, 3))$objective)))
  keep <- pl >= fit[["ll"]] - qchisq(0.95, 1) / 2
  range(grid[keep])
}
rows <- list()
for (q in c(0.90, 0.95, 0.975, 0.99)) {
  u <- unname(quantile(x, q)); exc <- x[x > u] - u
  m <- fevd(x, threshold = u, type = "GP", method = "MLE")
  own <- mle_own(exc); pci <- profile_ci(exc, own)
  fp <- evd::fpot(x, threshold = u, std.err = FALSE)
  pwm <- list(results = pwm_gpd(exc))
  ci <- try(ci(m, type = "parameter", which.par = 2, method = "proflik", xrange = c(0.3, 2), verbose = FALSE), silent = TRUE)
  pr <- if (inherits(ci, "try-error")) c(NA, NA, NA) else as.numeric(ci)
  ll <- function(xi, sg) sum(evd::dgpd(exc, loc = 0, scale = sg, shape = xi, log = TRUE))
  ll_mle <- ll(m$results$par[["shape"]], m$results$par[["scale"]]); ll_fp <- ll(fp$estimate["shape"], fp$estimate["scale"])
  ll_pw <- ll(pwm$results[["shape"]], pwm$results[["scale"]])
  rows[[length(rows) + 1]] <- data.frame(xi_mle_own = own[["shape"]], sigma_mle_own = own[["scale"]], ll_mle_own = own[["ll"]], xi_prof_own_lo = pci[1], xi_prof_own_hi = pci[2], ll_mle = ll_mle, ll_fpot = ll_fp, ll_pwm = ll_pw, quantile = q, threshold = u, n_above = length(exc),
    xi_mle = m$results$par[["shape"]], sigma_mle = m$results$par[["scale"]],
    xi_evd_fpot = unname(fp$estimate["shape"]), xi_pwm = pwm$results[["shape"]], sigma_pwm = pwm$results[["scale"]],
    xi_profile_lo = pr[1], xi_profile_hi = pr[3])
}
res <- do.call(rbind, rows)
py <- read.csv(data_dir("large_claims_threshold.csv"))
res$xi_python <- py$xi[match(res$quantile, py$quantile)]
res$sigma_python <- py$sigma[match(res$quantile, py$quantile)]
res$ll_python <- mapply(function(q, xi, sg) { u <- unname(quantile(x, q)); sum(evd::dgpd(x[x > u] - u, loc = 0, scale = sg, shape = xi, log = TRUE)) }, res$quantile, res$xi_python, res$sigma_python)
res$diff_mle_python <- res$xi_mle - res$xi_python
res$diff_own_python <- res$xi_mle_own - res$xi_python
res$diff_evd_python <- res$xi_evd_fpot - res$xi_python
write.csv(res, res_dir("C_gpd.csv"), row.names = FALSE)
print(res)
