# Module A: reserving. Chain ladder, Mack and ODP bootstrap in ChainLadder against the Python implementation (loss-reserving).
suppressPackageStartupMessages({library(ChainLadder)})
source("R/common.R")
set.seed(1)
N <- 10
lines <- c(ppauto = "Private passenger auto", comauto = "Commercial auto", wkcomp = "Workers' compensation", othliab = "Other liability (occurrence)")
n_boot <- as.integer(Sys.getenv("N_BOOT", "1000"))

lognormal_q <- function(mean, sd, q) { if (mean <= 0 || sd <= 0) return(mean); s2 <- log(1 + (sd / mean)^2); mean * exp(-0.5 * s2 + sqrt(s2) * qnorm(q)) }

load_line <- function(line) {
  d <- read.csv(data_dir("cas", paste0(line, "_pos.csv")))
  suf <- sub("CumPaidLoss_", "", grep("^CumPaidLoss_", names(d), value = TRUE)[1])
  pc <- paste0("CumPaidLoss_", suf); prem <- paste0("EarnedPremNet_", suf)
  out <- list()
  for (code in unique(d$GRCODE)) {
    g <- d[d$GRCODE == code, ]
    if (nrow(g) != N * N) next
    g <- g[order(g$AccidentYear, g$DevelopmentLag), ]
    paid <- matrix(g[[pc]], N, N, byrow = TRUE)
    pr <- tapply(g[[prem]], g$AccidentYear, function(x) x[1])
    out[[as.character(code)]] <- list(company = code, line = line, paid = paid, premium = as.numeric(pr))
  }
  out
}
upper <- function(m) { u <- m; for (i in 1:N) if (i > 1) u[i, (N - i + 2):N] <- NA; u }
usable <- function(t) { up <- upper(t$paid); all(t$premium >= 500) && all(up[, 1] > 0) && all(t$paid[, N] > 0) }

rows <- list()
for (ln in names(lines)) for (t in load_line(ln)) {
  if (!usable(t)) next
  tri <- as.triangle(upper(t$paid), origin = "ay", dev = "dev")
  mk <- try(MackChainLadder(tri, weights = 1, alpha = 1, est.sigma = "Mack", tail = FALSE), silent = TRUE)
  if (inherits(mk, "try-error")) next
  full <- mk$FullTriangle
  latest <- apply(tri, 1, function(r) tail(r[!is.na(r)], 1))
  res_by <- full[, N] - latest
  reserve <- sum(res_by[2:N]); if (!is.finite(reserve) || reserve <= 0) next
  se <- mk$Total.Mack.S.E
  b <- BootChainLadder(tri, R = n_boot, process.distr = "od.pois", seed = as.integer(t$company))
  bs <- b$IBNR.Totals
  q <- quantile(bs, c(0.05, 0.25, 0.5, 0.75, 0.95, 0.995))
  actual <- sum((t$paid[, N] - latest)[2:N])
  rows[[length(rows) + 1]] <- data.frame(line = ln, company = t$company, r_cl_paid = reserve, r_mack_se = se, r_boot_mean = mean(bs), r_boot_sd = sd(bs),
    r_boot_q05 = q[1], r_boot_q95 = q[5], r_boot_q995 = q[6], actual_reserve = actual,
    r_mack_q05 = lognormal_q(reserve, se, 0.05), r_mack_q95 = lognormal_q(reserve, se, 0.95))
}
r <- do.call(rbind, rows)
py <- read.csv(data_dir("python_backtest_triangles.csv"))
m <- merge(r, py[, c("line", "company", "pred_cl_paid", "mack_se", "boot_mean", "boot_sd", "boot_q05", "boot_q95", "mack_q05", "mack_q95", "actual_reserve")], by = c("line", "company"), suffixes = c("", "_py"))
m$d_reserve <- m$r_cl_paid / m$pred_cl_paid - 1
m$d_mack_se <- m$r_mack_se / m$mack_se - 1
m$d_boot_sd <- m$r_boot_sd / m$boot_sd - 1
m$d_boot_mean <- m$r_boot_mean / m$boot_mean - 1
write.csv(m, res_dir("A_reserving_by_triangle.csv"), row.names = FALSE)

cov <- function(lo, hi, act) mean(act >= lo & act <= hi)
summ <- do.call(rbind, lapply(split(m, m$line), function(d) data.frame(
  line = d$line[1], triangles_R = nrow(d), triangles_Python = sum(py$line == d$line[1]),
  max_abs_diff_reserve = max(abs(d$d_reserve)), max_abs_diff_mack_se = max(abs(d$d_mack_se)),
  median_diff_boot_sd = median(d$d_boot_sd), p90_abs_diff_boot_sd = quantile(abs(d$d_boot_sd), 0.9),
  mack90_coverage_R = cov(d$r_mack_q05, d$r_mack_q95, d$actual_reserve), mack90_coverage_py = cov(d$mack_q05, d$mack_q95, d$actual_reserve_py %||% d$actual_reserve),
  boot90_coverage_R = cov(d$r_boot_q05, d$r_boot_q95, d$actual_reserve), boot90_coverage_py = cov(d$boot_q05, d$boot_q95, d$actual_reserve))))
write.csv(summ, res_dir("A_reserving_summary.csv"), row.names = FALSE)
print(summ)

# Taylor-Ashe: published reserve 18,680,856, Mack total s.e. 2,447,095
ta <- MackChainLadder(GenIns, est.sigma = "Mack")
ta_out <- data.frame(item = c("reserve", "mack_se_total"), R = c(sum(ta$FullTriangle[, 10] - apply(GenIns, 1, function(r) tail(r[!is.na(r)], 1))), ta$Total.Mack.S.E), published = c(18680856, 2447095))
write.csv(ta_out, res_dir("A_taylor_ashe.csv"), row.names = FALSE)
print(ta_out)
