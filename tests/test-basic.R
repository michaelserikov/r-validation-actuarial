library(testthat)
library(ChainLadder)
source("../R/common.R")

test_that("Taylor-Ashe chain ladder reserve and Mack standard error equal the published values", {
  ta <- MackChainLadder(GenIns, est.sigma = "Mack")
  latest <- apply(GenIns, 1, function(r) tail(r[!is.na(r)], 1))
  reserve <- sum(ta$FullTriangle[, 10] - latest)
  expect_equal(reserve, 18680856, tolerance = 1e-6)
  expect_equal(unname(ta$Total.Mack.S.E), 2447095, tolerance = 5e-4)
})

test_that("Gini of a perfect ranking is positive and of a constant prediction is zero", {
  set.seed(1)
  y <- rpois(2000, 0.1); w <- rep(1, 2000)
  expect_lt(abs(lorenz_gini(y, rep(1, 2000), w)), 0.01)
  expect_gt(lorenz_gini(y, y + runif(2000, 0, 0.1), w), 0.5)
})

test_that("Poisson deviance is zero at the observation and positive otherwise", {
  y <- c(0, 1, 2, 3); w <- rep(1, 4)
  expect_equal(poisson_dev(y[-1], y[-1], w[-1]), 0)
  expect_gt(poisson_dev(y[-1], rep(2, 3), w[-1]), 0)
})
