# Do the Python reserving, pricing and tail models agree with the R reference packages? (v0.1)

[![tests](https://github.com/michaelserikov/r-validation-actuarial/actions/workflows/tests.yml/badge.svg)](https://github.com/michaelserikov/r-validation-actuarial/actions/workflows/tests.yml)

Three earlier projects fitted standard insurance models in Python: chain ladder, Mack and over-dispersed Poisson bootstrap on the CAS Loss Reserving Database ([loss-reserving](https://github.com/michaelserikov/loss-reserving)), a Poisson GLM against gradient boosting on the French motor liability data freMTPL2 ([motor-pricing-glm](https://github.com/michaelserikov/motor-pricing-glm)), and generalised Pareto fits to large claims in the same project. This project repeats them in R with the established packages and writes down where the two implementations agree and where they do not. It is an independent re-implementation check on open data, not an actuarial qualification or a pricing or reserving opinion.

**Main findings**
- **Reserving agrees to the sixth digit.** On all 229 usable company triangles of four lines of business the chain ladder reserve from `ChainLadder::MackChainLadder` equals the Python reserve (largest relative difference 4.7e-6, the rounding of the stored Python results) and Mack's standard error agrees (largest relative difference 4.6e-4). The Taylor and Ashe triangle gives the published 18,680,856 and 2,447,095 in both languages. The 90% interval coverage on the lower triangle is identical in both languages for Mack (53% to 79% by line) and for the bootstrap (57% to 79%).
- **The bootstrap agrees on average, not triangle by triangle.** The median relative difference of the bootstrap standard deviation is below 1% in every line, and 207 of 229 triangles lie within 20%. The 22 other triangles are small reserves or triangles with negative increments, where resampling residuals is unstable: the two random number streams and the treatment of negative increments give standard deviations that differ by a factor of two or more (listed in `results/A_reserving_by_triangle.csv`).
- **The Poisson GLM is reproduced, a GAM with vehicle age bands improves on it.** On the Python train and test split (542,410 and 135,603 policies) the R GLM has a Poisson deviance of 0.59898 and a Gini of 0.2792, the Python GLM 0.59901 and 0.2788. A GAM with smooth driver age and bonus-malus effects and vehicle age bands reaches 0.59689 and 0.2868. It stays well below gradient boosting (0.57674, Gini 0.3389). With a smooth for vehicle age instead of bands the GAM was worse than the GLM (0.6047, Gini 0.2681), because the effect of new vehicles is a step that a smooth cannot follow.
- **Default optimisers disagree on the tail, and the likelihood shows which one is right.** For the generalised Pareto fit to large claims over four thresholds (212 to 2,114 exceedances) the default maximum likelihood in `extRemes` stops 0.7 to 4.1 log-likelihood units below the optimum and `evd::fpot` 66 to 1,036 units below, with shape estimates of 0.84 to 0.95 (`extRemes`) and 0.34 to 0.40 (`fpot`) against 0.88 to 1.07. An own maximum likelihood in R (Nelder-Mead then BFGS) reproduces the scipy estimate of the Python project to 4e-5 in the shape. The likelihood is flat along a ridge of shape and scale, so the profile likelihood interval is wide (shape 0.85 to 1.08 at the 95% threshold) and the estimates from probability weighted moments (0.73 to 0.82) lie below the maximum likelihood values, inside the interval only at the 0.90 and 0.99 thresholds.

## Results

**Reserving** (`R/A_reserving.R`, `results/A_reserving_summary.csv`): R and Python on the 229 triangles of the backtest (1,000 bootstrap draws each)

| Line | Triangles | Largest difference in reserve | Largest difference in Mack s.e. | Median difference in bootstrap s.d. | Mack 90% coverage R / Python | Bootstrap 90% coverage R / Python |
|---|---|---|---|---|---|---|
| Private passenger auto | 70 | 4e-6 | 4e-6 | +0.7% | 60% / 60% | 57% / 57% |
| Commercial auto | 60 | 5e-6 | 5e-4 | -0.3% | 70% / 70% | 75% / 75% |
| Workers' compensation | 47 | 4e-6 | 4e-6 | +0.4% | 53% / 53% | 57% / 57% |
| Other liability | 52 | 4e-6 | 3e-6 | -0.1% | 79% / 79% | 79% / 79% |

**Frequency** (`R/B_pricing.R`, `results/B_frequency.csv`; test set, claims per year of exposure)

| Model | Poisson deviance | Explained deviance | Gini |
|---|---|---|---|
| Poisson GLM, R | 0.59898 | 5.02% | 0.2792 |
| Poisson GLM, Python | 0.59901 | 5.01% | 0.2788 |
| GAM, R (`mgcv`, smooth driver age and bonus-malus, vehicle age bands) | 0.59689 | 5.35% | 0.2868 |
| Gradient boosting, Python | 0.57674 | 8.54% | 0.3389 |

**Tail** (`R/C_tails.R`, `results/C_gpd.csv`): generalised Pareto fit to the individual claims of the training policies, shape and log-likelihood by threshold quantile

| Threshold quantile | Exceedances | Shape, scipy (Python) | Shape, own ML in R | Shape, `extRemes` ML | Shape, `evd::fpot` | Profile 95% interval (own) | Log-likelihood scipy / own / `extRemes` / `fpot` |
|---|---|---|---|---|---|---|---|
| 0.900 | 2,114 | 0.884 | 0.884 | 0.858 | 0.371 | 0.81 to 0.97 | -20,112.9 / -20,112.9 / -20,113.8 / -21,148.6 |
| 0.950 | 1,057 | 0.961 | 0.961 | 0.900 | 0.398 | 0.85 to 1.08 | -10,666.7 / -10,666.7 / -10,669.0 / -11,196.7 |
| 0.975 | 529 | 1.069 | 1.069 | 0.951 | 0.399 | 0.91 to 1.26 | -5,673.7 / -5,673.7 / -5,677.8 / -5,928.9 |
| 0.990 | 212 | 0.936 | 0.936 | 0.842 | 0.340 | 0.69 to 1.25 | -2,488.7 / -2,488.7 / -2,489.5 / -2,554.3 |

## Method

- **Same data.** Python writes the exact policies, split and claims of the earlier projects to `data/` (`scripts/export_from_python.py`); the CAS triangles and the stored Python results are copies from the loss-reserving repository. The triangle filter (net premium of at least 0.5 million USD in every accident year, positive first column and positive paid claims in the last) is the Python rule re-implemented.
- **Reserving.** `MackChainLadder` with volume-weighted factors and Mack's variance extrapolation, `BootChainLadder` with the over-dispersed Poisson process distribution, 1,000 draws, seed set to the company code. The Python bootstrap uses its own random number generator, so only moments and coverage are comparable.
- **Frequency.** Poisson with the log of exposure as offset. GLM with the same tariff cells as Python (age bands, vehicle age bands, capped power, area, brand, fuel, region, bonus-malus and its logarithm, log density) and no penalty. GAM by `bam` with fast REML and exact covariates.
- **Tail.** Threshold at four quantiles of the individual claims, exceedances above it. Maximum likelihood in `extRemes::fevd`, `evd::fpot` and an own optimiser; probability weighted moments (Hosking and Wallis) implemented directly; profile likelihood interval from a grid of the shape.
- **Tolerances used in the statements above.** Reserve and Mack standard error: relative difference below 1e-3 in every triangle. Bootstrap: only the median difference of the standard deviation below 2% per line and the coverage within two triangles. Frequency: deviance within 1e-4 and Gini within 2e-3 of Python for the GLM. Tail: shape within 1e-3 of scipy and log-likelihood not below scipy by more than 0.01 for the own optimiser.

## Limits

- One run per module, one random split (seed 7), one bootstrap seed. The bootstrap comparison says nothing about which implementation is closer to the truth.
- The GAM specification was changed once after the first fit (vehicle age as bands); both fits are reported above.
- The profile interval uses a grid of 0.005 and a one-dimensional optimiser for the scale.
- No panel regression (fixest) module in this version, and no actuarial qualification behind these checks.
- R 4.3.3 with `ChainLadder`, `mgcv`, `extRemes`, `evd` at the versions of the Posit package manager snapshot of 24 November 2026; no `renv.lock` is committed yet.

## Reproduce

```
python scripts/export_from_python.py        # needs the two Python projects next to this folder
Rscript R/A_reserving.R                      # about 9 minutes
Rscript R/B_pricing.R                        # about 3 to 12 minutes for the GAM
Rscript R/C_tails.R
Rscript -e 'testthat::test_dir("tests")'
```
