# Module B: motor pricing on freMTPL2 with the exact split of the Python project. Poisson GLM and GAM (mgcv) against Python numbers.
suppressPackageStartupMessages({library(mgcv); library(data.table)})
source("R/common.R")
pol <- fread(data_dir("fremtpl2_policies.csv"))
pol[, `:=`(Area = factor(Area), VehBrand = factor(VehBrand), VehGas = factor(VehGas), Region = factor(Region))]
pol[, ClaimNb := pmin(ClaimNb, 4)]; pol[, Exposure := pmin(Exposure, 1)]
cut_f <- function(x, br) cut(x, br, right = TRUE)
pol[, DrivAgeBand := cut_f(DrivAge, c(17, 20, 25, 30, 40, 50, 60, 70, 120))]
pol[, VehAgeBand := cut_f(VehAge, c(-1, 0, 1, 5, 10, 15, 120))]
pol[, VehPowerCap := factor(pmin(VehPower, 9))]
pol[, BM := pmin(BonusMalus, 150)]; pol[, LogBM := log(BM)]; pol[, LogDens := log(Density)]
tr <- pol[split == "train"]; te <- pol[split == "test"]
ytr <- tr$ClaimNb / tr$Exposure; yte <- te$ClaimNb / te$Exposure
cat("train", nrow(tr), "test", nrow(te), "\n")

# GLM: same tariff cells as the Python pipeline (bands, categoricals, bonus-malus, log density), Poisson, exposure weights
f_glm <- ClaimNb ~ DrivAgeBand + VehAgeBand + VehPowerCap + Area + VehBrand + VehGas + Region + BM + LogBM + LogDens + offset(log(Exposure))
glm_fit <- glm(f_glm, family = poisson(), data = tr)
p_glm <- predict(glm_fit, te, type = "response") / te$Exposure

# GAM: smooth effects of driver age, vehicle age, bonus-malus, density and power; categoricals as in the GLM (bam, fast REML, exact covariates)
f_gam <- ClaimNb ~ s(DrivAge, k = 20) + VehAgeBand + s(BM, k = 12) + s(LogDens, k = 6) + VehPowerCap + Area + VehBrand + VehGas + Region + offset(log(Exposure))
t0 <- Sys.time()
gam_fit <- bam(f_gam, family = poisson(), data = as.data.frame(tr), method = "fREML", discrete = FALSE, nthreads = 4)
cat("gam fit seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "\n")
p_gam <- predict(gam_fit, as.data.frame(te), type = "response") / te$Exposure

w <- te$Exposure
null <- weighted.mean(ytr, tr$Exposure)
dev0 <- poisson_dev(yte, rep(null, nrow(te)), w)
row <- function(name, p) data.frame(model = name, poisson_deviance = poisson_dev(yte, p, w), d2 = 1 - poisson_dev(yte, p, w) / dev0,
                                    gini = lorenz_gini(yte, p, w), predicted_claims = sum(p * w), observed_claims = sum(te$ClaimNb))
out <- rbind(row("mean_R", rep(null, nrow(te))), row("glm_R", p_glm), row("gam_R", p_gam))
out$gini[out$model == "mean_R"] <- NA   # a constant prediction has no ranking (ties)
py <- read.csv(data_dir("frequency.csv")); py$model <- paste0(py$model, "_python")
py <- py[, c("model", "poisson_deviance", "d2", "gini", "predicted_claims", "observed_claims")]
res <- rbind(out, py)
write.csv(res, res_dir("B_frequency.csv"), row.names = FALSE)
print(res)
# agreement of the two GLMs is limited by the penalty (alpha 1e-4) and scaling in scikit-learn: report the rank correlation too
cat("GLM vs GAM Spearman:", cor(p_glm, p_gam, method = "spearman"), "\n")
sm <- summary(gam_fit)$s.table
write.csv(as.data.frame(sm), res_dir("B_gam_smooths.csv"))
saveRDS(list(p_glm = p_glm, p_gam = p_gam), res_dir("B_predictions.rds"))
