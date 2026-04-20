################################################################################
## Empirical inference for the 1974-2023 temperature data
## Estimates the SFSAR model with functional rho(t) and applies
##
## Optional environment overrides:
##   APP_NCORES, APP_K, APP_M_EST, APP_M_BOOT, APP_L_BOOT,
##   APP_B_LM, APP_B_F, APP_B_CI, APP_ALPHA, APP_TOL,
##   APP_CI_METHOD, APP_CI_MULTIPLIER, APP_F_MULTIPLIER,
##   APP_LM_MULTIPLIER, APP_RHO_CLIP
################################################################################

library(fda)
library(parallel)

get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, winslash = "/")))
  }
  normalizePath(getwd(), winslash = "/")
}

SCRIPT_DIR <- get_script_dir()
ROOT_DIR   <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
options(sfsar_root_dir = ROOT_DIR)
setwd(SCRIPT_DIR)

source(file.path(ROOT_DIR, "bootstrap_inference.R"))

get_env_integer <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(as.integer(default))
  }
  value <- suppressWarnings(as.integer(raw))
  if (is.na(value)) {
    stop("Environment variable ", name, " must be an integer.")
  }
  value
}

get_env_numeric <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(as.numeric(default))
  }
  value <- suppressWarnings(as.numeric(raw))
  if (!is.finite(value)) {
    stop("Environment variable ", name, " must be numeric.")
  }
  value
}

expand_ylim <- function(values, pad = 0.06) {
  rng <- range(values, na.rm = TRUE)
  if (!all(is.finite(rng))) {
    stop("Cannot build plot limits from non-finite values.")
  }
  span <- diff(rng)
  if (span <= 0) {
    base <- max(1, abs(rng[1]) * 0.1)
    return(rng + c(-base, base))
  }
  rng + c(-pad, pad) * span
}

plot_curve_with_ci <- function(grids, ci_obj, main, ylab, ylim = NULL,
                               ci_col = "blue") {
  month_days  <- c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)
  month_ticks <- cumsum(month_days)[c(3, 6, 9, 12)]
  est         <- ci_obj$estimate 
  low         <- ci_obj$lower 
  upp         <- ci_obj$upper 
  if (is.null(ylim)) {
    ylim <- expand_ylim(c(est, low, upp))
  }

  plot(
    grids, est, type = "l", lwd = 1.8, col = "black",
    main = main, xlab = "", ylab = ylab,
    xlim = c(0, 365), ylim = ylim, xaxt = "n"
  )
  axis(
    1,
    at = month_ticks,
    labels = c("Mar", "Jun", "Sep", "Dec")
  )
  abline(v = month_ticks, col = "gray", lty = 3)
  lines(grids, low, lty = 2, lwd = 1.3, col = ci_col)
  lines(grids, upp, lty = 2, lwd = 1.3, col = ci_col)
}

################################################################################
## 1. Analysis configuration
################################################################################

config <- list(
  K         = get_env_integer("APP_K", 13L),
  M_est     = get_env_integer("APP_M_EST", 2920L),
  M_boot    = get_env_integer("APP_M_BOOT", 500L),
  L_boot    = get_env_integer("APP_L_BOOT", 500L),
  B_lm      = get_env_integer("APP_B_LM", 199L),
  B_f       = get_env_integer("APP_B_F", 199L),
  B_ci      = get_env_integer("APP_B_CI", 199L),
  alpha     = get_env_numeric("APP_ALPHA", 0.05),
  tolerance = get_env_numeric("APP_TOL", 1e-4),
  rho_clip  = get_env_numeric("APP_RHO_CLIP", 0.99),
  ncores    = get_env_integer(
    "APP_NCORES",
    if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)
  )
)

ci_method     <- Sys.getenv("APP_CI_METHOD", unset = "wald")
ci_method     <- match.arg(ci_method, c("wald", "normal", "bc", "percentile", "basic"))
ci_multiplier <- Sys.getenv("APP_CI_MULTIPLIER", unset = "gaussian")
ci_multiplier <- match.arg(ci_multiplier, c("gaussian", "rademacher", "mammen"))
f_multiplier  <- Sys.getenv("APP_F_MULTIPLIER", unset = "mammen")
f_multiplier  <- match.arg(f_multiplier, c("mammen", "rademacher", "gaussian"))
lm_multiplier <- Sys.getenv("APP_LM_MULTIPLIER", unset = "mammen")
lm_multiplier <- match.arg(lm_multiplier, c("mammen", "rademacher", "gaussian"))

config$ci_method     <- ci_method
config$ci_multiplier <- ci_multiplier
config$f_multiplier  <- f_multiplier
config$lm_multiplier <- lm_multiplier

################################################################################
## 2. Load data and build the empirical design
################################################################################

load("data_1974to2023_3h.RData")
load("sites.RData")
load("weightmat_inv.RData")

data <- data_1974to2023_3h
temp <- array(data, c(dim(data)[1], dim(data)[2] * dim(data)[3]))
dimnames(temp)[[2]] <- c(t(outer(dimnames(data)[[3]], dimnames(data)[[2]], paste, sep = "_")))

tbasis  <- create.bspline.basis(c(0, 365), config$K)
seq_obs <- seq(from = 0, by = 1 / 8, length.out = 2920)
N <- ncol(temp)

station_id <- sub("_.*", "", colnames(temp))
obs_year   <- sub(".*_", "", colnames(temp))

intercept <- rep(1, N)
elevation <- sites$elevation[match(station_id, sites$id)]
if (anyNA(elevation)) {
  stop("Missing elevation values after matching station ids.")
}

year_levels    <- as.character(1974:2023)
year           <- matrix(0, nrow = N, ncol = length(year_levels))
colnames(year) <- year_levels
for (i in seq_along(year_levels)) {
  year[, i] <- as.integer(obs_year == year_levels[i])
}

y      <- t(temp)
W      <- w_matrix / sum(w_matrix) * nrow(w_matrix)
w_year <- colSums(year) / sum(year)

constraint_row <- matrix(c(0, 0, w_year), nrow = 1)
z1             <- cbind(intercept, elevation, year)
z1             <- rbind(z1, constraint_row)
y1             <- rbind(y, rep(0, ncol(y)))
W1             <- matrix(0, nrow = nrow(y1), ncol = nrow(y1))
W1[seq_len(N), seq_len(N)] <- W
p1             <- ncol(z1)

################################################################################
## 3. Re-estimate the full SFSAR model
################################################################################

cat("===== Full-model estimation =====\n")
fit_full <- fsar_func_fast(
  y = y1, Z = z1, W = W1, seq = seq_obs,
  basis = tbasis, M = config$M_est, tolerance = config$tolerance
)

D_hat   <- matrix(fit_full$D, nrow = p1, ncol = config$K, byrow = FALSE)
r_hat   <- matrix(fit_full$r, nrow = 1)
eta_hat <- matrix(fit_full$eta, nrow = 1)

cat("Converged:", fit_full$converged, "\n")
cat("Iterations:", fit_full$iter, "\n")
cat("beta diff:", fit_full$beta_dif, "\n")
cat("rho diff:", fit_full$rho_dif, "\n")

################################################################################
## 4. LM test for H0: rho(t) = 0
################################################################################

cat("===== LM test for rho(t) =====\n")
lm_result <- bootstrap_LM_test(
  y = y1, Z = z1, W = W1, seq_obs = seq_obs,
  basis = tbasis, M = config$M_est, B = config$B_lm,
  ncores = config$ncores, multiplier = config$lm_multiplier
)

cat("LM statistic:", lm_result$LM, "\n")
cat("Bootstrap p-value:", lm_result$p_value, "\n")

################################################################################
## 5. F-type tests
################################################################################

cat("===== F-test: intercept =====\n")
ftest_intercept <- bootstrap_Ftest(
  y = y1, Z = z1, W = W1, seq_obs = seq_obs,
  basis = tbasis, M = config$M_est, B = config$B_f, j_test = 1,
  D_FM       = D_hat, r_FM = r_hat, eta_FM = eta_hat,
  ncores     = config$ncores, L_boot = config$L_boot, M_boot = config$M_boot,
  multiplier = config$f_multiplier
)
cat("F statistic (intercept):", ftest_intercept$Fn, "\n")
cat("Bootstrap p-value:", ftest_intercept$p_value, "\n")

cat("===== F-test: elevation =====\n")
ftest_elevation <- bootstrap_Ftest(
  y = y1, Z = z1, W = W1, seq_obs = seq_obs,
  basis = tbasis, M = config$M_est, B = config$B_f, j_test = 2,
  D_FM = D_hat, r_FM = r_hat, eta_FM = eta_hat,
  ncores = config$ncores, L_boot = config$L_boot, M_boot = config$M_boot,
  multiplier = config$f_multiplier
)
cat("F statistic (elevation):", ftest_elevation$Fn, "\n")
cat("Bootstrap p-value:", ftest_elevation$p_value, "\n")

cat("===== F-test: year effects =====\n")
ftest_year <- bootstrap_Ftest(
  y = y1, Z = z1, W = W1, seq_obs = seq_obs,
  basis = tbasis, M = config$M_est, B = config$B_f, j_test = 3:p1,
  D_FM = D_hat, r_FM = r_hat, eta_FM = eta_hat,
  ncores = config$ncores, L_boot = config$L_boot, M_boot = config$M_boot,
  multiplier = config$f_multiplier
)
cat("F statistic (year effects):", ftest_year$Fn, "\n")
cat("Bootstrap p-value:", ftest_year$p_value, "\n")

################################################################################
## 6. Pointwise confidence intervals for intercept, elevation, and rho(t)
################################################################################

cat("===== Pointwise confidence intervals =====\n")
ci_result <- bootstrap_CI(
  y = y1, Z = z1, W = W1, seq_obs = seq_obs,
  basis = tbasis, M = config$M_est, B = config$B_ci, alpha = config$alpha,
  D_hat = D_hat, r_hat = r_hat, eta_hat = eta_hat,
  ncores = config$ncores, L_boot = config$L_boot, M_boot = config$M_boot,
  ci_method = config$ci_method, multiplier = config$ci_multiplier,
  rho_clip = config$rho_clip
)

################################################################################
## 7. Save numerical results
################################################################################

summary_table <- data.frame(
  test = c("LM_rho", "F_intercept", "F_elevation", "F_year"),
  statistic = c(
    lm_result$LM,
    ftest_intercept$Fn,
    ftest_elevation$Fn,
    ftest_year$Fn
  ),
  p_value = c(
    lm_result$p_value,
    ftest_intercept$p_value,
    ftest_elevation$p_value,
    ftest_year$p_value
  ),
  B_valid = c(
    lm_result$B_valid,
    ftest_intercept$B_valid,
    ftest_elevation$B_valid,
    ftest_year$B_valid
  )
)

results <- list(
  config = config,
  fit_full = fit_full,
  D_hat = D_hat,
  r_hat = r_hat,
  eta_hat = eta_hat,
  lm_result = lm_result,
  ftest_intercept = ftest_intercept,
  ftest_elevation = ftest_elevation,
  ftest_year = ftest_year,
  ci_result = ci_result,
  summary_table = summary_table
)

save(results, file = "empirical_inference.RData")
write.csv(summary_table, file = "empirical_test_summary.csv", row.names = FALSE)

################################################################################
## 8. Plot estimated curves with blue dashed confidence intervals
################################################################################

eval_grids <- ci_result$eval_grids

pdf("fig_ci.pdf", width = 12, height = 4)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))
plot_curve_with_ci(
  eval_grids, ci_result$beta_ci[[1]],
  main = "intercept",
  ylab = expression("Temperature (" * degree * "C)")
)
plot_curve_with_ci(
  eval_grids, ci_result$beta_ci[[2]],
  main = "elevation effect",
  ylab = expression("Temperature (" * degree * "C per 100 meters)"),
  scale = 100
)
plot_curve_with_ci(
  eval_grids, ci_result$rho_ci,
  main = "spatial autocorrelation coefficient",
  ylab = expression(rho(t))
)
dev.off()

png("fig_ci.png", width = 1800, height = 600, res = 200)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))
plot_curve_with_ci(
  eval_grids, ci_result$beta_ci[[1]],
  main = "intercept",
  ylab = expression("Temperature (" * degree * "C)")
)
plot_curve_with_ci(
  eval_grids, ci_result$beta_ci[[2]],
  main = "elevation effect",
  ylab = expression("Temperature (" * degree * "C per 100 meters)"),
  scale = 100
)
plot_curve_with_ci(
  eval_grids, ci_result$rho_ci,
  main = "spatial autocorrelation coefficient",
  ylab = expression(rho(t))
)
dev.off()

cat("\n===== Summary =====\n")
print(summary_table)
cat("Saved: empirical_inference.RData\n")
cat("Saved: empirical_test_summary.csv\n")
cat("Saved: fig_ci.pdf\n")
cat("Saved: fig_ci.png\n")

