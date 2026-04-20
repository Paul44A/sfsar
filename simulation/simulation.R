################################################################################
## simulation.R
##
##   1. Estimation accuracy of FoS, SFSAR with scalar rho, SFSAR with rho(t)
##   2. LM test, F-type test, and pointwise CI for SFSAR with rho(t)
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(fda)
  library(Matrix)
  library(parallel)
})

get_script_dir <- function() {
  cmd_args     <- commandArgs(trailingOnly = FALSE)
  file_arg     <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
  }
  normalizePath(getwd())
}

get_env_integer <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(as.integer(default))
  }
  as.integer(raw)
}

get_env_numeric <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(as.numeric(default))
  }
  as.numeric(raw)
}

get_env_char <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(default)
  }
  raw
}

get_env_vector <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) {
    return(default)
  }
  trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
}

SCRIPT_DIR <- get_script_dir()
ROOT_DIR   <- normalizePath(file.path(SCRIPT_DIR, ".."), winslash = "/", mustWork = TRUE)
options(sfsar_root_dir = ROOT_DIR)
setwd(SCRIPT_DIR)

source(file.path(ROOT_DIR, "fos.R"))
source(file.path(ROOT_DIR, "fsar_func.R"))
source(file.path(ROOT_DIR, "fsar_scalar.R"))
source(file.path(ROOT_DIR, "y_gen.R"))
source(file.path(ROOT_DIR, "bootstrap_inference.R"))

OUT_DIR <- file.path(SCRIPT_DIR, "results")
CSV_DIR <- file.path(OUT_DIR, "csv")
RDS_DIR <- file.path(OUT_DIR, "rds")
CHK_DIR <- file.path(OUT_DIR, "checkpoints")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CSV_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(CHK_DIR, showWarnings = FALSE, recursive = TRUE)

CONFIG <- list(
  seed          = get_env_integer("SIM_SEED", 20260313L),
  n_sim         = get_env_integer("SIM_NSIM", 200L),
  B_boot        = get_env_integer("SIM_BBOOT", 199L),
  ncores        = get_env_integer("SIM_NCORES", 24L),
  K             = get_env_integer("SIM_K", 9L),
  L_obs         = get_env_integer("SIM_LOBS", 101L),
  M_est         = get_env_integer("SIM_MEST", 101L),
  M_boot        = get_env_integer("SIM_MBOOT", 101L),
  alpha         = get_env_numeric("SIM_ALPHA", 0.05),
  rho_clip      = get_env_numeric("SIM_RHO_CLIP", 0.99),
  tolerance     = get_env_numeric("SIM_TOL", 1e-4),
  ci_method     = get_env_char("SIM_CI_METHOD", "wald"),
  ci_multiplier = get_env_char("SIM_CI_MULTIPLIER", "mammen"),
  stages        = get_env_vector("SIM_STAGES", c("estimation", "lm", "f", "ci")),
  N_list        = as.integer(get_env_vector("SIM_N_LIST", c("100", "400", "1600"))),
  rho_labels    = get_env_vector("SIM_RHO_LIST", c("rho0", "rho_const075", "rho_sin075")),
  f_beta_cases  = get_env_vector("SIM_F_BETA_CASES", c("null", "alt")),
  eval_grid_n   = get_env_integer("SIM_EVAL_GRID_N", 401L),
  method_ncores = 1L,
  j_test = 4L
)

CONFIG$seq_obs      <- seq(0, 1, length.out = CONFIG$L_obs)
CONFIG$eval_grid    <- seq(0, 1, length.out = CONFIG$eval_grid_n)
CONFIG$beta_names   <- c("beta0", "beta1", "beta2", "beta3")
CONFIG$method_names <- c("fos", "sfsar_rho", "sfsar_rhot")
CONFIG$side_map     <- c("100" = 10L, "400" = 20L, "1600" = 40L)

W_CACHE <- new.env(parent = emptyenv())

make_basis <- function() {
  create.bspline.basis(c(0, 1), CONFIG$K)
}

make_rook_torus_W <- function(N) {
  key <- as.character(N)
  if (exists(key, envir = W_CACHE, inherits = FALSE)) {
    return(get(key, envir = W_CACHE, inherits = FALSE))
  }

  side <- CONFIG$side_map[[key]]
  if (is.null(side) || !is.finite(side)) {
    stop("Unsupported N for torus rook graph: ", N)
  }
  side <- as.integer(side)
  idx <- function(r, c) ((r %% side) * side) + (c %% side) + 1L
  W <- matrix(0, N, N)
  for (r in 0:(side - 1L)) {
    for (c in 0:(side - 1L)) {
      cur <- idx(r, c)
      nbrs <- c(
        idx(r - 1L, c),
        idx(r + 1L, c),
        idx(r, c - 1L),
        idx(r, c + 1L)
      )
      W[cur, nbrs] <- 1
    }
  }
  W <- W / rowSums(W)
  assign(key, W, envir = W_CACHE)
  W
}

scale_curve <- function(x, amplitude = 1) {
  amplitude * x / max(abs(x))
}

make_beta_curves <- function(seq_grid, beta_case = c("alt", "null")) {
  beta_case <- match.arg(beta_case)
  beta0     <- 0.4 + scale_curve(dnorm(seq_grid, mean = 0.50, sd = 0.18), amplitude = 0.9)
  beta1     <- 0.6 * sin(2 * pi * seq_grid)
  beta2     <- 0.7 * cos(pi * seq_grid)
  beta3_alt <- scale_curve(
    dnorm(seq_grid, mean = 0.25, sd = 0.10) - 0.8 * dnorm(seq_grid, mean = 0.75, sd = 0.10),
    amplitude = 0.5
  )
  beta3 <- if (beta_case == "alt") beta3_alt else rep(0, length(seq_grid))

  rbind(beta0, beta1, beta2, beta3)
}

make_rho_curve <- function(seq_grid, rho_label) {
  switch(
    rho_label,
    rho0 = rep(0, length(seq_grid)),
    rho_const075 = rep(0.75, length(seq_grid)),
    rho_sin075 = 0.75 * sin(2 * pi * seq_grid),
    stop("Unknown rho label: ", rho_label)
  )
}

curve_matrix_to_coef <- function(curve_mat, basis, seq_grid) {
  coef_list <- lapply(seq_len(nrow(curve_mat)), function(j) {
    smooth.basis(seq_grid, curve_mat[j, ], basis)$fd$coefs
  })
  t(do.call(cbind, coef_list))
}

curve_to_coef_row <- function(curve_vec, basis, seq_grid) {
  matrix(smooth.basis(seq_grid, curve_vec, basis)$fd$coefs, nrow = 1)
}

build_truth <- function(rho_label, beta_case = c("alt", "null")) {
  beta_case   <- match.arg(beta_case)
  basis       <- make_basis()
  beta_curves <- make_beta_curves(CONFIG$seq_obs, beta_case = beta_case)
  rho_curve   <- make_rho_curve(CONFIG$seq_obs, rho_label)
  list(
    basis            = basis,
    beta_curves_obs  = beta_curves,
    rho_curve_obs    = rho_curve,
    D_true           = curve_matrix_to_coef(beta_curves, basis, CONFIG$seq_obs),
    r_true           = curve_to_coef_row(rho_curve, basis, CONFIG$seq_obs),
    beta_curves_eval = make_beta_curves(CONFIG$eval_grid, beta_case = beta_case),
    rho_curve_eval   = make_rho_curve(CONFIG$eval_grid, rho_label)
  )
}

generate_design_matrix <- function(N) {
  x_bin       <- rbinom(N, size = 1, prob = 0.5)
  x_g1        <- as.numeric(scale(rnorm(N), center = TRUE, scale = TRUE))
  x_g2        <- as.numeric(scale(rnorm(N), center = TRUE, scale = TRUE))
  Z           <- cbind(1, x_bin, x_g1, x_g2)
  colnames(Z) <- c("intercept", "x_bin", "x_gauss1", "x_gauss2")
  Z
}

trapz_integral <- function(x, y) {
  sum((y[-1] + y[-length(y)]) * diff(x) / 2)
}

curve_metrics <- function(est, truth, grid) {
  sq_err <- (est - truth)^2
  c(
    imse = trapz_integral(grid, sq_err),
    rmse = sqrt(mean(sq_err))
  )
}

task_seed <- function(task, stage_code) {
  rho_index  <- match(task$rho_label, CONFIG$rho_labels)
  N_index    <- match(task$N, CONFIG$N_list)
  beta_index <- if (!is.null(task$beta_case)) match(task$beta_case, CONFIG$f_beta_cases) else 0L
  CONFIG$seed + stage_code * 1000000L + N_index * 10000L + rho_index * 1000L + beta_index * 100L + task$rep_id
}

extract_fit_curves <- function(method, fit_obj, basis, grid, p) {
  K <- basis$nbasis
  D_hat <- switch(
    method,
    fos = matrix(fit_obj[[1]], nrow = p, ncol = K, byrow = FALSE),
    sfsar_rho = matrix(fit_obj$D, nrow = p, ncol = K, byrow = FALSE),
    sfsar_rhot = matrix(fit_obj$D, nrow = p, ncol = K, byrow = FALSE),
    stop("Unknown method: ", method)
  )
  Phi <- eval.basis(grid, basis)
  beta_hat <- D_hat %*% t(Phi)

  rho_hat <- switch(
    method,
    fos        = rep(0, length(grid)),
    sfsar_rho  = rep(as.numeric(fit_obj$r)[1], length(grid)),
    sfsar_rhot = {
      if (!is.null(fit_obj$eta) && length(as.numeric(fit_obj$eta)) > 0) {
        eval_rho_curve(grid, basis, eta = as.numeric(fit_obj$eta), rho_bound = fit_obj$rho_bound)
      } else {
        eval_rho_curve(grid, basis, r = matrix(fit_obj$r, nrow = 1))
      }
    }
  )

  list(beta_hat = beta_hat, rho_hat = as.numeric(rho_hat))
}

fit_method <- function(method, y, Z, W, basis) {
  t_start <- proc.time()[["elapsed"]]
  fit_obj <- switch(
    method,
    fos = fos(y = y, Z = Z, seq = CONFIG$seq_obs, basis = basis,
              M = CONFIG$M_est, ncores = CONFIG$method_ncores),
    sfsar_rho = fsar_scalar(y = y, Z = Z, W = W, seq = CONFIG$seq_obs, basis = basis,
                            M = CONFIG$M_est, tolerance = CONFIG$tolerance,
                            ncores = CONFIG$method_ncores),
    sfsar_rhot = fsar_func_fast(y = y, Z = Z, W = W, seq = CONFIG$seq_obs, basis = basis,
                                M = CONFIG$M_est, tolerance = CONFIG$tolerance),
    stop("Unknown method: ", method)
  )
  elapsed <- proc.time()[["elapsed"]] - t_start
  list(fit = fit_obj, elapsed_sec = elapsed)
}

estimation_worker <- function(task) {
  set.seed(task_seed(task, stage_code = 1L))
  truth <- build_truth(task$rho_label, beta_case = "alt")
  W     <- make_rook_torus_W(task$N)
  Z     <- generate_design_matrix(task$N)
  y     <- y_gen(
    D = truth$D_true, r = truth$r_true, Z = Z, W = W, basis = truth$basis,
    mean = 0, sd = 0.25, L = CONFIG$L_obs
  )

  out_rows <- vector("list", length(CONFIG$method_names))
  for (m_idx in seq_along(CONFIG$method_names)) {
    method <- CONFIG$method_names[m_idx]
    fit_res <- tryCatch(fit_method(method, y, Z, W, truth$basis), error = function(e) e)
    if (inherits(fit_res, "error")) {
      out_rows[[m_idx]] <- data.frame(
        stage = "estimation",
        N                = task$N,
        rho_label        = task$rho_label,
        rep_id           = task$rep_id,
        method           = method,
        converged        = FALSE,
        elapsed_sec      = NA_real_,
        error_message    = conditionMessage(fit_res),
        stringsAsFactors = FALSE
      )
      next
    }

    curves <- extract_fit_curves(method, fit_res$fit, truth$basis, CONFIG$eval_grid, p = nrow(truth$D_true))
    metric_list <- lapply(seq_len(nrow(truth$beta_curves_eval)), function(j) {
      curve_metrics(curves$beta_hat[j, ], truth$beta_curves_eval[j, ], CONFIG$eval_grid)
    })
    rho_metric <- curve_metrics(curves$rho_hat, truth$rho_curve_eval, CONFIG$eval_grid)

    row <- data.frame(
      stage     = "estimation",
      N         = task$N,
      rho_label = task$rho_label,
      rep_id    = task$rep_id,
      method    = method,
      converged = if (method == "fos") TRUE else isTRUE(fit_res$fit$converged),
      iter      = if (!is.null(fit_res$fit$iter)) fit_res$fit$iter else NA_integer_,
      elapsed_sec      = fit_res$elapsed_sec,
      error_message    = NA_character_,
      stringsAsFactors = FALSE
    )
    for (j in seq_along(metric_list)) {
      row[[paste0("imse_", CONFIG$beta_names[j])]] <- metric_list[[j]][["imse"]]
      row[[paste0("rmse_", CONFIG$beta_names[j])]] <- metric_list[[j]][["rmse"]]
    }
    row[["imse_rho"]] <- rho_metric[["imse"]]
    row[["rmse_rho"]] <- rho_metric[["rmse"]]
    row[["imse_beta_mean"]] <- mean(vapply(metric_list, `[[`, numeric(1), "imse"))
    row[["rmse_beta_mean"]] <- mean(vapply(metric_list, `[[`, numeric(1), "rmse"))
    out_rows[[m_idx]] <- row
  }
  bind_rows(out_rows)
}

lm_worker <- function(task) {
  set.seed(task_seed(task, stage_code = 2L))
  truth <- build_truth(task$rho_label, beta_case = "alt")
  W     <- make_rook_torus_W(task$N)
  Z     <- generate_design_matrix(task$N)
  y     <- y_gen(
    D = truth$D_true, r = truth$r_true, Z = Z, W = W, basis = truth$basis,
    mean = 0, sd = 0.25, L = CONFIG$L_obs
  )

  lm_res <- tryCatch(
    bootstrap_LM_test(
      y = y, Z = Z, W = W, seq_obs = CONFIG$seq_obs,
      basis = truth$basis, M = CONFIG$M_est, B = CONFIG$B_boot, ncores = 1
    ),
    error = function(e) e
  )

  if (inherits(lm_res, "error")) {
    return(data.frame(
      stage = "lm",
      N = task$N,
      rho_label = task$rho_label,
      rep_id = task$rep_id,
      reject = NA_integer_,
      p_value = NA_real_,
      LM = NA_real_,
      B_valid = 0L,
      error_message = conditionMessage(lm_res),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    stage            = "lm",
    N                = task$N,
    rho_label        = task$rho_label,
    rep_id           = task$rep_id,
    reject           = as.integer(lm_res$p_value < CONFIG$alpha),
    p_value          = lm_res$p_value,
    LM               = lm_res$LM,
    B_valid          = lm_res$B_valid,
    stat_name        = lm_res$stat_name,
    error_message    = NA_character_,
    stringsAsFactors = FALSE
  )
}

f_worker <- function(task) {
  set.seed(task_seed(task, stage_code = 3L))
  truth <- build_truth(task$rho_label, beta_case = task$beta_case)
  W     <- make_rook_torus_W(task$N)
  Z     <- generate_design_matrix(task$N)
  y     <- y_gen(
    D = truth$D_true, r = truth$r_true, Z = Z, W = W, basis = truth$basis,
    mean = 0, sd = 0.25, L = CONFIG$L_obs
  )

  fit_res <- tryCatch(
    fsar_func_fast(y = y, Z = Z, W = W, seq = CONFIG$seq_obs,
                   basis = truth$basis, M = CONFIG$M_est, tolerance = CONFIG$tolerance),
    error = function(e) e
  )
  if (inherits(fit_res, "error")) {
    return(data.frame(
      stage     = "f",
      N         = task$N,
      rho_label = task$rho_label,
      beta_case = task$beta_case,
      rep_id    = task$rep_id,
      reject    = NA_integer_,
      p_value   = NA_real_,
      statistic = NA_real_,
      B_valid   = 0L,
      error_message = paste("initial fit:", conditionMessage(fit_res)),
      stringsAsFactors = FALSE
    ))
  }

  D_hat   <- matrix(fit_res$D, nrow = nrow(truth$D_true), ncol = CONFIG$K, byrow = FALSE)
  r_hat   <- matrix(fit_res$r, nrow = 1)
  eta_hat <- matrix(fit_res$eta, nrow = 1)

  f_res <- tryCatch(
    bootstrap_Ftest(
      y = y, Z = Z, W = W, seq_obs = CONFIG$seq_obs,
      basis = truth$basis, M = CONFIG$M_est, B = CONFIG$B_boot, j_test = CONFIG$j_test,
      D_FM = D_hat, r_FM = r_hat, eta_FM = eta_hat,
      ncores = 1, L_boot = CONFIG$L_obs, M_boot = CONFIG$M_boot
    ),
    error = function(e) e
  )

  if (inherits(f_res, "error")) {
    return(data.frame(
      stage            = "f",
      N                = task$N,
      rho_label        = task$rho_label,
      beta_case        = task$beta_case,
      rep_id           = task$rep_id,
      reject           = NA_integer_,
      p_value          = NA_real_,
      statistic        = NA_real_,
      B_valid          = 0L,
      error_message    = conditionMessage(f_res),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    stage              = "f",
    N                  = task$N,
    rho_label          = task$rho_label,
    beta_case          = task$beta_case,
    rep_id             = task$rep_id,
    scenario_role      = if (task$beta_case == "null") "size" else "power",
    reject             = as.integer(f_res$p_value < CONFIG$alpha),
    p_value            = f_res$p_value,
    statistic          = f_res$Fn,
    delta_ISE          = f_res$delta_ISE,
    legacy_Ftype       = f_res$legacy_Ftype,
    B_valid            = f_res$B_valid,
    mean_boot_attempts = mean(f_res$boot_attempts),
    fm_converged_rate  = mean(f_res$boot_fm_converged),
    rm_converged_rate  = mean(f_res$boot_rm_converged),
    error_message      = NA_character_,
    stringsAsFactors   = FALSE
  )
}

ci_worker <- function(task) {
  set.seed(task_seed(task, stage_code = 4L))
  truth <- build_truth(task$rho_label, beta_case = "alt")
  W     <- make_rook_torus_W(task$N)
  Z     <- generate_design_matrix(task$N)
  y     <- y_gen(
    D = truth$D_true, r = truth$r_true, Z = Z, W = W, basis = truth$basis,
    mean = 0, sd = 0.25, L = CONFIG$L_obs
  )

  fit_res <- tryCatch(
    fsar_func_fast(y = y, Z = Z, W = W, seq = CONFIG$seq_obs,
                   basis = truth$basis, M = CONFIG$M_est, tolerance = CONFIG$tolerance),
    error = function(e) e
  )
  if (inherits(fit_res, "error")) {
    return(data.frame(
      stage = "ci",
      N = task$N,
      rho_label = task$rho_label,
      rep_id = task$rep_id,
      boot_converged_rate = NA_real_,
      error_message = paste("initial fit:", conditionMessage(fit_res)),
      stringsAsFactors = FALSE
    ))
  }

  D_hat <- matrix(fit_res$D, nrow = nrow(truth$D_true), ncol = CONFIG$K, byrow = FALSE)
  r_hat <- matrix(fit_res$r, nrow = 1)
  eta_hat <- matrix(fit_res$eta, nrow = 1)

  ci_res <- tryCatch(
    bootstrap_CI(
      y = y, Z = Z, W = W, seq_obs = CONFIG$seq_obs,
      basis     = truth$basis, M = CONFIG$M_est, B = CONFIG$B_boot, alpha = CONFIG$alpha,
      D_hat     = D_hat, r_hat = r_hat, eta_hat = eta_hat, ncores = 1,
      L_boot    = CONFIG$L_obs, M_boot = CONFIG$M_boot,
      ci_method = CONFIG$ci_method, multiplier = CONFIG$ci_multiplier,
      rho_clip  = CONFIG$rho_clip
    ),
    error = function(e) e
  )

  if (inherits(ci_res, "error")) {
    return(data.frame(
      stage               = "ci",
      N                   = task$N,
      rho_label           = task$rho_label,
      rep_id              = task$rep_id,
      boot_converged_rate = NA_real_,
      error_message       = conditionMessage(ci_res),
      stringsAsFactors    = FALSE
    ))
  }

  eval_grid       <- ci_res$eval_grids
  Phi_eval        <- eval.basis(eval_grid, truth$basis)
  beta_truth_eval <- truth$beta_curves_eval
  if (length(eval_grid) != ncol(beta_truth_eval)) {
    beta_truth_eval <- make_beta_curves(eval_grid, beta_case = "alt")
  }
  rho_truth_eval <- if (length(eval_grid) == length(truth$rho_curve_eval)) {
    truth$rho_curve_eval
  } else {
    make_rho_curve(eval_grid, task$rho_label)
  }

  row <- data.frame(
    stage               = "ci",
    N                   = task$N,
    rho_label           = task$rho_label,
    rep_id              = task$rep_id,
    ci_method           = ci_res$ci_method,
    multiplier          = ci_res$multiplier,
    boot_converged_rate = mean(ci_res$boot_converged),
    median_boot_iter    = median(ci_res$boot_iter),
    max_boot_iter       = max(ci_res$boot_iter),
    error_message       = NA_character_,
    stringsAsFactors    = FALSE
  )

  for (j in seq_along(CONFIG$beta_names)) {
    lower  <- ci_res$beta_ci[[j]]$lower
    upper   <- ci_res$beta_ci[[j]]$upper
    truth_j <- beta_truth_eval[j, ]
    row[[paste0("cov_", CONFIG$beta_names[j])]] <- mean(truth_j >= lower & truth_j <= upper)
    row[[paste0("len_", CONFIG$beta_names[j])]] <- mean(upper - lower)
  }
  row[["cov_rho"]]  <- mean(rho_truth_eval >= ci_res$rho_ci$lower & rho_truth_eval <= ci_res$rho_ci$upper)
  row[["len_rho"]]  <- mean(ci_res$rho_ci$upper - ci_res$rho_ci$lower)
  row[["mean_cov"]] <- mean(unlist(row[paste0("cov_", c(CONFIG$beta_names, "rho"))]), na.rm = TRUE)
  row[["mean_len"]] <- mean(unlist(row[paste0("len_", c(CONFIG$beta_names, "rho"))]), na.rm = TRUE)
  row
}

parallel_task_apply <- function(tasks, worker_fun, seed_offset = 0L) {
  if (length(tasks) == 0L) {
    return(list())
  }
  if (CONFIG$ncores <= 1L) {
    return(lapply(tasks, worker_fun))
  }

  cl <- parallel::makeCluster(CONFIG$ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterSetRNGStream(cl, iseed = CONFIG$seed + seed_offset)
  parallel::clusterExport(
    cl,
    varlist = c(
      "CONFIG", "SCRIPT_DIR", "ROOT_DIR", "W_CACHE",
      "make_basis", "make_rook_torus_W", "scale_curve", "make_beta_curves",
      "make_rho_curve", "curve_matrix_to_coef", "curve_to_coef_row", "build_truth",
      "generate_design_matrix", "trapz_integral", "curve_metrics", "task_seed",
      "extract_fit_curves", "fit_method",
      "estimation_worker", "lm_worker", "f_worker", "ci_worker"
    ),
    envir = .GlobalEnv
  )
  parallel::clusterEvalQ(cl, {
    library(dplyr)
    library(fda)
    library(Matrix)
    library(parallel)
    setwd(SCRIPT_DIR)
    options(sfsar_root_dir = ROOT_DIR)
    source(file.path(ROOT_DIR, "fos.R"))
    source(file.path(ROOT_DIR, "fsar_func.R"))
    source(file.path(ROOT_DIR, "fsar_scalar.R"))
    source(file.path(ROOT_DIR, "y_gen.R"))
    source(file.path(ROOT_DIR, "bootstrap_inference.R"))
    NULL
  })
  parallel::parLapplyLB(cl, tasks, worker_fun)
}

update_stage_outputs <- function(stage_name, scenario_key, scenario_df, summary_group) {
  scenario_rds       <- file.path(CHK_DIR, paste0(stage_name, "_", scenario_key, ".rds"))
  saveRDS(scenario_df, scenario_rds)

  pattern     <- paste0("^", stage_name, "_.*\\.rds$")
  stage_files <- list.files(CHK_DIR, pattern = pattern, full.names = TRUE)
  stage_raw   <- bind_rows(lapply(stage_files, readRDS))
  stage_raw   <- stage_raw %>% arrange(across(any_of(names(stage_raw))))

  raw_csv     <- file.path(CSV_DIR, paste0(stage_name, "_raw.csv"))
  raw_rds     <- file.path(RDS_DIR, paste0(stage_name, "_raw.rds"))
  write.csv(stage_raw, raw_csv, row.names = FALSE)
  saveRDS(stage_raw, raw_rds)

  numeric_cols <- names(stage_raw)[vapply(stage_raw, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("rep_id", summary_group))
  if (length(numeric_cols) > 0L) {
    stage_summary <- stage_raw %>%
      group_by(across(all_of(summary_group))) %>%
      summarise(
        across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE), .names = "{.col}"),
        n_records = n(),
        .groups = "drop"
      )
    summary_csv <- file.path(CSV_DIR, paste0(stage_name, "_summary.csv"))
    summary_rds <- file.path(RDS_DIR, paste0(stage_name, "_summary.rds"))
    write.csv(stage_summary, summary_csv, row.names = FALSE)
    saveRDS(stage_summary, summary_rds)
  }
}

run_estimation_stage <- function() {
  for (N in CONFIG$N_list) {
    for (rho_label in CONFIG$rho_labels) {
      scenario_key <- paste0("N", N, "_", rho_label)
      cat("\n[estimation]", scenario_key, "\n")
      tasks <- lapply(seq_len(CONFIG$n_sim), function(rep_id) {
        list(N = N, rho_label = rho_label, rep_id = rep_id)
      })
      res_list <- parallel_task_apply(tasks, estimation_worker, seed_offset = 11L + N)
      scenario_df <- bind_rows(res_list)
      update_stage_outputs("estimation", scenario_key, scenario_df,
                           summary_group = c("N", "rho_label", "method"))
    }
  }
}

run_lm_stage <- function() {
  for (N in CONFIG$N_list) {
    for (rho_label in CONFIG$rho_labels) {
      scenario_key <- paste0("N", N, "_", rho_label)
      cat("\n[lm]", scenario_key, "\n")
      tasks <- lapply(seq_len(CONFIG$n_sim), function(rep_id) {
        list(N = N, rho_label = rho_label, rep_id = rep_id)
      })
      res_list <- parallel_task_apply(tasks, lm_worker, seed_offset = 21L + N)
      scenario_df <- bind_rows(res_list)
      update_stage_outputs("lm", scenario_key, scenario_df,
                           summary_group = c("N", "rho_label", "stat_name"))
    }
  }
}

run_f_stage <- function() {
  for (N in CONFIG$N_list) {
    for (rho_label in CONFIG$rho_labels) {
      for (beta_case in CONFIG$f_beta_cases) {
        scenario_key <- paste0("N", N, "_", rho_label, "_", beta_case)
        cat("\n[f]", scenario_key, "\n")
        tasks <- lapply(seq_len(CONFIG$n_sim), function(rep_id) {
          list(N = N, rho_label = rho_label, beta_case = beta_case, rep_id = rep_id)
        })
        res_list <- parallel_task_apply(tasks, f_worker, seed_offset = 31L + N)
        scenario_df <- bind_rows(res_list)
        update_stage_outputs("f", scenario_key, scenario_df,
                             summary_group = c("N", "rho_label", "beta_case", "scenario_role"))
      }
    }
  }
}

run_ci_stage <- function() {
  for (N in CONFIG$N_list) {
    for (rho_label in CONFIG$rho_labels) {
      scenario_key <- paste0("N", N, "_", rho_label)
      cat("\n[ci]", scenario_key, "\n")
      tasks <- lapply(seq_len(CONFIG$n_sim), function(rep_id) {
        list(N = N, rho_label = rho_label, rep_id = rep_id)
      })
      res_list <- parallel_task_apply(tasks, ci_worker, seed_offset = 41L + N)
      scenario_df <- bind_rows(res_list)
      update_stage_outputs("ci", scenario_key, scenario_df,
                           summary_group = c("N", "rho_label", "ci_method", "multiplier"))
    }
  }
}

saveRDS(CONFIG, file.path(RDS_DIR, "simulation_config.rds"))
write.csv(
  data.frame(
    name = names(CONFIG),
    value = vapply(CONFIG, function(x) {
      if (length(x) == 0L) {
        return("")
      }
      paste(as.character(x), collapse = ",")
    }, character(1)),
    stringsAsFactors = FALSE
  ),
  file.path(CSV_DIR, "simulation_config.csv"),
  row.names = FALSE
)

cat("Simulation configuration:\n")
print(CONFIG)

if ("estimation" %in% CONFIG$stages) {
  run_estimation_stage()
}
if ("lm" %in% CONFIG$stages) {
  run_lm_stage()
}
if ("f" %in% CONFIG$stages) {
  run_f_stage()
}
if ("ci" %in% CONFIG$stages) {
  run_ci_stage()
}

