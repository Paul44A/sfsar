
################################################################################
## Bootstrap inference for SFSAR models
## Implements Algorithm 1 (LM test), Algorithm 2 (F-type test), and Algorithm 3 (Pointwise CI) 
################################################################################

library(fda)
library(parallel)
if (requireNamespace("Matrix", quietly = TRUE)) {
  W_can_be_sparse <- TRUE
} else {
  W_can_be_sparse <- FALSE
}

get_current_script_path <- function(default_name = "bootstrap_inference.R") {
  opt_root <- getOption("sfsar_root_dir", default = NULL)
  if (!is.null(opt_root) && dir.exists(opt_root)) {
    candidate <- file.path(opt_root, default_name)
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- normalizePath(
      sub("^--file=", "", file_arg[1]),
      winslash = "/",
      mustWork = FALSE
    )
    if (basename(script_path) == default_name) {
      return(script_path)
    }
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile)) {
      script_path <- normalizePath(ofile, winslash = "/", mustWork = FALSE)
      if (basename(script_path) == default_name) {
        return(script_path)
      }
    }
  }

  normalizePath(file.path(getwd(), default_name), winslash = "/", mustWork = FALSE)
}

BOOTSTRAP_INFERENCE_FILE <- get_current_script_path()

get_rho_stability_bound <- function(W, margin = 0.99) {
  row_abs_sum_max <- max(rowSums(abs(W)))
  if (!is.finite(row_abs_sum_max) || row_abs_sum_max <= 0) {
    return(Inf)
  }
  margin / row_abs_sum_max
}

.fsar_w_cache <- new.env(parent = emptyenv())

get_fsar_w_cache <- function(W, margin = 0.99) {
  if (exists("W_ref", envir = .fsar_w_cache, inherits = FALSE) &&
      identical(W, .fsar_w_cache$W_ref)) {
    return(.fsar_w_cache$data)
  }

  eigvals         <- eigen(W, only.values = TRUE)$values
  spectral_radius <- max(Mod(eigvals))
  rho_row_bound   <- get_rho_stability_bound(W, margin = margin)
  rho_spec_bound  <- if (is.finite(spectral_radius) && spectral_radius > 0) {
    margin / spectral_radius
  } else {
    0
  }
  rho_bound <- if (is.finite(rho_row_bound) && rho_row_bound > 0) {
    min(rho_row_bound, rho_spec_bound)
  } else {
    rho_spec_bound
  }

  data <- list(
    Ws        = W + t(W),
    WtW       = crossprod(W),
    eigvals   = eigvals,
    rho_bound = rho_bound
  )
  .fsar_w_cache$W_ref <- W
  .fsar_w_cache$data  <- data
  data
}

project_vals_to_basis_coef <- function(vals, Phi) {
  as.numeric(qr.solve(Phi, vals))
}

normalize_optional_coef <- function(x, expected_len = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  x_num <- as.numeric(x)
  if (length(x_num) == 0) {
    return(NULL)
  }
  if (!is.null(expected_len) && length(x_num) != expected_len) {
    stop("Coefficient length mismatch: expected ", expected_len,
         ", got ", length(x_num), ".")
  }
  x_num
}

normalize_constraint_coef <- function(constraint_c, q, K) {
  if (is.null(constraint_c)) {
    return(matrix(0, nrow = q, ncol = K))
  }
  if (inherits(constraint_c, "fd")) {
    coef_mat <- t(constraint_c$coefs)
  } else if (is.matrix(constraint_c)) {
    coef_mat <- if (nrow(constraint_c) == q && ncol(constraint_c) == K) {
      constraint_c
    } else if (nrow(constraint_c) == K && ncol(constraint_c) == q) {
      t(constraint_c)
    } else {
      stop("constraint_c matrix must have dimension q×K or K×q.")
    }
  } else {
    coef_num <- as.numeric(constraint_c)
    if (length(coef_num) != q * K) {
      stop("constraint_c must contain q*K coefficients.")
    }
    coef_mat <- matrix(coef_num, nrow = q, ncol = K, byrow = FALSE)
  }
  if (anyNA(coef_mat) || any(!is.finite(coef_mat))) {
    stop("constraint_c must be finite.")
  }
  coef_mat
}

prepare_linear_constraint_info <- function(p, K, j_test = NULL,
                                           constraint_C = NULL, constraint_c = NULL) {
  if (is.null(constraint_C)) {
    j_test <- sort(unique(as.integer(j_test)))
    if (length(j_test) == 0L || anyNA(j_test)) {
      stop("bootstrap_Ftest requires either j_test or constraint_C.")
    }
    if (any(j_test < 1L | j_test > p)) {
      stop("j_test contains indices outside 1:p.")
    }
    C_mat <- diag(p)[j_test, , drop = FALSE]
  } else {
    C_mat <- as.matrix(constraint_C)
    if (ncol(C_mat) != p) {
      stop("constraint_C must have p columns.")
    }
  }

  q <- nrow(C_mat)
  if (q <= 0L) {
    stop("constraint_C must have at least one row.")
  }
  if (qr(C_mat)$rank != q) {
    stop("constraint_C must have full row rank.")
  }

  c_coef       <- normalize_constraint_coef(constraint_c, q = q, K = K)
  D_particular <- matrix(0, nrow = p, ncol = K)
  for (k in seq_len(K)) {
    D_particular[, k] <- qr.solve(C_mat, c_coef[, k])
  }

  qr_Ct      <- qr(t(C_mat))
  Q_full     <- qr.Q(qr_Ct, complete = TRUE)
  null_basis <- Q_full[, seq.int(q + 1L, p), drop = FALSE]
  p_reduced  <- ncol(null_basis)
  if (p_reduced <= 0L) {
    stop("Constraints leave no free regression coefficients; this case is not supported.")
  }

  list(
    C            = C_mat,
    c_coef       = c_coef,
    q            = q,
    D_particular = D_particular,
    null_basis   = null_basis,
    p_reduced    = p_reduced
  )
}

estimate_restricted_sfsar <- function(y, Z, W, seq_obs, basis, M, tolerance,
                                      constraint_info,
                                      D_init_reduced = NULL,
                                      r_init = NULL, eta_init = NULL) {
  Phi_obs     <- t(eval.basis(seq_obs, basis))
  offset_eval <- Z %*% constraint_info$D_particular %*% Phi_obs
  y_tilde     <- y - offset_eval
  Z_reduced   <- Z %*% constraint_info$null_basis

  fit <- fsar_func_fast(
    y = y_tilde, Z = Z_reduced, W = W, seq = seq_obs,
    basis = basis, M = M, tolerance = tolerance,
    D_init = D_init_reduced, r_init = r_init, eta_init = eta_init
  )

  D_reduced <- matrix(fit$D, nrow = constraint_info$p_reduced, ncol = basis$nbasis,
                      byrow = FALSE)
  D_full <- constraint_info$D_particular + constraint_info$null_basis %*% D_reduced

  list(
    D_full    = D_full,
    D_reduced = D_reduced,
    r         = matrix(fit$r, nrow = 1),
    eta       = matrix(fit$eta, nrow = 1),
    converged = isTRUE(fit$converged),
    iter      = fit$iter,
    beta_dif  = fit$beta_dif,
    rho_dif   = fit$rho_dif
  )
}

compute_projection_info <- function(Z, include_residualizer = TRUE) {
  qr_Z     <- qr(Z)
  Q        <- qr.Q(qr_Z, complete = FALSE)
  hat_diag <- rowSums(Q^2)
  M <- NULL
  if (include_residualizer) {
    M <- diag(nrow(Z)) - tcrossprod(Q)
  }
  list(qr = qr_Z, Q = Q, hat_diag = hat_diag, M = M)
}

compute_leverage_row_scale <- function(hat_diag, hc_type = c("none", "hc2", "hc3"),
                                       floor = 1e-8) {
  hc_type <- match.arg(hc_type)
  base    <- pmax(1 - as.numeric(hat_diag), floor)
  switch(
    hc_type,
    none = rep(1, length(base)),
    hc2  = sqrt(base),
    hc3  = base
  )
}

run_bootstrap_jobs <- function(B, worker_fun, ncores = 1, progress_label = "bootstrap") {
  idx <- seq_len(B)
  if (ncores > 1L) {
    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      wd          <- getwd()
      script_file <- BOOTSTRAP_INFERENCE_FILE
      parallel::clusterCall(cl, function(path) setwd(path), wd)
      parallel::clusterExport(cl, varlist = "script_file", envir = environment())
      parallel::clusterEvalQ(cl, {
        library(fda)
        library(parallel)
        source(script_file)
        NULL
      })
      return(parallel::parLapplyLB(cl, idx, worker_fun))
    }
    return(parallel::mclapply(idx, worker_fun, mc.cores = ncores))
  }

  out <- vector("list", B)
  for (b in idx) {
    cat("\r ", progress_label, " ", b, "/", B, sep = "")
    out[[b]] <- worker_fun(b)
  }
  cat("\n")
  out
}

eval_bounded_rho_from_eta <- function(Phi, eta_coef, rho_bound) {
  eta_coef <- as.numeric(eta_coef)
  rho_bound * tanh(as.numeric(Phi %*% eta_coef))
}

eta_from_rho_coef <- function(rho_coef, Phi, rho_bound, eps = 1e-6) {
  rho_vals <- as.numeric(Phi %*% as.numeric(rho_coef))
  rho_safe <- pmin(pmax(rho_vals, -(rho_bound - eps)), rho_bound - eps)
  eta_vals <- atanh(rho_safe / rho_bound)
  matrix(project_vals_to_basis_coef(eta_vals, Phi), nrow = 1)
}

eval_rho_curve <- function(grids, basis, r = NULL, eta = NULL, rho_bound = NULL) {
  Phi <- eval.basis(grids, basis)
  if (!is.null(eta)) {
    if (is.null(rho_bound)) {
      stop("eval_rho_curve requires rho_bound when eta is supplied.")
    }
    return(eval_bounded_rho_from_eta(Phi, eta, rho_bound))
  }
  as.numeric(Phi %*% as.numeric(r))
}

optimize_bounded_rho_step <- function(Phi, itv, q0, q1, q2, rho_bound, eigvals,
                                      N, eta_init, maxit = 100) {
  eta_init  <- as.numeric(eta_init)
  sse_floor <- max(.Machine$double.eps,
                   1e-10 * mean(q0[q0 > 0], na.rm = TRUE))
  if (!is.finite(sse_floor)) {
    sse_floor <- 1e-8
  }

  eig_mat   <- matrix(eigvals, nrow = nrow(Phi), ncol = length(eigvals),byrow = TRUE)
  last_par  <- NULL
  last_eval <- NULL

  evaluate_once <- function(par) {
    if (!is.null(last_par) && max(abs(par - last_par)) < 1e-12) {
      return(last_eval)
    }

    eta_vals     <- as.numeric(Phi %*% par)
    tanh_vals    <- tanh(eta_vals)
    rho_vals     <- rho_bound * tanh_vals
    drho_deta    <- rho_bound * (1 - tanh_vals^2)
    sse_vals     <- q0 - rho_vals * q1 + rho_vals^2 * q2
    sse_safe     <- pmax(sse_vals, sse_floor)
    denom        <- 1 - outer(rho_vals, eigvals, `*`)
    denom_mod    <- pmax(Mod(denom), sse_floor)
    logdet_vals  <- rowSums(log(denom_mod))
    logdet_score <- rowSums(Re(eig_mat / denom))

    obj_score <- (-q1 + 2 * rho_vals * q2) / sse_safe +
      (2 / N) * logdet_score

    last_par <<- par
    last_eval <<- list(
      objective = sum(log(sse_safe / N) - (2 / N) * logdet_vals) * itv,
      gradient  = as.numeric(crossprod(Phi, obj_score * drho_deta) * itv),
      rho_vals  = rho_vals
    )
    last_eval
  }

  objective <- function(par) {
    evaluate_once(par)$objective
  }
  gradient <- function(par) {
    evaluate_once(par)$gradient
  }

  opt <- optim(
    par    = eta_init,
    fn     = objective,
    gr     = gradient,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-8)
  )

  opt_eval <- evaluate_once(opt$par)
  list(
    eta_coef = matrix(opt$par, nrow = 1),
    rho_vals = opt_eval$rho_vals,
    rho_coef = matrix(project_vals_to_basis_coef(opt_eval$rho_vals, Phi), nrow = 1),
    objective = opt$value,
    convergence = opt$convergence
  )
}


fsar_func_fast <- function(y, Z, W, seq, basis, M, tolerance,
                           D_init = NULL, r_init = NULL, eta_init = NULL,
                           max_iter = 1000) {
  N       <- dim(W)[1];  p <- dim(Z)[2];  K <- basis$nbasis
  C       <- t(smooth.basis(seq, t(y), basis)$fd$coef)  # N×K
  W_cache <- get_fsar_w_cache(W)
  Ws      <- W_cache$Ws
  WtW     <- W_cache$WtW
  eigvals <- W_cache$eigvals
  ZtZ     <- crossprod(Z)              
  ZtWsZ   <- crossprod(Z, Ws  %*% Z)   
  ZtWtWZ  <- crossprod(Z, WtW %*% Z)   
  ZtC     <- crossprod(Z, C)            
  ZtWsC   <- crossprod(Z, Ws  %*% C) 
  ZtWtWC  <- crossprod(Z, WtW %*% C)   

  grids        <- seq(basis$rangeval[1], basis$rangeval[2], length.out = M)
  itv          <- grids[2] - grids[1]
  Phi          <- eval.basis(grids, basis)    # M × K
  rho_bound    <- W_cache$rho_bound
  D_init_num   <- normalize_optional_coef(D_init, expected_len = K * p)
  r_init_num   <- normalize_optional_coef(r_init, expected_len = K)
  eta_init_num <- normalize_optional_coef(eta_init, expected_len = K)

  D <- if (!is.null(D_init_num)) D_init_num else rep(0, K * p)
  r <- if (!is.null(r_init_num)) {
    matrix(r_init_num, 1, K)
  } else if (!is.null(eta_init_num)) {
    matrix(
      project_vals_to_basis_coef(
        eval_bounded_rho_from_eta(Phi, eta_init_num, rho_bound),
        Phi
      ),
      nrow = 1
    )
  } else {
    matrix(0, 1, K)
  }
  eta <- if (!is.null(eta_init_num) && rho_bound > 0) {
    matrix(eta_init_num, nrow = 1)
  } else if (!is.null(r_init_num) && rho_bound > 0) {
    eta_from_rho_coef(r, Phi, rho_bound)
  } else {
    matrix(0, nrow = 1, ncol = K)
  }

  converged <- FALSE
  beta_dif  <- Inf
  rho_dif   <- Inf
  for (iter in 1:max_iter) {
    rho_v  <- if (rho_bound > 0) {
      eval_bounded_rho_from_eta(Phi, eta, rho_bound)
    } else {
      as.numeric(Phi %*% t(r))
    }
    G0 <- crossprod(Phi) * itv                    
    G1 <- crossprod(Phi, Phi * rho_v) * itv       
    G2 <- crossprod(Phi, Phi * rho_v^2) * itv    

    BkA   <- kronecker(G0, ZtZ) - kronecker(G1, ZtWsZ) + kronecker(G2, ZtWtWZ)
    E_sum <- ZtC %*% G0 - ZtWsC %*% G1 + ZtWtWC %*% G2   # p×K
    D_new <- as.numeric(solve(BkA, as.numeric(E_sum)))

    D_mat <- matrix(D_new, nrow = p, ncol = K, byrow = FALSE)
    Res   <- C - Z %*% D_mat                      
    Q0    <- crossprod(Res)                        
    WRes  <- W %*% Res                            
    Q1    <- crossprod(WRes)                       
    Q2    <- crossprod(Res, Ws %*% Res)           

    q0 <- rowSums(Phi * (Phi %*% Q0))             
    q1 <- rowSums(Phi * (Phi %*% Q2))             
    q2 <- rowSums(Phi * (Phi %*% Q1))              

    if (rho_bound > 0) {
      rho_opt <- optimize_bounded_rho_step(
        Phi = Phi, itv = itv,
        q0 = q0, q1 = q1, q2 = q2,
        rho_bound = rho_bound, eigvals = eigvals, N = N,
        eta_init = eta, maxit = 100
      )
      eta_new      <- rho_opt$eta_coef
      r_new        <- rho_opt$rho_coef
      rho_new_vals <- rho_opt$rho_vals
    } else {
      eta_new      <- matrix(0, nrow = 1, ncol = K)
      r_new        <- matrix(0, nrow = 1, ncol = K)
      rho_new_vals <- rep(0, length(grids))
    }

    # Check convergence
    beta_dif <- mean(abs(eval.fd(grids,
      fd(t(matrix(D_new - D, p, K, byrow = FALSE)), basis)))) / p
    rho_dif  <- mean(abs(rho_new_vals - rho_v))

    D   <- D_new
    r   <- r_new
    eta <- eta_new
    if (beta_dif < tolerance && rho_dif < tolerance) {
      converged <- TRUE; break
    }
  }
  return(list(
    D         = D,
    r         = r,
    iter      = iter,
    converged = converged,
    beta_dif  = beta_dif,
    rho_dif   = rho_dif,
    eta       = eta,
    rho_bound = rho_bound
  ))
}

################################################################################
## 1. LM Test for Spatial Correlation 
################################################################################

compute_LM_components <- function(C, Z, W, D, basis, projection_info = NULL,
                                  G = NULL, hc_type = c("none", "hc2", "hc3"),
                                  project_w = TRUE, energy_eval = NULL) {
  # C:     N×K coefficient matrix of response basis expansion
  # Z:     N×p design matrix
  # W:     N×N spatial weight matrix
  # D:     p×K coefficient matrix of beta(t)
  # basis: basis object
  N <- nrow(C)
  if (is.null(G)) {
    G <- inprod(basis, basis)
  }

  E             <- C - Z %*% D
  numerator     <- sum((t(E) %*% W %*% E) * G)
  residual_norm <- sum(crossprod(E) * G)
  info_scale    <- sum(t(W) * W) / N^2
  lm_scale      <- max(residual_norm * info_scale, .Machine$double.eps)

  list(
    LM            = numerator / lm_scale,
    numerator     = numerator,
    residual_norm = residual_norm,
    info_scale    = info_scale,
    denominator   = lm_scale
  )
}

compute_LM_stat <- function(C, Z, W, D, basis, projection_info = NULL,
                            G = NULL, hc_type = c("none", "hc2", "hc3"),
                            project_w = TRUE, energy_eval = NULL) {
  compute_LM_components(
    C = C, Z = Z, W = W, D = D, basis = basis,
    projection_info = projection_info, G = G,
    hc_type = hc_type, project_w = project_w,
    energy_eval = energy_eval
  )$LM
}

estimate_fos_coef_qr <- function(Z, C, projection_info = NULL) {
  D_hat <- if (is.null(projection_info)) {
    qr.solve(Z, C)
  } else {
    qr.coef(projection_info$qr, C)
  }
  if (anyNA(D_hat) || any(!is.finite(D_hat))) {
    stop("FoS QR projection failed to produce finite coefficients.")
  }
  D_hat
}

compute_fos_bootstrap_residuals <- function(y, Z, D, basis, grids, projection_info = NULL,
                                            hc_type = c("none", "hc2", "hc3")) {
  Phi         <- t(eval.basis(grids, basis))             
  fitted_eval <- Z %*% D %*% Phi                 
  V_hat       <- y - fitted_eval                       
  V_star      <- sweep(V_hat, 2, colMeans(V_hat))      
  list(
    V_star      = V_star,
    fitted_eval = fitted_eval,
    grids       = grids,
    Phi         = Phi
  )
}

bootstrap_LM_test <- function(y, Z, W, seq_obs, basis, M, B, ncores = 1,
                              multiplier = c("mammen", "rademacher", "gaussian"),
                              hc_type = c("none", "hc2", "hc3"),
                              project_w = FALSE) {
  # y: N×L response matrix (stations × time points)
  # Z: N×p design matrix
  # W: N×N spatial weight matrix
  # Returns: LM stat, bootstrap p-value, bootstrap distribution
  multiplier <- match.arg(multiplier)
  hc_type <- match.arg(hc_type)
  if (hc_type != "none") {
    warning("Algorithm 1 uses raw centered residuals; hc_type is ignored.")
  }
  if (isTRUE(project_w)) {
    warning("Algorithm 1 uses W directly; project_w is ignored.")
  }
  projection_info <- compute_projection_info(Z, include_residualizer = FALSE)
  G               <- inprod(basis, basis)

  C        <- t(smooth.basis(seq_obs, t(y), basis)$fd$coef)     # N×K
  D_fos    <- estimate_fos_coef_qr(Z, C, projection_info = projection_info)
  res_info <- compute_fos_bootstrap_residuals(
    y, Z, D_fos, basis, seq_obs
  )

  LM_obs <- compute_LM_stat(
    C, Z, W, D_fos, basis, G = G
  )

  # Step 3-5: Bootstrap
  one_LM_boot <- function(b) {
    eta    <- draw_bootstrap_multipliers(nrow(y), method = multiplier)
    y_boot <- res_info$fitted_eval + res_info$V_star * eta
    tryCatch({
      C_boot <- t(smooth.basis(seq_obs, t(y_boot), basis)$fd$coef)
      D_boot <- estimate_fos_coef_qr(Z, C_boot, projection_info = projection_info)
      stat_b <- compute_LM_stat(
        C_boot, Z, W, D_boot, basis, G = G
      )
      list(stat = stat_b, valid = is.finite(stat_b), error = NULL)
    }, error = function(e) {
      list(stat = NA_real_, valid = FALSE, error = conditionMessage(e))
    })
  }
  LM_boot_list <- run_bootstrap_jobs(B, one_LM_boot, ncores = ncores, progress_label = "LM bootstrap")
  LM_boot      <- vapply(LM_boot_list, function(x) x$stat, numeric(1))
  boot_valid   <- vapply(LM_boot_list, function(x) isTRUE(x$valid), logical(1))
  n_valid      <- sum(boot_valid)
  if (n_valid == 0L) {
    stop("LM bootstrap produced no valid replicates.")
  }
  p_value <- mean(LM_boot[boot_valid] >= LM_obs)
  return(list(
    LM         = LM_obs,
    p_value    = p_value,
    LM_boot    = LM_boot,
    boot_valid = boot_valid,
    B_valid    = n_valid,
    multiplier = multiplier,
    stat_name  = "LM_eq13"
  ))
}

################################################################################
## 2. ISE computation 
################################################################################

compute_ISE <- function(E, W, r, basis, M, eta = NULL, rho_bound = NULL) {
  A0 <- crossprod(E)
  A1 <- t(E) %*% (W + t(W)) %*% E
  A2 <- t(E) %*% crossprod(W) %*% E

  grids <- seq(basis$rangeval[1], basis$rangeval[2], length.out = M)
  itv   <- grids[2] - grids[1]
  Phi   <- eval.basis(grids, basis)              # M × K
  rho_vals <- if (!is.null(eta)) {
    if (is.null(rho_bound)) {
      rho_bound <- get_fsar_w_cache(W)$rho_bound
    }
    eval_bounded_rho_from_eta(Phi, eta, rho_bound)
  } else {
    as.numeric(r %*% t(Phi))
  }
  PhiA0 <- Phi %*% A0
  PhiA1 <- Phi %*% A1
  PhiA2 <- Phi %*% A2
  a0    <- rowSums(Phi * PhiA0)
  a1    <- rowSums(Phi * PhiA1)
  a2    <- rowSums(Phi * PhiA2)
  ISE   <- sum(a0 - rho_vals * a1 + rho_vals^2 * a2) * itv
  return(ISE)
}

################################################################################
## 3. Bootstrap helpers: residuals & response generation
################################################################################

compute_residuals_on_grid <- function(D, r, Z, W, C, basis, grids,
                                      eta = NULL, rho_bound = NULL) {
  N      <- nrow(C)
  L_grid <- length(grids)
  Phi    <- t(eval.basis(grids, basis))           

  E       <- C - Z %*% D                        
  WE      <- W %*% E                           
  E_eval  <- E  %*% Phi                          
  WE_eval <- WE %*% Phi                       
  rho_vals <- if (!is.null(eta)) {
    if (is.null(rho_bound)) {
      rho_bound <- get_fsar_w_cache(W)$rho_bound
    }
    eval_bounded_rho_from_eta(t(Phi), eta, rho_bound)
  } else {
    as.numeric(r %*% Phi)
  }

  V_eval <- E_eval - WE_eval * matrix(rho_vals, nrow = N, ncol = L_grid, byrow = TRUE)
  V_star <- sweep(V_eval, 2, colMeans(V_eval))

  return(list(V_star = V_star, rho_vals = rho_vals, grids = grids, Phi = Phi))
}

draw_bootstrap_multipliers <- function(N, method = c("gaussian", "rademacher", "mammen")) {
  method <- match.arg(method)
  if (method == "gaussian") {
    return(rnorm(N))
  }
  if (method == "mammen") {
    sqrt5   <- sqrt(5)
    val_lo  <- (1 - sqrt5) / 2
    val_hi  <- (1 + sqrt5) / 2
    prob_hi <- (sqrt5 - 1) / (2 * sqrt5)
    return(ifelse(runif(N) < prob_hi, val_hi, val_lo))
  }
  sample(c(-1, 1), N, replace = TRUE)
}

stabilize_rho_path <- function(rho_vals, rho_clip = 0.99) {
  if (is.null(rho_clip) || !is.finite(rho_clip)) {
    return(as.numeric(rho_vals))
  }
  rho_vals <- as.numeric(rho_vals)
  max_abs  <- max(abs(rho_vals))
  if (!is.finite(max_abs) || max_abs <= rho_clip) {
    return(rho_vals)
  }
  rho_vals * (rho_clip / max_abs)
}

compute_pointwise_bootstrap_ci <- function(est, boot_mat, alpha = 0.05,
                                           ci_method = c("wald", "normal", "bc", "percentile", "basic")) {
  ci_method   <- match.arg(ci_method)
  B           <- nrow(boot_mat)
  G           <- ncol(boot_mat)
  eps         <- max(1 / (B + 1), 1e-6)
  clamp_probs <- function(x) pmin(pmax(x, eps), 1 - eps)
  q_low       <- alpha / 2
  q_upp       <- 1 - alpha / 2
  perc_low    <- apply(boot_mat, 2, quantile, probs = q_low, names = FALSE, type = 8)
  perc_upp    <- apply(boot_mat, 2, quantile, probs = q_upp, names = FALSE, type = 8)
  basic_low   <- 2 * est - perc_upp
  basic_upp   <- 2 * est - perc_low

  boot_mean   <- colMeans(boot_mat)
  boot_sd     <- apply(boot_mat, 2, sd)
  boot_sd[!is.finite(boot_sd)] <- 0
  zcrit      <- qnorm(q_upp)
  wald_low   <- est - zcrit * boot_sd
  wald_upp   <- est + zcrit * boot_sd
  normal_ctr <- 2 * est - boot_mean
  normal_low <- normal_ctr - zcrit * boot_sd
  normal_upp <- normal_ctr + zcrit * boot_sd

  est_mat      <- matrix(est, nrow = B, ncol = G, byrow = TRUE)
  prop_less    <- (colSums(boot_mat < est_mat) + 0.5 * colSums(boot_mat == est_mat)) / B
  z0           <- qnorm(clamp_probs(prop_less))
  z_alpha      <- qnorm(c(q_low, q_upp))
  bc_probs_low <- clamp_probs(pnorm(2 * z0 + z_alpha[1]))
  bc_probs_upp <- clamp_probs(pnorm(2 * z0 + z_alpha[2]))
  bc_low       <- vapply(seq_len(G), function(g) {
    quantile(boot_mat[, g], probs = bc_probs_low[g], names = FALSE, type = 8)
  }, numeric(1))
  bc_upp <- vapply(seq_len(G), function(g) {
    quantile(boot_mat[, g], probs = bc_probs_upp[g], names = FALSE, type = 8)
  }, numeric(1))

  bounds <- switch(
    ci_method,
    wald       = list(lower = wald_low,  upper = wald_upp),
    percentile = list(lower = perc_low,   upper = perc_upp),
    basic      = list(lower = basic_low,  upper = basic_upp),
    bc         = list(lower = bc_low,     upper = bc_upp),
    normal     = list(lower = normal_low, upper = normal_upp)
  )

  list(
    estimate = est,
    lower = bounds$lower,
    upper = bounds$upper,
    method = ci_method,
    lower_wald = wald_low,
    upper_wald = wald_upp,
    lower_percentile = perc_low,
    upper_percentile = perc_upp,
    lower_basic = basic_low,
    upper_basic = basic_upp,
    lower_bc = bc_low,
    upper_bc = bc_upp,
    lower_normal = normal_low,
    upper_normal = normal_upp,
    boot_mean = boot_mean,
    boot_sd = boot_sd
  )
}

gen_bootstrap_y <- function(D_mean, r_hat, Z, W, V_star,
                            rho_vals, Phi, grids, eta,
                            rho_clip = 0.99) {
  N <- nrow(Z);  L <- length(grids)
  rho_safe         <- stabilize_rho_path(rho_vals, rho_clip = rho_clip)
  V_boot           <- V_star * eta
  XD_eval          <- Z %*% D_mean %*% Phi

  if (W_can_be_sparse && N >= 50) {
    W_sp   <- Matrix::Matrix(W, sparse = TRUE)
    I_N    <- Matrix::Diagonal(N)
    y_boot <- matrix(0, N, L)
    for (j in 1:L) {
      A_j         <- I_N - rho_safe[j] * W_sp
      y_boot[, j] <- as.vector(XD_eval[, j] + Matrix::solve(A_j, V_boot[, j]))
    }
  } else {
    y_boot <- matrix(0, N, L)
    for (j in 1:L) {
      S_j <- diag(N) - rho_safe[j] * W
      y_boot[, j] <- XD_eval[, j] + solve(S_j, V_boot[, j])
    }
  }
  return(y_boot)
}

compute_F_test_components <- function(E_FM, E_RM, W, basis, M,
                                      r_FM, r_RM,
                                      eta_FM = NULL, eta_RM = NULL,
                                      rho_bound = NULL, q = NULL,
                                      N = nrow(E_FM), p = NULL) {
  ISE       <- compute_ISE(E_FM, W, r_FM, basis, M, eta = eta_FM, rho_bound = rho_bound)
  ISEc      <- compute_ISE(E_RM, W, r_RM, basis, M, eta = eta_RM, rho_bound = rho_bound)
  delta_ISE <- ISEc - ISE
  F_stat    <- if (!is.null(q) && !is.null(p) && q > 0) {
    denom <- pmax(ISE / (N - p), .Machine$double.eps)
    (delta_ISE / q) / denom
  } else {
    NA_real_
  }
  list(
    statistic = F_stat,
    Fn = F_stat,
    delta_ISE = delta_ISE,
    legacy_Ftype = F_stat,
    ISE = ISE,
    ISEc = ISEc
  )
}

################################################################################
## 4. F-type Test for Linear Hypotheses
################################################################################

bootstrap_Ftest <- function(y, Z, W, seq_obs, basis, M, B, j_test = NULL,
                            constraint_C = NULL, constraint_c = NULL,
                            D_FM = NULL, r_FM = NULL, eta_FM = NULL, ncores = 1,
                            L_boot = 500, M_boot = 500,
                            multiplier = c("mammen", "rademacher", "gaussian"),
                            max_boot_retry = 5L) {
  multiplier     <- match.arg(multiplier)
  max_boot_retry <- as.integer(max_boot_retry)
  N              <- nrow(y)
  p              <- ncol(Z)
  K              <- basis$nbasis

  C <- t(smooth.basis(seq_obs, t(y), basis)$fd$coef)
  constraint_info <- prepare_linear_constraint_info(
    p = p, K = K, j_test = j_test,
    constraint_C = constraint_C, constraint_c = constraint_c
  )
  q <- constraint_info$q

  rho_bound <- get_fsar_w_cache(W)$rho_bound

  if (is.null(D_FM) || is.null(r_FM)) {
    cat("Estimating full model...\n")
    res_FM <- fsar_func_fast(y = y, Z = Z, W = W, seq = seq_obs,
                             basis = basis, M = M, tolerance = 1e-4)
    D_FM   <- matrix(res_FM$D, nrow = p, ncol = K, byrow = FALSE)
    r_FM   <- matrix(res_FM$r, nrow = 1)
    eta_FM <- matrix(res_FM$eta, nrow = 1)
  }
  if (is.null(eta_FM) || length(as.numeric(eta_FM)) == 0) {
    eta_FM <- eta_from_rho_coef(r_FM, eval.basis(seq_obs, basis), rho_bound)
  }

  cat("Estimating reduced model...\n")
  res_RM <- estimate_restricted_sfsar(
    y = y, Z = Z, W = W, seq_obs = seq_obs,
    basis = basis, M = M, tolerance = 1e-4,
    constraint_info = constraint_info
  )
  D_RM         <- res_RM$D_full
  D_RM_reduced <- res_RM$D_reduced
  r_RM         <- res_RM$r
  eta_RM       <- res_RM$eta
  if (is.null(eta_RM) || length(as.numeric(eta_RM)) == 0) {
    eta_RM <- eta_from_rho_coef(r_RM, eval.basis(seq_obs, basis), rho_bound)
  }

  E_FM  <- C - Z %*% D_FM
  E_RM  <- C - Z %*% D_RM
  obs_components <- compute_F_test_components(
    E_FM, E_RM, W, basis, M,
    r_FM = r_FM, r_RM = r_RM,
    eta_FM = eta_FM, eta_RM = eta_RM,
    rho_bound = rho_bound, q = q, N = N, p = p
  )
  Fn_obs <- obs_components$statistic

  res_info <- compute_residuals_on_grid(
    D_RM, r_RM, Z, W, C, basis, grids = seq_obs,
    eta = eta_RM, rho_bound = rho_bound
  )

  one_F_boot <- function(b) {
    attempts <- 0L
    repeat {
      attempts <- attempts + 1L
      eta <- draw_bootstrap_multipliers(N, method = multiplier)
      boot_try <- tryCatch({
        y_boot <- gen_bootstrap_y(D_RM, r_RM, Z, W, res_info$V_star,
                                  res_info$rho_vals, res_info$Phi,
                                  res_info$grids, eta)
        C_boot <- t(smooth.basis(res_info$grids, t(y_boot), basis)$fd$coef)
        rb_FM  <- fsar_func_fast(y = y_boot, Z = Z,    W = W, seq = res_info$grids,
                                basis = basis, M = M_boot, tolerance = 5e-4,
                                D_init = as.numeric(D_FM), r_init = as.numeric(r_FM),
                                eta_init = eta_FM)
        rb_RM <- estimate_restricted_sfsar(
          y = y_boot, Z = Z, W = W, seq_obs = res_info$grids,
          basis = basis, M = M_boot, tolerance = 5e-4,
          constraint_info = constraint_info,
          D_init_reduced = as.numeric(D_RM_reduced),
          r_init = as.numeric(r_RM), eta_init = eta_RM
        )
        Db_FM   <- matrix(rb_FM$D, nrow = p,     ncol = K, byrow = FALSE)
        rb_FM_r <- matrix(rb_FM$r, nrow = 1)
        Db_RM   <- rb_RM$D_full
        rb_RM_r <- rb_RM$r
        Eb_FM   <- C_boot - Z %*% Db_FM
        Eb_RM   <- C_boot - Z %*% Db_RM
        comp_b  <- compute_F_test_components(
          Eb_FM, Eb_RM, W, basis, M_boot,
          r_FM = rb_FM_r, r_RM = rb_RM_r,
          eta_FM = rb_FM$eta, eta_RM = rb_RM$eta,
          rho_bound = rho_bound, q = q, N = N, p = p
        )
        stat_b <- comp_b$statistic
        valid_b <- isTRUE(rb_FM$converged) && isTRUE(rb_RM$converged) &&
          is.finite(stat_b)
        list(
          stat         = stat_b,
          valid        = valid_b,
          attempts     = attempts,
          error        = NULL,
          fm_converged = isTRUE(rb_FM$converged),
          rm_converged = isTRUE(rb_RM$converged),
          delta_ISE    = comp_b$delta_ISE,
          legacy_Ftype = comp_b$legacy_Ftype
        )
      }, error = function(e) {
        list(
          stat         = NA_real_,
          valid        = FALSE,
          attempts     = attempts,
          error        = conditionMessage(e),
          fm_converged = FALSE,
          rm_converged = FALSE,
          delta_ISE    = NA_real_,
          legacy_Ftype = NA_real_
        )
      })

      if (isTRUE(boot_try$valid) || attempts >= max_boot_retry) {
        return(boot_try)
      }
    }
  }

  F_boot_list       <- run_bootstrap_jobs(B, one_F_boot, ncores = ncores, progress_label = "F-test bootstrap")
  Fn_boot           <- vapply(F_boot_list, function(x) x$stat, numeric(1))
  boot_valid        <- vapply(F_boot_list, function(x) isTRUE(x$valid), logical(1))
  boot_attempts     <- vapply(F_boot_list, function(x) x$attempts, integer(1))
  boot_delta_ISE    <- vapply(F_boot_list, function(x) x$delta_ISE, numeric(1))
  boot_legacy_Ftype <- vapply(F_boot_list, function(x) x$legacy_Ftype, numeric(1))
  boot_fm_converged <- vapply(F_boot_list, function(x) x$fm_converged, logical(1))
  boot_rm_converged <- vapply(F_boot_list, function(x) x$rm_converged, logical(1))
  n_valid           <- sum(boot_valid)
  if (n_valid == 0L) {
    stop("F-test bootstrap produced no valid replicates.")
  }

  p_value <- mean(Fn_boot[boot_valid] >= Fn_obs)
  return(list(Fn = Fn_obs, p_value = p_value, Fn_boot = Fn_boot,
              statistic         = Fn_obs, stat_name = "F_eq14",
              delta_ISE         = obs_components$delta_ISE,
              legacy_Ftype      = obs_components$legacy_Ftype,
              ISE_full          = obs_components$ISE,
              ISE_restricted    = obs_components$ISEc,
              D_FM = D_FM, r_FM = r_FM, eta_FM = eta_FM,
              D_RM = D_RM, r_RM = r_RM, eta_RM = eta_RM,
              constraint_C      = constraint_info$C,
              constraint_c      = constraint_info$c_coef,
              rho_bound         = rho_bound,
              multiplier        = multiplier,
              B_valid           = n_valid,
              boot_valid        = boot_valid,
              boot_attempts     = boot_attempts,
              boot_delta_ISE    = boot_delta_ISE,
              boot_legacy_Ftype = boot_legacy_Ftype,
              boot_fm_converged = boot_fm_converged,
              boot_rm_converged = boot_rm_converged,
              max_boot_retry    = max_boot_retry))
}

################################################################################
## 5. Pointwise Confidence Intervals 
################################################################################

bootstrap_CI <- function(y, Z, W, seq_obs, basis, M, B, alpha = 0.05,
                         D_hat = NULL, r_hat = NULL, eta_hat = NULL, ncores = 1,
                         L_boot = 500, M_boot = 500,
                         ci_method = c("wald", "normal", "bc", "percentile", "basic"),
                         multiplier = c("gaussian", "rademacher", "mammen"),
                         rho_clip = 0.99) {
  ci_method  <- match.arg(ci_method)
  multiplier <- match.arg(multiplier)
  N          <- nrow(y)
  p          <- ncol(Z)
  K          <- basis$nbasis
  rho_bound  <- get_fsar_w_cache(W)$rho_bound

  C <- t(smooth.basis(seq_obs, t(y), basis)$fd$coef)

  if (is.null(D_hat) || is.null(r_hat)) {
    cat("Estimating SFSAR model...\n")
    result <- fsar_func_fast(y = y, Z = Z, W = W, seq = seq_obs,
                             basis = basis, M = M, tolerance = 1e-4)
    D_hat   <- matrix(result$D, nrow = p, ncol = K, byrow = FALSE)
    r_hat   <- matrix(result$r, nrow = 1)
    eta_hat <- matrix(result$eta, nrow = 1)
  }
  if (is.null(eta_hat) || length(as.numeric(eta_hat)) == 0) {
    eta_hat <- eta_from_rho_coef(r_hat, eval.basis(seq_obs, basis), rho_bound)
  }


  res_info <- compute_residuals_on_grid(
    D_hat, r_hat, Z, W, C, basis, grids = seq_obs,
    eta = eta_hat, rho_bound = rho_bound
  )

  # One bootstrap replicate using fsar_func_fast with warm-start
  one_CI_boot <- function(b) {
    tryCatch({
      eta    <- draw_bootstrap_multipliers(N, method = multiplier)
      y_boot <- gen_bootstrap_y(D_hat, r_hat, Z, W, res_info$V_star,
                                res_info$rho_vals, res_info$Phi,
                                res_info$grids, eta, rho_clip = rho_clip)
      res_b  <- fsar_func_fast(y = y_boot, Z = Z, W = W, seq = res_info$grids,
                              basis = basis, M = M_boot, tolerance = 5e-4,
                              D_init = as.numeric(D_hat), r_init = as.numeric(r_hat),
                              eta_init = eta_hat)
      list(
        D_b       = matrix(res_b$D, nrow = p, ncol = K, byrow = FALSE),
        r_b       = as.numeric(res_b$r),
        eta_b     = as.numeric(res_b$eta),
        converged = isTRUE(res_b$converged),
        iter      = res_b$iter,
        beta_dif  = res_b$beta_dif,
        rho_dif   = res_b$rho_dif,
        error     = NULL
      )
    }, error = function(e) {
      list(
        D_b       = matrix(NA_real_, nrow = p, ncol = K),
        r_b       = rep(NA_real_, K),
        eta_b     = rep(NA_real_, K),
        converged = FALSE,
        iter      = NA_integer_,
        beta_dif  = NA_real_,
        rho_dif   = NA_real_,
        error     = conditionMessage(e)
      )
    })
  }

  D_boot_all     <- array(0, dim = c(B, p, K))
  r_boot_all     <- matrix(0, nrow = B, ncol = K)
  eta_boot_all   <- matrix(0, nrow = B, ncol = K)
  boot_converged <- logical(B)
  boot_iter      <- integer(B)
  boot_beta_dif  <- rep(NA_real_, B)
  boot_rho_dif   <- rep(NA_real_, B)
  boot_list      <- run_bootstrap_jobs(B, one_CI_boot, ncores = ncores, progress_label = "CI bootstrap")
  for (b in 1:B) {
    D_boot_all[b, , ] <- boot_list[[b]]$D_b
    r_boot_all[b, ]   <- boot_list[[b]]$r_b
    eta_boot_all[b, ] <- boot_list[[b]]$eta_b
    boot_converged[b] <- boot_list[[b]]$converged
    boot_iter[b]      <- boot_list[[b]]$iter
    boot_beta_dif[b]  <- boot_list[[b]]$beta_dif
    boot_rho_dif[b]   <- boot_list[[b]]$rho_dif
  }
  boot_valid <- boot_converged &
    apply(D_boot_all, 1, function(x) all(is.finite(x))) &
    apply(eta_boot_all, 1, function(x) all(is.finite(x)))
  if (sum(boot_valid) < 2L) {
    stop("CI bootstrap produced fewer than two valid replicates.")
  }

  eval_grids <- seq(basis$rangeval[1], basis$rangeval[2], length.out = 1000)
  Phi_eval   <- eval.basis(eval_grids, basis)       # 1000 × K

  beta_ci <- list()
  for (j in 1:p) {
    est_j        <- as.numeric(Phi_eval %*% D_hat[j, ])
    boot_j       <- tcrossprod(D_boot_all[boot_valid, j, ], Phi_eval)   # B_valid × 1000
    beta_ci[[j]] <- compute_pointwise_bootstrap_ci(
      est = est_j, boot_mat = boot_j, alpha = alpha, ci_method = ci_method
    )
  }

  rho_est <- if (!is.null(eta_hat)) {
    eval_bounded_rho_from_eta(Phi_eval, eta_hat, rho_bound)
  } else {
    as.numeric(Phi_eval %*% t(r_hat))
  }
  rho_boot <- if (sum(boot_valid) > 0) {
    rho_bound * tanh(eta_boot_all[boot_valid, , drop = FALSE] %*% t(Phi_eval))
  } else {
    r_boot_all[boot_valid, , drop = FALSE] %*% t(Phi_eval)
  }
  rho_ci <- compute_pointwise_bootstrap_ci(
    est = rho_est, boot_mat = rho_boot, alpha = alpha, ci_method = ci_method
  )

  return(list(beta_ci = beta_ci, rho_ci = rho_ci,
              D_hat          = D_hat, r_hat = r_hat, eta_hat = eta_hat,
              D_boot         = D_boot_all, r_boot = r_boot_all, eta_boot = eta_boot_all,
              boot_valid     = boot_valid,
              B_valid        = sum(boot_valid),
              boot_converged = boot_converged,
              boot_iter      = boot_iter,
              boot_beta_dif  = boot_beta_dif,
              boot_rho_dif   = boot_rho_dif,
              eval_grids     = eval_grids,
              ci_method      = ci_method,
              multiplier     = multiplier,
              rho_clip       = rho_clip,
              rho_bound      = rho_bound))
}
