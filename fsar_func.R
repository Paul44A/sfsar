#### fsar_func(): Functional Simultaneous Auto-regressive Regression ###########
## INPUT:
# y:         observed response points
# Z:         design matrix
# W:         spatial weight matrix 
# seq:       observed points
# basis:     b-spline basis/fourier basis
# M:         divide the range into M points to approximate integrals
# tolerance: allowed to identify convergence
# ncores:    the number of cores to use, must be exactly 1 on Windows
## OUTPUT:
# list(D,r,iter), where
# D:        vec( p*K coefficient matrix of basis of functional slope beta(t) )
# r:        1*K coefficient matrix of basis of functional spatial autocorrelation coefficient rho(t)
# iter:     number of iterations
####


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

optimize_bounded_rho_step <- function(Phi, itv, q0, q1, q2, rho_bound, eigvals,
                                      N, eta_init, maxit = 100) {
  eta_init <- as.numeric(eta_init)
  sse_floor <- max(.Machine$double.eps,
                   1e-10 * mean(q0[q0 > 0], na.rm = TRUE))
  if (!is.finite(sse_floor)) {
    sse_floor <- 1e-8
  }

  eig_mat <- matrix(eigvals, nrow = nrow(Phi), ncol = length(eigvals),
                    byrow = TRUE)
  last_par <- NULL
  last_eval <- NULL

  evaluate_once <- function(par) {
    if (!is.null(last_par) && max(abs(par - last_par)) < 1e-12) {
      return(last_eval)
    }

    eta_vals  <- as.numeric(Phi %*% par)
    tanh_vals <- tanh(eta_vals)
    rho_vals  <- rho_bound * tanh_vals
    drho_deta <- rho_bound * (1 - tanh_vals^2)

    sse_vals  <- q0 - rho_vals * q1 + rho_vals^2 * q2
    sse_safe  <- pmax(sse_vals, sse_floor)

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

  objective <- function(par) evaluate_once(par)$objective
  gradient  <- function(par) evaluate_once(par)$gradient

  opt <- optim(
    par     = eta_init,
    fn      = objective,
    gr      = gradient,
    method  = "BFGS",
    control = list(maxit = maxit, reltol = 1e-8)
  )

  opt_eval <- evaluate_once(opt$par)
  list(
    eta_coef    = matrix(opt$par, nrow = 1),
    rho_vals    = opt_eval$rho_vals,
    rho_coef    = matrix(project_vals_to_basis_coef(opt_eval$rho_vals, Phi), nrow = 1),
    objective   = opt$value,
    convergence = opt$convergence
  )
}

fsar_func_fast <- function(y, Z, W, seq, basis, M, tolerance,
                           D_init = NULL, r_init = NULL, eta_init = NULL,
                           max_iter = 1000) {
  N <- dim(W)[1];  p <- dim(Z)[2];  K <- basis$nbasis
  C <- t(smooth.basis(seq, t(y), basis)$fd$coef)  

  W_cache <- get_fsar_w_cache(W)
  Ws      <- W_cache$Ws
  WtW     <- W_cache$WtW
  eigvals <- W_cache$eigvals

  ZtZ    <- crossprod(Z)              
  ZtWsZ  <- crossprod(Z, Ws  %*% Z)  
  ZtWtWZ <- crossprod(Z, WtW %*% Z)   
  ZtC    <- crossprod(Z, C)           
  ZtWsC  <- crossprod(Z, Ws  %*% C)  
  ZtWtWC <- crossprod(Z, WtW %*% C)  

  grids        <- seq(basis$rangeval[1], basis$rangeval[2], length.out = M)
  itv          <- grids[2] - grids[1]
  Phi          <- eval.basis(grids, basis)    
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
    E_sum <- ZtC %*% G0 - ZtWsC %*% G1 + ZtWtWC %*% G2
    D_new <- as.numeric(solve(BkA, as.numeric(E_sum)))

    D_mat <- matrix(D_new, nrow = p, ncol = K, byrow = FALSE)
    Res   <- C - Z %*% D_mat
    Q0    <- crossprod(Res)
    WRes  <- W %*% Res
    Q1    <- crossprod(WRes)
    Q2    <- crossprod(Res, Ws %*% Res)

    q0    <- rowSums(Phi * (Phi %*% Q0))
    q1    <- rowSums(Phi * (Phi %*% Q2))
    q2    <- rowSums(Phi * (Phi %*% Q1))

    if (rho_bound > 0) {
      rho_opt <- optimize_bounded_rho_step(
        Phi = Phi, itv = itv,
        q0 = q0, q1 = q1, q2 = q2,
        rho_bound = rho_bound, eigvals = eigvals, N = N,
        eta_init  = eta, maxit = 100
      )
      eta_new      <- rho_opt$eta_coef
      r_new        <- rho_opt$rho_coef
      rho_new_vals <- rho_opt$rho_vals
    } else {
      eta_new      <- matrix(0, nrow = 1, ncol = K)
      r_new        <- matrix(0, nrow = 1, ncol = K)
      rho_new_vals <- rep(0, length(grids))
    }

    beta_dif <- mean(abs(eval.fd(grids,
      fd(t(matrix(D_new - D, p, K, byrow = FALSE)), basis)))) / p
    rho_dif  <- mean(abs(rho_new_vals - rho_v))

    D <- D_new
    r <- r_new
    eta <- eta_new
    if (beta_dif < tolerance && rho_dif < tolerance) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    message("fsar_func_fast: did not fully converge within max_iter; last diff = ",
            round(beta_dif, 4), ", ", round(rho_dif, 4))
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



fsar_func = function(y, Z, W, seq, basis, M, tolerance, ncores){ 
  t1      = Sys.time()
  N       = dim(W)[1]
  p       = dim(Z)[2]
  K       = basis$nbasis
  ySmooth = smooth.basis(seq, t(y), basis) 
  yfd     = ySmooth$fd 
  coef    = yfd$coef   
  C       = t(coef)    #
  D       = matrix(0, nrow = K*p, ncol = 1) 
  r       = matrix(0, nrow = 1,   ncol = K) 
  D_new   = matrix(0, nrow = K*p, ncol = 1) 
  r_new   = matrix(0, nrow = 1,   ncol = K) 
  gc()
  grids    = seq(basis$rangeval[1],basis$rangeval[2],  length.out = M)
  max_iter = 1000; converged = FALSE

  for (iter in 1:max_iter) {
    D_estimate = function(t,r){
      phi   = t(eval.basis(t, basis))
      I_minus_Wrphi = diag(N) - W * as.numeric(r %*% phi)
      A_t   = t(Z) %*% t(I_minus_Wrphi) %*% I_minus_Wrphi %*% Z
      B_t   = phi  %*% t(phi)
      E_t   = t(Z) %*% t(I_minus_Wrphi) %*% I_minus_Wrphi %*% C %*% B_t
      BkA_t = kronecker(t(B_t), A_t)
      list(BkA_t = BkA_t, E_t = E_t)
    }
    D_result = mclapply(grids, D_estimate, r = r_new, mc.cores = ncores)
    BkA      = mclapply(D_result, function(x) x[[1]], mc.cores = ncores)
    E        = mclapply(D_result, function(x) x[[2]], mc.cores = ncores)
    BkA      = Reduce(`+`, BkA) * as.numeric(grids[2] - grids[1])
    E        = Reduce(`+`, E)   * as.numeric(grids[2] - grids[1])
    E_vec    = as.vector(E)
    BkA_inv  = solve(BkA)
    D_new    = BkA_inv %*% E_vec
    
    C_estimate = function(t,d){
      phi        = t(eval.basis(t, basis))
      C_minus_ZD = C - Z %*% matrix(d, nrow = p, ncol = K, byrow = F)
      C1_t       = phi %*% t(phi) %*% t(C_minus_ZD) %*% t(W) %*% W %*% (C_minus_ZD %*% phi %*% t(phi))
      C2_t       = phi %*% t(phi) %*% t(C_minus_ZD) %*% W    %*% (C_minus_ZD %*% phi)
      list(C1_t, C2_t)
    }
    C_result = mclapply(grids, C_estimate, d = D_new, mc.cores = ncores)
    C1_t     = mclapply(C_result, function(x) x[[1]], mc.cores = ncores)
    C2_t     = mclapply(C_result, function(x) x[[2]], mc.cores = ncores)
    C1       = Reduce(`+`, C1_t) * as.numeric(grids[2] - grids[1])
    C2       = Reduce(`+`, C2_t) * as.numeric(grids[2] - grids[1])
    C1_inv   = solve(C1)
    r_new    = t(C1_inv %*% C2)
    
    beta_dif = fd(t(matrix(D_new - D, nrow = p, ncol = K, byrow = F)), basis)
    beta_dif = mean(abs(eval.fd(grids, beta_dif))) / p
    rho_dif  = fd(t(r_new - r), basis)
    rho_dif  = mean(abs(eval.fd(grids, rho_dif)))
    if (beta_dif < tolerance && rho_dif < tolerance) {
      converged = TRUE
      t2        = Sys.time()
      print("---------------------------------------------------------------------")
      print(paste("The",iter,"th", "iteration"))
      print(paste("Time cost:",t2-t1))
      print(paste("The difference of beta(t) from the previous iteration:", beta_dif))
      print(paste("The difference of rho(t)  from the previous iteration:", rho_dif))
      D = D_new
      r = r_new
      gc()
      break
    }
    
    # print and update D & r
    print("---------------------------------------------------------------------")
    print(paste("The",iter,"th", "iteration"))
    print(paste("The difference of beta(t) from the previous iteration:", beta_dif))
    print(paste("The difference of rho(t)  from the previous iteration:", rho_dif))
    D = D_new
    r = r_new
    gc()
    
  }
  
  if (converged) {
    cat("Algorithm converged.\n")
  } else {
    cat("Algorithm did not converge within the maximum number of iterations.\n")
  }
  return(list(D, r, iter))
}
