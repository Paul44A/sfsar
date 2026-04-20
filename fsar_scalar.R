
#### fsar_scalar(): Functional Simultaneous Auto-regressive Regression with ####
## INPUT:
# y:         observed response points
# Z:         design matrix
# W:         spatial weight matrix 
# seq:       observed points
# basis:     b-spline basis/fourier basis
# M:         divide the range into M points to approximate integrals
# tolerance: allowed to identify convergence
## OUTPUT:
# list(D,r,iter), where
# D:        vec( p*K coefficient matrix of basis of functional slope beta(t) )
# r:        1*1  rho
# iter:     number of iterations
################################################################################

fsar_scalar = function(y, Z, W, seq, basis, M, tolerance, ncores){ 
  N       = dim(W)[1]
  p       = dim(Z)[2]
  K       = basis$nbasis
  C       = t(smooth.basis(seq, t(y), basis)$fd$coef)  # N x K
  G       = inprod(basis, basis)                       # K x K Gram matrix

  eigvals         <- eigen(W, only.values = TRUE)$values
  spectral_radius <- max(Mod(eigvals))
  row_abs_sum_max <- max(rowSums(abs(W)))
  rho_row_bound   <- if (is.finite(row_abs_sum_max) && row_abs_sum_max > 0) {
    0.99 / row_abs_sum_max
  } else {
    Inf
  }
  rho_spec_bound <- if (is.finite(spectral_radius) && spectral_radius > 0) {
    0.99 / spectral_radius
  } else {
    0
  }
  rho_bound <- if (is.finite(rho_row_bound) && rho_row_bound > 0) {
    min(rho_row_bound, rho_spec_bound)
  } else {
    rho_spec_bound
  }

  ZtZ    <- crossprod(Z)
  Ws     <- W + t(W)
  WtW    <- crossprod(W)
  ZtWsZ  <- crossprod(Z, Ws %*% Z)
  ZtWtWZ <- crossprod(Z, WtW %*% Z)
  ZtC    <- crossprod(Z, C)
  ZtWsC  <- crossprod(Z, Ws %*% C)
  ZtWtWC <- crossprod(Z, WtW %*% C)

  sse_floor <- max(.Machine$double.eps,
                   1e-10 * sum(abs(crossprod(C) * G)))
  if (!is.finite(sse_floor) || sse_floor <= 0) {
    sse_floor <- 1e-8
  }

  profile_for_rho <- function(rho) {
    A_r <- ZtZ - rho * ZtWsZ + rho^2 * ZtWtWZ
    B_r <- ZtC - rho * ZtWsC + rho^2 * ZtWtWC
    D_r <- tryCatch(solve(A_r, B_r), error = function(e) NULL)
    if (is.null(D_r)) {
      return(NULL)
    }
    E_r <- C - Z %*% D_r
    sse <- sum((crossprod(E_r) - rho * (t(E_r) %*% Ws %*% E_r) +
                  rho^2 * (t(E_r) %*% WtW %*% E_r)) * G)
    if (!is.finite(sse)) {
      return(NULL)
    }
    list(D = D_r, E = E_r, sse = sse)
  }

  objective <- function(rho) {
    prof <- profile_for_rho(rho)
    if (is.null(prof)) {
      return(1e12)
    }
    sse_safe   <- max(prof$sse, sse_floor)
    logdet_abs <- sum(log(pmax(Mod(1 - rho * eigvals), sse_floor)))
    log(sse_safe / N) - (2 / N) * logdet_abs
  }

  if (rho_bound <= 1e-8) {
    rho_hat <- 0
    prof_hat <- profile_for_rho(rho_hat)
    opt_obj <- objective(rho_hat)
  } else {
    rho_eps <- min(1e-6, rho_bound / 1000)
    opt_res <- optimize(
      f        = objective,
      interval = c(-rho_bound + rho_eps, rho_bound - rho_eps),
      tol      = max(tolerance, 1e-6)
    )
    rho_hat  <- opt_res$minimum
    prof_hat <- profile_for_rho(rho_hat)
    opt_obj  <- opt_res$objective
  }

  if (is.null(prof_hat)) {
    stop("fsar_scalar QMLE failed to obtain a valid profiled estimate.")
  }

  return(list(
    D = as.numeric(prof_hat$D),
    r = as.numeric(rho_hat),
    iter = 1,
    converged = TRUE,
    objective = opt_obj,
    rho_bound = rho_bound
  ))
}
