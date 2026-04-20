
#### fos(): Function-on-scalar Regression#################
## INPUT:
# y:         observed response points
# Z:         design matrix
# seq:       observed points
# basis:     b-spline basis/fourier basis
# M:         divide the range into M points to approximate integrals
# ncores:    the number of cores to use, must be exactly 1 on Windows
## OUTPUT:
# list(D), where
# D:        vec( p*K coefficient matrix of basis of functional slope beta(t) )
####

fos = function(y, Z, seq, basis, M, ncores){ 
  t1      = Sys.time()
  N       = dim(Z)[1]
  p       = dim(Z)[2]
  K       = basis$nbasis
  ySmooth = smooth.basis(seq, t(y), basis) # smoothed y(t), class=fdSmooth
  yfd     = ySmooth$fd # functional object y(t)
  coef    = yfd$coef   # K*N coefficient matrix of basis of response y(t)
  C       = t(coef)    # N*K coefficient matrix of basis of response y(t)
  # parameter to be estimated
  D       = matrix(0,   nrow = K*p, ncol = 1) # vectorized p*K coefficient matrix of basis of functional slope beta(t)

  grids    = seq(basis$rangeval[1],basis$rangeval[2],  length.out = M)
  
  D_estimate = function(t){
    phi   = t(eval.basis(t, basis))
    A_t   = t(Z) %*%  Z
    B_t   = phi  %*% t(phi)
    E_t   = t(Z) %*%  C %*% B_t
    BkA_t = kronecker(t(B_t), A_t)
    list(BkA_t = BkA_t, E_t = E_t)
  }
  D_result = mclapply(grids, D_estimate, mc.cores = ncores)
  BkA      = mclapply(D_result, function(x) x[[1]], mc.cores = ncores)
  E        = mclapply(D_result, function(x) x[[2]], mc.cores = ncores)
  BkA      = Reduce(`+`, BkA) * as.numeric(grids[2] - grids[1])
  E        = Reduce(`+`, E)   * as.numeric(grids[2] - grids[1])
  E_vec    = as.vector(E)
  BkA_inv  = solve(BkA)
  D_new    = BkA_inv %*% E_vec
  D        = D_new
  return(list(D))
}

