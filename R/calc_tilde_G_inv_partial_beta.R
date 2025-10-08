#' Calculates the derivative of the inverse of G (the CAR model matrix)
#'
#' @param M     matrix used in the definition of the CAR model
#' @param A     matrix used in the definition of the CAR model
#' @param beta  autocorrelation parameter of the CAR model
#' @keywords internal
#' @importFrom Matrix         rankMatrix
#' @importFrom pracma         nullspace
#' @importFrom stats          cov2cor
#'
#' @returns the derivative of the inverse of G (the CAR model matrix)
calc_tilde_G_inv_partial_beta = function(M, A, beta){

  M = as.matrix(M)
  A = as.matrix(A)

  # set diagonal to 1 for isolated nodes (A can't be defined for those)
  no_islands_id = diag(M)!=0
  G_inv_partial_beta = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])
  G_inv = diag(dim(M)[1])

  # Matrix is only positive semidefinite due to the graph not being connected
  M_no_islands = M[no_islands_id,no_islands_id]
  A_no_islands = A[no_islands_id,no_islands_id]
  dimension = dim(A_no_islands)[1]
  rankA = Matrix::rankMatrix(A_no_islands)[1]
  U_real = Re(eigen(A_no_islands)$vectors[,1:rankA])
  N_real = pracma::nullspace(A_no_islands)
  U_full = cbind(U_real,N_real)
  eigen_real = Re(eigen(A_no_islands)$values[1:rankA])
  D_beta_full = diag(c(eigen_real * beta / (1 - eigen_real * beta),
                       rep(0,dimension-rankA)))
  D_beta_full_partial_beta = diag(c(eigen_real/ ((1 - eigen_real * beta)^2),
                                    rep(0,dimension-rankA)))
  solve_U = solve(U_full)
  solve_M = solve(M_no_islands)

  G_inv_partial_beta_no_islands =
    (U_full%*%D_beta_full_partial_beta%*%solve_U)%*%solve_M
  G_inv_partial_beta[no_islands_id,no_islands_id] =
    G_inv_partial_beta_no_islands
  G_inv_no_islands = (diag(dimension) +
                        U_full%*%D_beta_full%*%solve_U)%*%solve_M
  G_inv[no_islands_id, no_islands_id] = G_inv_no_islands

  S_inv_sqrt_no_islands = diag(1/sqrt(diag(G_inv_no_islands)))
  S_inv_sqrt_partial_beta_no_islands = (-1/2)*
    diag(1/(diag(G_inv_no_islands)^(3/2))*diag(G_inv_partial_beta_no_islands))

  tilde_G_inv_partial_beta_no_islands = S_inv_sqrt_partial_beta_no_islands %*%
    G_inv_no_islands %*%
    S_inv_sqrt_no_islands +
    S_inv_sqrt_no_islands %*%
    G_inv_partial_beta_no_islands %*%
    S_inv_sqrt_no_islands +
    S_inv_sqrt_no_islands %*%
    G_inv_no_islands %*%
    S_inv_sqrt_partial_beta_no_islands

  #islands are constant, so their derivatives are all 0
  tilde_G_inv_partial_beta = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])
  tilde_G_inv_partial_beta[no_islands_id, no_islands_id] =
    tilde_G_inv_partial_beta_no_islands

  # the diagonal values are also constant, so the diagonal derivative is also 0
  diag(tilde_G_inv_partial_beta) = 0

  # return the correlation matrix at the end
  G_inv[no_islands_id, no_islands_id] = stats::cov2cor(G_inv_no_islands)
  return(list(tilde_G_inv = G_inv,
              tilde_G_inv_partial_beta=tilde_G_inv_partial_beta))
}
