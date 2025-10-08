#' Computes the inverse of the correlation matrix of the CAR model
#'
#' @param M                   used to define the CAR correlation matrix
#' @param A                   used to define the CAR correlation matrix
#' @param beta                autocorrelation parameter of the CAR model
#' @param U_full              (optional) can be used to fasten computation
#' @param solve_U_full        (optional) can be used to fasten computation
#' @param solve_M_no_islands  (optional) can be used to fasten computation
#' @param eigen_real          (optional) can be used to fasten computation
#' @param return_U_D_M        (optoinal) can be used to fasten computation
#'
#' @returns the inverse of the correlation matrix of the CAR model
#' @keywords internal
#' @importFrom stats          cov2cor
#' @importFrom pracma         nullspace
#' @importFrom Matrix         rankMatrix
tilde_G_inv = function(M=NULL, A=NULL, beta,
                       U_full=NULL, solve_U_full=NULL,
                       solve_M_no_islands=NULL, eigen_real=NULL,
                       return_U_D_M=FALSE){
  # round down if too close to the edge
  if(beta==1){
    beta=.999
  }
  beta = beta*(beta<1-1e-4)+1e-4*(beta>=1-1e-4)

  # set diagonal to 1 for isolated nodes (A can't be defined for those)
  no_islands_id = diag(M)!=0
  G = matrix(0,ncol=dim(M)[1],nrow=dim(M)[1])

  # it is possible to save U_full to quicken the computation
  M = as.matrix(M)
  A = as.matrix(A)

  M_no_islands = M[no_islands_id,no_islands_id]
  A_no_islands = A[no_islands_id,no_islands_id]
  dimension = dim(A_no_islands)[1]

  # Matrix is only positive semidefinite (the are eigenvalues equal to 0),
  # due to the graph not being connected
  rankA = Matrix::rankMatrix(A_no_islands)[1]

  if(is.null(U_full)){
    U_real = Re(eigen(A_no_islands)$vectors[,1:rankA])
    N_real = pracma::nullspace(A_no_islands)
    U_full = cbind(U_real,N_real)
    eigen_real = Re(eigen(A_no_islands)$values[1:rankA])

    D_beta_full = diag(c(eigen_real * beta / (1 - eigen_real * beta),
                         rep(0,dimension-rankA)))

    G_inv_no_islands = (diag(dimension) +
                          U_full %*%
                          D_beta_full %*%
                          solve(U_full)) %*%
      solve(M_no_islands)
  } else{
    D_beta_full = diag(c(eigen_real * beta / (1 - eigen_real * beta),
                         rep(0,dimension-rankA)))

    G_inv_no_islands = (diag(dimension) + U_full %*%
                          D_beta_full %*%
                          solve_U_full) %*%
      solve_M_no_islands

  }

  G[no_islands_id,no_islands_id] = stats::cov2cor(G_inv_no_islands)
  diag(G)=1 # islands are independent of all other nodes

  if(return_U_D_M){ # to compute the real value
    return(list(U_full=U_full, solve_U_full=solve(U_full),
                solve_M_no_islands=solve(M_no_islands), eigen_real=eigen_real))
  } else{
    return(G)
  }
}
