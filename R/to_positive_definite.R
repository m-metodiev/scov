#' Maps the pairwise covariates to symmetric, positive definite matrices
#'
#' @param pairwise_covariate_matrices list of pairwise covariate matrices
#' @param eig.tol treat EV a as negative if (a / largest eigenvalue) <= eig.tol
#'
#' @returns list of pairwise covariate matrices
#' @keywords internal
#' @importFrom Matrix           nearPD
#' @importFrom Matrix           rankMatrix
to_positive_definite = function(pairwise_covariate_matrices, eig.tol=1e-06){
  dimension = nrow(pairwise_covariate_matrices[[1]])
  for(i in seq_along(pairwise_covariate_matrices)){

    # diagonal values always have to be 1 (correlation matrices)
    diag(pairwise_covariate_matrices[[i]])=1

    # check if matrix is symmetric
    if(!isSymmetric(as.matrix(pairwise_covariate_matrices[[i]]))){
      pairwise_covariate_matrices[[i]]=
        Matrix::nearPD(pairwise_covariate_matrices[[i]],
                       corr=TRUE,conv.tol=100000)$mat
    }

    # check if matrix is positive semidefinite
    sorted_eigenvalues = sort(eigen(pairwise_covariate_matrices[[i]])$values)
    # ignore eigenvalues that are exactly 0 (more stable numerically)
    rank_of_matrix = Matrix::rankMatrix(pairwise_covariate_matrices[[i]])[1]
    normalize_zeroes = c(rep(FALSE,dimension-rank_of_matrix),
                        rep(TRUE,rank_of_matrix))
    max_eigen_val = max(sorted_eigenvalues)
    if(min(sorted_eigenvalues[normalize_zeroes]/max_eigen_val) < eig.tol){
      pairwise_covariate_matrices[[i]] =
        as.matrix(Matrix::nearPD(pairwise_covariate_matrices[[i]],
                       corr=TRUE,conv.tol=100000)$mat)
    }
  }
  return(pairwise_covariate_matrices)
}
