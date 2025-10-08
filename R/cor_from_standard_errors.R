#' Computes correlation matrix from a normalized dataset (=standard errors)
#'
#' @param varepsilon n x d matrix (n = number of observations, d = dimension)
#'
#' @returns the correlation matrix
#' @keywords internal
#' @importFrom Matrix          nearPD
cor_from_standard_errors = function(varepsilon){
  dimension=dim(varepsilon)[2]
  num_observations = dim(varepsilon)[1]
  correlation_matrix = matrix(0,ncol=dimension,nrow=dimension)
  for(t in 1:num_observations){
    correlation_matrix = correlation_matrix +
      t(t(varepsilon[t,]))%*%t(varepsilon[t,])
  }
  correlation_matrix = correlation_matrix/(num_observations-1)
  diag(correlation_matrix)=1

  #Find nearest positive definite matrix
  correlation_matrix = Matrix::nearPD(correlation_matrix,
                                      corr=TRUE,
                                      conv.tol=100000)$mat
  return(as.matrix(correlation_matrix))
}
