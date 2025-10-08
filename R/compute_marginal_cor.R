#' Computes non-positive-semidefinite approximation of correlation matrix
#'
#' @param dataset n x d matrix (n = number of observations, d = dimension)
#'
#' @returns non-positive-semidefinite approximation of correlation matrix
#' @keywords internal
compute_marginal_cor = function(dataset){
  dimension = ncol(dataset)
  correlation_matrix = matrix(ncol=dimension,nrow=dimension)

  #use pairwise correlation estimates
  for(i in (1:dimension)){
    for(j in (i:dimension)){
      dataseti_notmissing = !is.na(dataset[,i])
      datasetj_notmissing = !is.na(dataset[,j])
      correlation_matrix_ij = cor(dataset[dataseti_notmissing&
                                            datasetj_notmissing,i],
                                  dataset[dataseti_notmissing&
                                            datasetj_notmissing,j])
      correlation_matrix[i,j] = correlation_matrix_ij
      correlation_matrix[j,i] = correlation_matrix_ij
    }
  }
  return(correlation_matrix)
}
