#' Estimates the correlation matrix of the dataset
#'
#' @param dataset     n x d matrix (n = number of observations, d = dimension)
#' @param mean_estim  estimate of the mean vector of the dataset
#' @param sd_estim    estimate of the standard deviation vector of the dataset
#'
#' @returns an estimate of the correlation matrix
#' @keywords internal
correlation_matrix = function(dataset, mean_estim, sd_estim){
  num_observations = nrow(dataset)
  dimension = ncol(dataset)
  if(sum(is.na(dataset))==0){
    normalized_dataset = (dataset - t(matrix(rep(mean_estim,num_observations),
                                             ncol=num_observations)))%*%
      diag(1/sd_estim)

    return(cor_from_standard_errors(normalized_dataset))
  } else{
    normalized_dataset = sapply(seq_along(dataset[1,]),
                                function(s) (dataset[,s] - mean_estim[s])/
                                  sd_estim[s])
    correlation_matrix = compute_marginal_cor(normalized_dataset)
    return(correlation_matrix)
  }
}
