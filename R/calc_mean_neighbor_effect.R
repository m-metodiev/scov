#' Calculates the average correlation for the spatial effect
#'
#' @inheritParams avg_effect
#' @param beta                autocorrelation parameter of the CAR model
#' @param delta               linear effect of the CAR model correlation matrix
#' @inheritParams scov
#'
#' @returns the average correlation for the spatial effect
#' @keywords internal
#'
calc_mean_neighbor_effect = function(matList, beta, delta, adj_positions){
  G_inv = tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],beta)[adj_positions,
                                                            adj_positions]
  A = matList$Al[[1]][adj_positions,adj_positions]
  A[is.na(A)] = 0
  return(mean(delta*G_inv[as.matrix(A!=0)]))
}
