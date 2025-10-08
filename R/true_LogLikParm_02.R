#' Computes the "true" (i.e., not translated) log-likelihood (needed for BIC)
#'
#' @param adj_positions         positions of spatial effect (if embedded)
#' @param parm                  the parameter
#' @param matList               list of matrices (pairwise + spatial)
#' @param dataset               an n (observations) x d (dimension) matrix
#' @param interaction_effects   list of interaction effects (vector of names)
#'
#' @returns the "true" (i.e., not translated) log-likelihood
#' @keywords internal
#' @importFrom mvtnorm          dmvnorm
true_LogLikParm_02 <- function(adj_positions, parm, matList, dataset,
                               interaction_effects=list()){
  num_observations = nrow(dataset)
  has_missing_values = sum(is.na(dataset))>0
  res = 0

  for(t in (1:num_observations)){
    adj_pos_t = !is.na(dataset[t,])
    res =  res +
      mvtnorm::dmvnorm(
        dataset[t,adj_pos_t],
        sigma=as.matrix(CovMat_03(
          adj_positions=adj_positions,
          parm=parm,
          matList=matList,
          interaction_effects=interaction_effects)$Sigma)[adj_pos_t,
                                                          adj_pos_t], log=TRUE)
  }
  return(res)
}
