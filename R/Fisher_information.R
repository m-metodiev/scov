#' Computes the Fisher information matrix
#'
#' @param adj_positions       positions of the spatial effect (if embedded)
#' @param parm                the parameter
#' @param matList             list of the matrices (pairwise+spatial)
#' @param interaction_effects list of pairwise effects (vectors of named pairs)
#'
#' @returns the Fisher information matrix
#' @keywords internal
Fisher_information = function(adj_positions, parm, matList,
                              interaction_effects=list()){
  Sigma = CovMat_03(parm,
                    matList,
                    adj_positions=adj_positions,
                    interaction_effects=interaction_effects)$Sigma
  Sigma_inv = solve(Sigma)
  Sigma_der = GradLogLikParm_02(adj_positions,
                                parm,
                                matList,
                                dataset=matrix(0,ncol=dim(Sigma)[1],
                                               nrow=dim(Sigma)[1]),
                                interaction_effects=interaction_effects,
                                return_Sigma_der=TRUE)
  Fisher_mat = matrix(0, ncol=length(parm), nrow=length(parm))
  for(i in seq_along(parm)){
    for(j in seq_along(parm)){
      Fisher_mat[i,j] = (1/2)*sum(diag(Sigma_inv%*%
                                         Sigma_der[[i]]%*%
                                         Sigma_inv%*%
                                         Sigma_der[[j]]))
    }
  }
  return(Fisher_mat)
}
