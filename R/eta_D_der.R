#' Calculates the derivative of the spatial average effect
#'
#' @param parm                the parameter
#' @param matList             a list of the matrices (pairwise+spatial)
#' @param adj_positions       the positions of the spatial effect (if embedded)
#' @param interaction_effects the interaction effect pairs (2-pair-name-vectors)
#' @param index               position of the spatial effect within the list
#'
#' @returns the derivative of the spatial average effect
#' @keywords internal
eta_D_der = function(parm, matList, adj_positions,
                     #link_der_beta,
                     interaction_effects=list(),index=4){
  #browser()
  parm=sapply(c(parm),function(p) p)
  covMatstuff = CovMat_03(parm=parm ,
                          matList=matList,
                          adj_positions=adj_positions,
                          interaction_effects=interaction_effects)
  ml_combined = covMatstuff$matList_combined
  Sigma = CovMat_03(parm, matList,adj_positions=adj_positions,
                    interaction_effects=interaction_effects)$Sigma
  Sigma_der = GradLogLikParm_02(adj_positions,
                                parm,
                                matList,
                                dataset=matrix(0,ncol=dim(Sigma)[1],
                                               nrow=dim(Sigma)[1]),
                                interaction_effects=interaction_effects,
                                return_Sigma_der=TRUE)

  # determine support
  matList_supp = matList

  #Set to 1 if countries are not neighbors
  matList_supp$Al[[1]][is.na(matList_supp$Al[[1]])] = 0
  matList_supp$Gl[[1]] = (matList_supp$Al[[1]][adj_positions,
                                               adj_positions] != 0) + 0
  ml_combined_supp = combined_matList(matList_supp,
                                      interaction_effects)$matList_full
  #diagonals are not in the sum
  diag(ml_combined[[index]]) = 0
  diag(ml_combined_supp[[index]]) = 0

  eta_D_der = c(sum(parm[index]*ml_combined[[index]]*ml_combined_supp[[index]]),
                sum(Sigma_der[[index]]*ml_combined_supp[[index]])) /
    sum(ml_combined_supp[[index]])
  return(eta_D_der)
}
