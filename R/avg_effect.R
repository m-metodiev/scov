#' Calculate average effects (the mean effect over the matrix support)
#'
#' @param parm                the parameter
#' @param matList             list of pairwise and spatial matrices
#' @inheritParams scov
#'
#' @returns the average effects (the mean effect over the matrix support)
#'
#' @keywords internal
avg_effect = function(parm, matList, adj_positions, interaction_effects=list()){

  #Determine support of the matrices
  matList_supp = matList

  if(!is.null(matList$Al)){
    #Set to 0 if countries are not neighbors
    matList_supp$Al[[1]][is.na(matList_supp$Al[[1]])] = 0
    matList_supp$Gl[[1]] =
      (matList_supp$Al[[1]][adj_positions,adj_positions] != 0) + 0
  }

  ml_combined_supp = combined_matList(matList_supp,
                                      interaction_effects)$matList_full

  covMatstuff = CovMat_03(parm=parm, matList=matList,
                          adj_positions=adj_positions,
                          interaction_effects=interaction_effects)
  ml_combined = covMatstuff$matList_combined
  alpha_delta = covMatstuff$alpha_delta

  for(Fk in ml_combined_supp){
    diag(Fk)=0
  }
  for(Fk in ml_combined){
    diag(Fk)=0
  }
  # comparing the diagonals makes no sense
  return(sapply(seq_along(alpha_delta),
                function(i) alpha_delta[i]*
                  sum((ml_combined_supp[[i]]!=0)*ml_combined[[i]])/
                  sum(ml_combined_supp[[i]]!=0)))
}
