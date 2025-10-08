#' Calculates the gradient of the function of the transformed parameter
#'
#' @param adj_positions       positions of the spatial effect (if embedded)
#' @param logParm             the transformed parameter
#' @param matList             the list of matrices (pairwise + spatial)
#' @param dataset             an n (observations) x d (dimensions) matrix
#' @param interaction_effects list of interaction effects (vectors of names)
#'
#' @returns the gradient of the function of the transformed parameter
#' @keywords internal
GradLogLikLogParm_02 <- function(adj_positions, logParm, matList, dataset,
                                 interaction_effects=list()){

  parm = backward_transform_param(logParm)

  if(names(parm)[length(names(parm))]=="beta"){
    # round down if parameter is on the edge
    parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]>=(1-1e-8)]=.99
    parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]<1e-8] =
      rep(0.001/(length(parm)-1),sum(parm[1:(length(parm)-1)]<1e-8))
    parm[length(parm)] = parm[length(parm)]*(parm[length(parm)]<1-1e-4) +
      1e-4*(parm[length(parm)]>=1-1e-4)
  } else{
    # round down if parameter is on the edge
    parm[parm>=(1-1e-8)]=.99
    parm[parm<1e-8] = rep(0.001/length(parm),sum(parm<1e-8))
  }


  jacobian = backward_transform_param_jacobian(logParm)

  if(sum(is.na(dataset))==0){
    gradient=GradLogLikParm_02(adj_positions=adj_positions,
                               parm=parm,
                               matList,
                               dataset,
                               interaction_effects=interaction_effects)
  } else{
    if(length(parm)>1){
      # use linearity of the determinant
      gradient =
        rowSums(
          sapply(
            1:nrow(dataset),
            function(t) GradLogLikParm_02(
              adj_positions=adj_positions,
              parm=parm,
              matList,
              t(dataset[t,]),
              interaction_effects=interaction_effects)))
    } else{
      # use linearity of the determinant
      gradient =
        sum(
          sapply(
            1:nrow(dataset),
            function(t) GradLogLikParm_02(
              adj_positions=adj_positions,
              parm=parm,
              matList,
              t(dataset[t,]),
              interaction_effects=interaction_effects)))
    }

  }
  return(gradient %*% jacobian)
}
