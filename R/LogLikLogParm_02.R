#' Computes (a translation of) the loglikelihood for the transformed parameter
#'
#' @param adj_positions       positions of spatial effect (if embedded)
#' @param logParm             the transformed parameter
#' @param matList             the list of matrices (pairwise + spatial)
#' @param dataset             an n (observations) x d (dimension) matrix
#' @param interaction_effects list of interaction effects (vectors of names)
#'
#' @returns (a translation of) the loglikelihood
#' @keywords internal
LogLikLogParm_02 <- function(adj_positions, logParm, matList, dataset,
                             interaction_effects=list()){
  parm = backward_transform_param(logParm)
  if(names(parm)[length(parm)]=="beta"){
    # "Round down" if parameters are too close to the edge
    parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]>=(1-1e-8)] = .99
    parm[1:(length(parm)-1)][parm[1:(length(parm)-1)]<1e-8] =
      rep(0.001/(length(parm)-1),sum(parm[1:(length(parm)-1)]<1e-8))
    if(sum(parm[1:(length(parm)-1)])>=(1-1e-8)){
      parm[1:(length(parm)-1)] =
        parm[1:(length(parm)-1)] / sum(parm[1:(length(parm)-1)]) * (1 - 1e-8)
    }
    parm[length(parm)] = parm[length(parm)] * (parm[length(parm)] < 1-1e-4) +
      1e-4*(parm[length(parm)]>=1-1e-4)
  } else{
    # "Round down" if parameters are too close to the edge
    parm[parm>=(1-1e-8)] = .99
    parm[parm<1e-8] = rep(0.001/length(parm),sum(parm<1e-8))
    if(sum(parm)>=(1-1e-8)){
      parm = parm / sum(parm) * (1 - 1e-8)
    }
  }

  LogLikParm_02(adj_positions=adj_positions,
                parm,
                matList,
                dataset,
                interaction_effects=interaction_effects)
}
