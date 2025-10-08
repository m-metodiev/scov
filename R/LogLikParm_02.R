#' Computes (a translation of) the loglikelihood
#'
#' @param adj_positions       positions of spatial effect (if embedded)
#' @param parm                the parameter
#' @param matList             the list of matrices (pairwise + spatial)
#' @param dataset             an n (observations) x d (dimension) matrix
#' @param interaction_effects list of interaction effects (vectors of names)
#'
#' @returns (a translation of) the loglikelihood
#' @keywords internal
LogLikParm_02 <- function(adj_positions, parm, matList, dataset,
                          interaction_effects=list()){
  num_observations = nrow(dataset)
  has_missing_values = sum(is.na(dataset))>0
  if(has_missing_values){
    res = 0

    for(t in (1:num_observations)){
      #browser()
      adj_pos_t = !is.na(dataset[t,])
      Sigma =
        as.matrix(
          CovMat_03(adj_positions=adj_positions,
                    parm=parm,
                    matList=matList,
                    interaction_effects=interaction_effects)$Sigma)[adj_pos_t,
                                                                    adj_pos_t]
      S = dataset[t,adj_pos_t]%*%t(dataset[t,adj_pos_t])
      res = res - sum(log(eigen(Sigma)$values)) - sum(diag(S%*%solve(Sigma)))
    }

  } else{
    Sigma = as.matrix(CovMat_03(adj_positions=adj_positions,
                                parm=parm, matList=matList,
                                interaction_effects=interaction_effects)$Sigma)
    S = dataset[1,]%*%t(dataset[1,])
    if(num_observations!=1){
      for(t in (2:num_observations)){
        S = S + dataset[t,]%*%t(dataset[t,])
      }
    }
    S = S/num_observations
    res = -sum(log(eigen(Sigma)$values))-sum(diag(S%*%solve(Sigma)))
  }
  return(res)
}
