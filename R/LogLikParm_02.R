#' Computes (a translation of) the loglikelihood
#'
#' @param adj_positions       positions of spatial effect (if embedded)
#' @param parm                the parameter
#' @param matList             the list of matrices (pairwise + spatial)
#' @param dataset             an n (observations) x d (dimension) matrix
#' @param interaction_effects list of interaction effects (vectors of names)
#' @param mus                 the mean parameters
#' @param sigmas              the standard deviation parameters
#'
#' @returns (a translation of) the loglikelihood
#' @keywords internal
LogLikParm_02 <- function(adj_positions, parm, matList, dataset,
                          interaction_effects=list(),
                          mus=NULL,
                          sigmas=NULL){
  num_observations = nrow(dataset)
  has_missing_values = sum(is.na(dataset))>0

  # the dataset is assumed to be standardized if mus and sigmas are missing
  if(is.null(mus)&is.null(sigmas)){
    if(has_missing_values){
      res = 0

      for(t in (1:num_observations)){
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
  } else{
    if(has_missing_values){
      res = 0
      for(t in (1:num_observations)){
        adj_pos_t = !is.na(dataset[t,])
        mus_t = mus[adj_pos_t]
        sigmas_t = sigmas[adj_pos_t]
        dataset_t = dataset[t,adj_pos_t]
        dataset_t_normalized = (dataset_t - mus_t) / sigmas_t
        Sigma =
          as.matrix(
            CovMat_03(adj_positions=adj_positions,
                      parm=parm,
                      matList=matList,
                      interaction_effects=interaction_effects)$Sigma)[adj_pos_t,
                                                                      adj_pos_t]
        S = dataset_t_normalized%*%t(dataset_t_normalized)
        res = res -
          sum(log(eigen(Sigma)$values)) -
          sum(diag(S%*%solve(Sigma))) -
          sum(2*log(sigmas_t))
      }
    } else{
      Sigma = as.matrix(CovMat_03(adj_positions=adj_positions,
                                  parm=parm, matList=matList,
                                  interaction_effects=interaction_effects)$Sigma)
      dataset_normalized = (dataset - t(matrix(rep(mus,num_observations),
                                               ncol=num_observations)))%*%diag(1/sigmas)
      S = dataset_normalized[1,]%*%t(dataset_normalized[1,])
      if(num_observations!=1){
        for(t in (2:num_observations)){
          S = S + dataset_normalized[t,]%*%t(dataset_normalized[t,])
        }
      }
      S = S/num_observations
      res = -sum(log(eigen(Sigma)$values))-
        sum(diag(S%*%solve(Sigma))) -
        sum(2*log(sigmas))
    }
  }

  return(res)
}
