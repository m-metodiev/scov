#' Calculates the gradient of the function of the transformed parameter
#'
#' @param adj_positions       positions of the spatial effect (if embedded)
#' @param logParm             the transformed parameter
#' @param matList             the list of matrices (pairwise + spatial)
#' @param dataset             an n (observations) x d (dimensions) matrix
#' @param interaction_effects list of interaction effects (vectors of names)
#' @param joint_estimation    estimates everything jointly if TRUE,
#'                            uses a 2 step procedure if FALSE
#'
#' @returns the gradient of the function of the transformed parameter
#' @keywords internal
GradLogLikLogParm_02 <- function(adj_positions, logParm, matList, dataset,
                                 interaction_effects=list(),
                                 joint_estimation=FALSE){
  num_observations = nrow(dataset)
  if(joint_estimation){
    mus = logParm[1:ncol(dataset)]
    logsigmas = logParm[(ncol(dataset)+1):(2*ncol(dataset))]
    sigmas = exp(logsigmas)
    parm = backward_transform_param(logParm[-(1:(2*ncol(dataset)))])
    logParm = logParm[-(1:(2*ncol(dataset)))]
  } else{
    # dataset is already normalized by the first step, so no need for mu, sigma
    mus = NULL
    sigmas = NULL
    parm = backward_transform_param(logParm)
  }

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
  if(joint_estimation){
    jacobian = diag(length(mus) + length(sigmas) + length(logParm))

    # mu wasn't transformed, so the derivative stays 1
    # sigma was transformed, so the derivative is the derivative of the exponential
    diag(jacobian[(ncol(dataset)+1):(2*ncol(dataset)),
                  (ncol(dataset)+1):(2*ncol(dataset))]) =
      sigmas
    jacobian[-(1:(2*ncol(dataset))),-(1:(2*ncol(dataset)))] =
      backward_transform_param_jacobian(logParm)
  } else{
    jacobian = backward_transform_param_jacobian(logParm)
  }
  if(sum(is.na(dataset))==0){
    if(joint_estimation){
      normalized_dataset = (dataset - t(matrix(rep(mus,num_observations),
                                               ncol=num_observations)))%*%diag(1/sigmas)
    } else{
      normalized_dataset = dataset
    }
    gradient=GradLogLikParm_02(adj_positions=adj_positions,
                               parm=parm,
                               matList,
                               normalized_dataset,
                               interaction_effects=interaction_effects)
    if(joint_estimation){
      Sigma =
        as.matrix(CovMat_03(parm,
                            matList,adj_positions=adj_positions,
                            interaction_effects=interaction_effects)$Sigma)

      # due to numeric issues, Sigma could be computationally singular
      # derivative is set to be very low if that happens
      if(sum(is.na(eigen(Sigma)$values))==0){
        Omega = solve(Sigma)
      } else{
        Omega = -exp(16)*diag(num_observations)
      }
      gradient_mus = c(2*diag(1/(sigmas)) %*% Omega %*% diag(1/(sigmas)) %*%
        (colSums(dataset)/num_observations-mus))
      bigmatrix = diag(1/sigmas)%*%
        diag(-normalized_dataset[1,])%*%
        Omega %*%
        matrix(normalized_dataset[1,],ncol=1)
      for(t in 2:num_observations){
        bigmatrix = bigmatrix + diag(1/sigmas)%*%
          diag(-normalized_dataset[t,])%*%
          Omega %*%
          matrix(normalized_dataset[t,],ncol=1)
      }
      gradient_sigmas = - 2/sigmas - (2/num_observations)*bigmatrix
      gradient = c(gradient_mus,gradient_sigmas,gradient)
    }
  } else{
    if(joint_estimation){
      normalized_dataset = sapply(seq_along(dataset[1,]),
                                 function(s) (dataset[,s] - mus[s])/sigmas[s])
    } else{
      normalized_dataset = dataset
    }
    if(length(parm)>1){
      # use linearity of the determinant
      gradient =
        rowSums(
          sapply(
            1:num_observations,
            function(t) GradLogLikParm_02(
              adj_positions=adj_positions,
              parm=parm,
              matList,
              t(normalized_dataset[t,]),
              interaction_effects=interaction_effects)))
      if(joint_estimation){
        Sigma =
          as.matrix(CovMat_03(parm,
                              matList,adj_positions=adj_positions,
                              interaction_effects=interaction_effects)$Sigma)
        gradient_mus = rep(0,length(dataset[1,]))
        gradient_sigmas = rep(0,length(dataset[1,]))
        for(s in seq_along(dataset[,1])){
          gradient_mus[!is.na(dataset[s,])] = gradient_mus[!is.na(dataset[s,])] +
            c(2*diag(1/(sigmas[!is.na(dataset[s,])])) %*%
                solve(Sigma[!is.na(dataset[s,]),!is.na(dataset[s,])]) %*%
                diag(1/(sigmas[!is.na(dataset[s,])])) %*%
                matrix((dataset[s,][!is.na(dataset[s,])] - mus[!is.na(dataset[s,])]),
                       ncol=1))
          bigmatrix = diag(1/sigmas[!is.na(dataset[s,])])%*%
            diag(-normalized_dataset[s,][!is.na(dataset[s,])])%*%
            solve(Sigma[!is.na(dataset[s,]),!is.na(dataset[s,])]) %*%
            matrix(normalized_dataset[s,][!is.na(dataset[s,])],ncol=1)
          gradient_sigmas[!is.na(dataset[s,])] = gradient_sigmas[!is.na(dataset[s,])] +
            - 2/sigmas[!is.na(dataset[s,])] - 2*bigmatrix
        }
        # gradient_mus = sum(sapply(seq_along(dataset[,1]),
        #                       function(s) c(2*diag(1/(sigmas[!is.na(dataset[s,])])) %*%
        #                                       solve(Sigma[!is.na(dataset[s,]),!is.na(dataset[s,])]) %*%
        #                                       diag(1/(sigmas[!is.na(dataset[s,])])) %*%
        #                                       matrix((dataset[s,][!is.na(dataset[s,])] - mus[!is.na(dataset[s,])])/num_observations,
        #                                              ncol=1))))
        gradient = c(gradient_mus,gradient_sigmas,gradient)
      }
    } else{
      # use linearity of the determinant
      gradient =
        sum(
          sapply(
            1:num_observations,
            function(t) GradLogLikParm_02(
              adj_positions=adj_positions,
              parm=parm,
              matList,
              t(normalized_dataset[t,]),
              interaction_effects=interaction_effects)))
    }

  }
  return(gradient %*% jacobian)
}
