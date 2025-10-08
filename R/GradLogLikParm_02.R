#' Calculates the gradient of the loglikelihood or the gradient of Sigma
#'
#' @param adj_positions         positions of the spatial effect (if embedded)
#' @param parm                  the parameter
#' @param matList               the list of matrices (pairwise + spatial)
#' @param dataset               an n (observations) x d (dimensions) matrix
#' @param interaction_effects   list of interaction effects (vectors of names)
#' @param return_Sigma_der      returns the gradient of Sigma if TRUE
#'
#' @returns the gradient of the loglikelihood or the gradient of Sigma
#' @keywords internal
GradLogLikParm_02 <- function(adj_positions, parm, matList, dataset,
                              interaction_effects=list(),
                              return_Sigma_der=FALSE){
  # These are actually 2 completely different functions:
  # One returns the gradient of the loglikelihood, one the gradient of Sigma
  # You can choose which one to use by setting the parameter return_Sigma_der
  # Warning: If the dataset contains missing values,
  # it is expected to be a matrix with only 1 row
  if(!is.null(matList$Fk)){
    dimension = dim(matList$Fk[[1]])[1]
  } else{
    dimension = length(adj_positions)
  }
  num_observations = nrow(dataset)
  if(num_observations==1){
    covY <- t(t(dataset[,!is.na(dataset[1,])])) %*%
      t(dataset[,!is.na(dataset[1,])]) / num_observations
  } else{
    covY <- t(dataset) %*% dataset / num_observations
  }

  if(names(parm)[length(parm)]=="beta"){
    # Calculate derivatives for matrix from the CAR model
    l=1
    betal = parm[length(parm)]
    G_inv_list = calc_tilde_G_inv_partial_beta(matList$Ml[[l]],
                                               matList$Al[[l]],
                                               betal[l])
    G_inv_list$tilde_G_inv = G_inv_list$tilde_G_inv[adj_positions,adj_positions]
    G_inv_list$tilde_G_inv_partial_beta =
      G_inv_list$tilde_G_inv_partial_beta[adj_positions,adj_positions]
    matList$Gl[[1]] = G_inv_list$tilde_G_inv
    gradLogLik.beta <- numeric(length(matList$Gl))
  }

  link_matList = combined_matList(matList,interaction_effects)$matList_full
  # all operations are performed on this matList

  gradLogLik.alpha_delta <- numeric(length(link_matList))

  # probably only works if numbers of Gl-matrices is equal to 1

  Sigma =
    as.matrix(CovMat_03(parm,
                        matList,adj_positions=adj_positions,
                        interaction_effects=interaction_effects)$Sigma)[!is.na(
                          dataset[1,]),!is.na(dataset[1,])]

  # due to numeric issues, Sigma could be computationally singular
  # derivative is set to be very low if that happens
  if(sum(is.na(eigen(Sigma)$values))==0){
    Omega = solve(Sigma)
  } else{
    Omega = -exp(16)*diag(num_observations)
  }

  for(k in seq_along(link_matList)){
    diag(link_matList[[k]]) = 0 # there is an identity matrix as the first F
    gradLogLik.alpha_delta[k] =
      - sum(diag((link_matList[[k]][!is.na(dataset[1,]),
                                    !is.na(dataset[1,])]) %*% Omega)) +
      sum(diag(covY %*%
                 Omega %*%
                 (link_matList[[k]][!is.na(dataset[1,]),
                                    !is.na(dataset[1,])]) %*%
                 Omega))
  }

  if(names(parm)[length(parm)]=="beta"){
    tilde_G_inv_partial_beta = G_inv_list$tilde_G_inv_partial_beta
    beta_der_list =
      combined_matList_partial_der(matList,
                                   link_matList,
                                   tilde_G_inv_partial_beta,
                                   interaction_effects=interaction_effects)

    Sigma_partial_beta = matrix(0,dimension,dimension)[!is.na(dataset[1,]),
                                                       !is.na(dataset[1,])]
    for(k in seq_along(beta_der_list)){
      Sigma_partial_beta = Sigma_partial_beta +
        parm[k] * beta_der_list[[k]][!is.na(dataset[1,]),!is.na(dataset[1,])]
    }
    gradLogLik.beta[l] = (-sum(diag(Sigma_partial_beta %*% Omega)) +
                            sum(diag(covY %*%
                                       Omega %*%
                                       Sigma_partial_beta %*%
                                       Omega)))
  } else{
    Sigma_partial_beta = NULL
    gradLogLik.beta = NULL
    # not defined if no spatial effect
  }

  if(return_Sigma_der){
    return(c(link_matList=link_matList,
             list(Sigma_partial_beta=Sigma_partial_beta)))
  } else{
    return(c(gradLogLik.alpha_delta, gradLogLik.beta))
  }
}
