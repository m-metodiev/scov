#' Computes the structured covariance matrix estimator (SCE)
#'
#' This function computes the SCE estimator for large covariances in the
#' presence of pairwise and spatial covariates from Metodiev et al. (2024).
#'
#' @param pairwise_covariate_matrices named list of square matrices
#' @param adj_matrix                  adjacency matrix of the spatial covariate
#' @param dataset                     the dataset given in matrix form
#' @param mean_estim                  mean vector estimate
#' @param sd_estim                    standard deviation vector estimate
#' @param grid_size                   grid-size for spatial effect
#' @param parallelize                 uses parallel-processing if TRUE
#' @param ncores                      number of cores for the parallelization
#' @param adj_positions               positions within the adjacency matrix
#' @param interaction_effects         list of interaction effects
#' @param init                        the initialization parameter vector
#' @param verbose                     prints progress if TRUE
#' @keywords internal
#' @importFrom stats          sd
#' @importFrom stats          optim
#'
#' @returns Returns a named list with the following elements:
#'
#'          parm,             estimated parameters of pairwise, spatial effects,
#'          average_effects,  average effects of the covariates,
#'          corrmat_estim,    estimator of the correlation matrix,
#'          covmat_estim,     estimator of the covariance matrix,
#'          bic,              the Bayesian information criterion (BIC)
#'
#'
#' @references Metodiev, M., Perrot-Dockès, M., Ouadah, S., Fosdick, B. K.,
#' Robin, S., Latouche, P., & Raftery, A. E. (2024). A Structured Estimator for
#' large Covariance Matrices in the Presence of Pairwise and Spatial Covariates.
#' arXiv preprint arXiv:2411.04520.
sce = function(pairwise_covariate_matrices, adj_matrix,
               dataset, mean_estim = NULL, sd_estim = NULL,
               grid_size=100, parallelize = TRUE, ncores=8,
               adj_positions=1:nrow(adj_matrix), interaction_effects=list(),
               init=NULL,verbose=TRUE){

  if(is.null(init)){
    ive_estim = ive(pairwise_covariate_matrices = pairwise_covariate_matrices,
                    adj_matrix = adj_matrix,
                    dataset = dataset,
                    mean_estim = mean_estim, sd_estim = sd_estim,
                    grid_size = grid_size, ncores = ncores,
                    adj_positions=adj_positions,
                    interaction_effects=interaction_effects)
    if(is.atomic(ive_estim)){
      return(-1)
    }
    init = ive_estim$parm
  }

  init = forward_transform_param(init)

  num_observations = nrow(dataset)
  # the SCE needs to be computed on the normalized dataset
  if(sum(is.na(dataset))==0){
    if(is.null(mean_estim)){
      mean_estim = colMeans(dataset)
    }
    if(is.null(sd_estim)){
      sd_estim = apply(dataset,2, stats::sd)
    }
    dataset = (dataset - t(matrix(rep(mean_estim,num_observations),
                                  ncol=num_observations)))%*%diag(1/sd_estim)
  } else{
    if(is.null(mean_estim)){
      mean_estim = apply(dataset,2,function(x) mean(x[!is.na(x)]))
    }
    if(is.null(sd_estim)){
      sd_estim = apply(dataset,2,function(x) stats::sd(x[!is.na(x)]))
    }
    dataset = sapply(seq_along(dataset[1,]),
                     function(s) (dataset[,s] - mean_estim[s])/sd_estim[s])
  }

  if((!is.null(pairwise_covariate_matrices))&
     is.null(names(pairwise_covariate_matrices))){
    names(pairwise_covariate_matrices) =
      LETTERS[1:length(pairwise_covariate_matrices)]
  }

  # define a list of relevant matrices:
  #   Ml: used to define the correlation matrix of the CAR model
  #   Al: used to define the correlation matrix of the CAR model
  #   Gl: the correlation matrix of the CAR model
  #   Fk: pairwise correlation matrices
  if(!is.null(adj_matrix)){
    Ml_Al = get_M_A(adj_matrix=adj_matrix)
    matList = list(Ml = Ml_Al$Ml,
                   Al = Ml_Al$Al,
                   Gl = list(tilde_G_inv(M=Ml_Al$Ml[[1]],
                                         A=Ml_Al$Al[[1]],
                                         beta=0.5)[adj_positions,
                                                   adj_positions]),
                   Fk = pairwise_covariate_matrices,
                   vois = adj_matrix)
  } else{
    matList = list(Fk = pairwise_covariate_matrices)
  }

  LogLikLogParm = function(x) LogLikLogParm_02(
    adj_positions=adj_positions,
    logParm=x,
    matList=matList,
    dataset=dataset,
    interaction_effects=interaction_effects)
  GradLogLikLogParm = function(x) GradLogLikLogParm_02(
    adj_positions=adj_positions,
    logParm=x,
    matList=matList,
    dataset=dataset,
    interaction_effects=interaction_effects)

  logLikInit <- LogLikLogParm(init)
  if(verbose){
    fit3 <- try(stats::optim(par=init,
                             fn=LogLikLogParm,
                             gr=GradLogLikLogParm,
                             control=list(fnscale=-1,
                                          trace=1,
                                          maxit=500),
                             method='BFGS'))
  } else{
    fit3 <- try(stats::optim(par=init,
                             fn=LogLikLogParm,
                             gr=GradLogLikLogParm,
                             control=list(fnscale=-1,
                                          trace=0,
                                          maxit=500),
                             method='BFGS'))
  }

  if(!is.character(fit3[1])){
    SigmaHat3 <- CovMat_03(adj_positions=adj_positions,
                           parm=backward_transform_param(fit3$par),
                           matList=matList,
                           interaction_effects=interaction_effects)$Sigma
    param_fit3 = backward_transform_param(fit3$par)
  } else{
    SigmaHat3 = NULL
    param_fit3 = NULL
  }

  bic = -2*true_LogLikParm_02(adj_positions, param_fit3, matList, dataset,
                              interaction_effects=interaction_effects) +
    length(init)*log(dim(dataset)[1])

  average_effects = avg_effect(parm=param_fit3,
                               matList=matList,
                               adj_positions=adj_positions,
                               interaction_effects=interaction_effects)
  comb_mat = combined_matList(matList,interaction_effects=interaction_effects)
  effect_names = names(comb_mat$matList_full)
  names(average_effects) = effect_names
  if(!is.null(adj_matrix)){
    names(param_fit3) = c(effect_names,"beta")
  } else{
    names(param_fit3) = effect_names
  }

  corrmat_estim = SigmaHat3
  return(list(parm=param_fit3,
              average_effects = average_effects,
              corrmat_estim = corrmat_estim,
              covmat_estim = diag(sd_estim)%*%corrmat_estim%*%diag(sd_estim),
              bic=bic))
}
