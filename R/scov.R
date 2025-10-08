#' Computes a structured estimator for covariance matrices
#'
#' This function computes the WSCE, SCE or IVE estimator for large covariances
#' in the presence of pairwise and spatial covariates from
#' Metodiev et al. (2024).
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
#' @param num_bootstrap_iters         number of bootstrap simulations
#' @param use_bootstrap               uses bootstrapping if TRUE
#' @param semiparametric              computes the IVE if TRUE, the SCE else
#' @param misspecification            computes the WSCE if TRUE, the WSCE else
#' @param seed                        a seed (can't be set to NULL)
#'
#' @returns Returns a named list with the following elements:
#'
#'          parm,             estimated parameters of pairwise, spatial effects
#'          average_effects,  average effects of the covariates
#'          corrmat_estim,    estimator of the correlation matrix
#'          covmat_estim,     estimator of the covariance matrix
#'          bic,              the Bayesian information criterion (BIC)
#'          lambda,           the asymptotically optimal weight of the WSCE
#'
#'
#' @references Metodiev, M., Perrot-Dockès, M., Ouadah, S., Fosdick, B. K.,
#' Robin, S., Latouche, P., & Raftery, A. E. (2024). A Structured Estimator for
#' large Covariance Matrices in the Presence of Pairwise and Spatial Covariates.
#' arXiv preprint arXiv:2411.04520.
#'
#' @export
#'
#' @examples
#'
#' intercept = matrix(1,ncol=4,nrow=4)
#' X1 = rbind(c(1,1,1,0),c(1,1,1,0),c(1,1,1,0),c(0,0,0,1))
#' X2 = rbind(c(1,0,0,0),c(0,1,1,1),c(0,1,1,1),c(0,1,1,1))
#' covar_mats = list(intercept=intercept,X1=X1,X2=X2)
#' adj_matrix = rbind(c(0,1,0,0),c(1,0,0,0),c(0,0,0,1),c(0,0,1,0))
#' mean = rep(0,4)
#' sigma = 0.05*intercept+0.2*X1+0.2*X2+0.1*X2*X1+0.4*(diag(4) + adj_matrix)
#' diag(sigma) = 1
#' dataset = mvtnorm::rmvnorm(1000,mean=mean,sigma=sigma)
#' scov(covar_mats, adj_matrix, dataset, interaction_effects=list(c("X1","X2")))
scov = function(pairwise_covariate_matrices, adj_matrix,
                dataset, mean_estim = NULL, sd_estim = NULL,
                grid_size=100, parallelize = FALSE, ncores=8,
                adj_positions=1:nrow(adj_matrix),
                interaction_effects=list(), init=NULL,
                use_bootstrap=FALSE, num_bootstrap_iters=100,
                semiparametric=FALSE, misspecification=FALSE, seed=0){


  if(!is.null(pairwise_covariate_matrices)){
    # transform to matrices if necessary
    pairwise_covariate_matrices_copy = list()
    for(i in seq_along(pairwise_covariate_matrices)){
      pairwise_covariate_matrices_copy[[i]] =
        as.matrix(pairwise_covariate_matrices[[i]])
    }
    names(pairwise_covariate_matrices_copy) = names(pairwise_covariate_matrices)
    pairwise_covariate_matrices = pairwise_covariate_matrices_copy

    # if not positive semidefinite,
    # the matrices are mapped to positive semidefinite matrix
    pairwise_covariate_matrices=
      to_positive_definite(pairwise_covariate_matrices)
  }

  if(!is.null(adj_matrix)){
    adj_matrix = as.matrix(adj_matrix)
  }

  if(semiparametric & (!misspecification)){

    ive_estim = ive(pairwise_covariate_matrices = pairwise_covariate_matrices,
                    adj_matrix = adj_matrix,
                    dataset = dataset,
                    mean_estim = mean_estim, sd_estim = sd_estim,
                    grid_size = grid_size, ncores = ncores,
                    adj_positions=adj_positions,
                    interaction_effects=interaction_effects,
                    parallelize = parallelize)
    return(ive_estim)
  } else{
    if(misspecification){
      wsce_estim = wsce(pairwise_covariate_matrices=pairwise_covariate_matrices,
                        adj_matrix=adj_matrix,
                        dataset=dataset,
                        mean_estim = mean_estim, sd_estim = sd_estim,
                        grid_size=grid_size, parallelize = parallelize,
                        ncores=ncores, adj_positions = adj_positions,
                        interaction_effects = interaction_effects,
                        init = init,
                        use_bootstrap=use_bootstrap,
                        num_bootstrap_iters=num_bootstrap_iters,
                        seed=seed)
      return(wsce_estim)
    } else{
      sce_estim = sce(pairwise_covariate_matrices=pairwise_covariate_matrices,
                      adj_matrix = adj_matrix,
                      dataset = dataset,
                      mean_estim = mean_estim, sd_estim = sd_estim,
                      grid_size = grid_size, ncores = ncores,
                      adj_positions = adj_positions,
                      init = init,
                      interaction_effects=interaction_effects,
                      parallelize = parallelize)
      return(sce_estim)
    }
  }
}
