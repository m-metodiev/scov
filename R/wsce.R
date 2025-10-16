#' Computes the weighted structured covariance matrix estimator (WSCE)
#'
#' This function computes the WSCE estimator for large covariances in the
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
#' @param sce_init                    the sce-initialization parameter vector
#' @param use_bootstrap               uses bootstrapping if TRUE
#' @param num_bootstrap_iters         number of bootstrap simulations
#' @param seed                        a seed
#' @param verbose                     prints progress if TRUE
#'
#' @returns Returns a named list with the following elements:
#'
#'          parm,             estimated parameters of pairwise, spatial effects,
#'          average_effects,  average effects of the covariates,
#'          corrmat_estim,    estimator of the correlation matrix,
#'          covmat_estim,     estimator of the covariance matrix,
#'          bic,              the Bayesian information criterion (BIC),
#'          lambda,           the asymptotically optimal weight of the WSCE
#'
#'
#' @references Metodiev, M., Perrot-Dockès, M., Ouadah, S., Fosdick, B. K.,
#' Robin, S., Latouche, P., & Raftery, A. E. (2024). A Structured Estimator for
#' large Covariance Matrices in the Presence of Pairwise and Spatial Covariates.
#' arXiv preprint arXiv:2411.04520.
#'
#' @keywords internal
#' @importFrom parallel       mclapply
#' @importFrom stats          cor
#' @importFrom stats          sd
#' @importFrom missMDA        imputePCA
#' @importFrom mvtnorm        rmvnorm
#' @importFrom withr          with_seed
#' @importFrom stats          cov
#' @importFrom stats          var
#' @importFrom purrr          quietly
wsce = function(pairwise_covariate_matrices, adj_matrix,
                dataset, mean_estim = NULL, sd_estim = NULL,
                grid_size=100, parallelize = FALSE, ncores=8,
                adj_positions=1:nrow(adj_matrix),
                interaction_effects=list(), init=NULL,
                sce_init=NULL, use_bootstrap=FALSE, num_bootstrap_iters=100,
                seed=0, verbose=TRUE){

  num_observations = nrow(dataset)
  if(!is.null(pairwise_covariate_matrices)){
    # if not positive semidefinite,
    # the matrices are mapped to positive semidefinite matrix
    pairwise_covariate_matrices=to_positive_definite(pairwise_covariate_matrices)
  }

  if(sum(is.na(dataset))==0){
    if(is.null(mean_estim)){
      use_bootstrap = TRUE # bootstrap is always used if the mean is estimated
      mean_estim = colMeans(dataset)
    }
    if(is.null(sd_estim)){
      use_bootstrap = TRUE # bootstrap is always used if the sd is estimated
      sd_estim = apply(dataset,2,stats::sd)
    }
    varepsilon = dataset
    varepsilon = (varepsilon - t(matrix(rep(mean_estim,num_observations),
                                        ncol=num_observations)))%*%
      diag(1/sd_estim)
    pearson_mat = cor_from_standard_errors(varepsilon)
  } else{
    if(is.null(mean_estim)){
      use_bootstrap = TRUE # bootstrap is always used if the mean is estimated
      mean_estim = apply(dataset,2,function(x) mean(x[!is.na(x)]))
    }
    if(is.null(sd_estim)){
      use_bootstrap = TRUE # bootstrap is always used if the sd is estimated
      sd_estim = apply(dataset,2,function(x) stats::sd(x[!is.na(x)]))
    }
    if(is.null(adj_matrix)){
      length_parm = length(pairwise_covariate_matrices)
    } else{
      # accounting for the two spatial variables
      length_parm = length(pairwise_covariate_matrices) + 2
    }

    varepsilon = missMDA::imputePCA(dataset,ncp=length_parm)$completeObs
    varepsilon = (varepsilon - t(matrix(rep(mean_estim,num_observations),
                                        ncol=num_observations)))%*%
      diag(1/sd_estim)
    pearson_mat = cor_from_standard_errors(varepsilon)
  }

  if(is.null(init)){
    ive_estim = ive(pairwise_covariate_matrices = pairwise_covariate_matrices,
                    adj_matrix = adj_matrix,
                    dataset = dataset,
                    mean_estim = mean_estim, sd_estim = sd_estim,
                    grid_size = grid_size, ncores = ncores,
                    adj_positions=adj_positions,
                    interaction_effects=interaction_effects,
                    parallelize = parallelize)
    if(is.atomic(ive_estim)){
      return(-1)
    }
    init = ive_estim$parm
  }

  if(is.null(sce_init)){
    sce_estim = sce(pairwise_covariate_matrices = pairwise_covariate_matrices,
                    adj_matrix = adj_matrix,
                    dataset = dataset,
                    mean_estim = mean_estim, sd_estim = sd_estim,
                    grid_size = grid_size, ncores = ncores,
                    adj_positions = adj_positions,
                    init = init,
                    interaction_effects = interaction_effects,
                    parallelize = parallelize,
                    verbose=verbose)
    sce_init = sce_estim$parm
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

  SCE_mat = CovMat_03(adj_positions=adj_positions,
                      parm=sce_init,
                      matList=matList,
                      interaction_effects=interaction_effects)$Sigma
  parm=sce_init

  if(!is.null(adj_matrix)){
    # numerical issues
    if(parm[length(parm)]<1e-1){
      parm[length(parm)]=1e-1
    }
    if(parm[length(parm)-1]<1e-4){
      parm[length(parm)-1]=1e-4
    }
  }

  if(!is.null(matList$Fk)){
    dimension = dim(matList$Fk[[1]])[1]
  } else{
    dimension = length(adj_positions)
  }

  # more expensive, but gives more accurate estimates
  try_again = TRUE
  while(try_again){
    if(sum(is.na(dataset))>0 | use_bootstrap){
      test_func = function(s){

        perc = round((s / num_bootstrap_iters) * 100)

        # make sure it is divisable by 10
        if(verbose){
          num_rounded = num_bootstrap_iters - (num_bootstrap_iters %% 10)
          if(num_rounded == 0){
            system(sprintf('echo "%s"', paste0(perc," percent", collapse="")))
          } else{
            if(s==1){
              system(sprintf('echo "%s"', paste0(perc," percent", collapse="")))
            }
            if(((s / num_bootstrap_iters) * 10) ==
               round((s / num_bootstrap_iters) * 10)){
              system(sprintf('echo "%s"', paste0(perc," percent", collapse="")))
            }
          }
        }

        # ensures replicability
        test_data = withr::with_seed(seed=s+seed,
                                     mvtnorm::rmvnorm(num_observations,
                                                      sigma=pearson_mat))
        # ensures replicability
        # test_data = mvtnorm::rmvnorm(num_observations, sigma=pearson_mat)
        test_data[is.na(dataset)]=NA

        # impute data if values are missing, use Pearson matrix otherwise
        if(sum(is.na(dataset))>0){
          sim_test = list(dataset=test_data,
                          correlation_matrix=compute_marginal_cor(test_data))
          pearson_test = cor_from_standard_errors(
            missMDA::imputePCA(test_data,ncp=length(parm))$completeObs)
        } else{
          sim_test =
            list(dataset=test_data,
                 correlation_matrix=cor_from_standard_errors(test_data))
          pearson_test = stats::cor(test_data)
        }

        compute_noisy =
          function() sce(pairwise_covariate_matrices, adj_matrix,
                         test_data, mean_estim = rep(0,dimension),
                         sd_estim = rep(1,dimension), grid_size=grid_size,
                         adj_positions=adj_positions,
                         interaction_effects=interaction_effects,
                         init=parm)$corrmat_estim
        compute_quietly = purrr::quietly(compute_noisy)
        SCE_test = compute_quietly()$result

        return(array(c(SCE_test, pearson_test), dim = c(dimension,
                                                        dimension,
                                                        2)))
      }

      pearson_mat=as.matrix(pearson_mat)
      if(verbose){
        message("loading ...", quote=FALSE)
      }
      if(parallelize){
        cores=detectCores()
        bootstrap_sample = parallel::mclapply(1:num_bootstrap_iters,
                                              FUN=function(s) test_func(s),
                                              mc.cores=min(cores[1]-1,ncores))
      } else{
        bootstrap_sample = lapply(1:num_bootstrap_iters,
                                  FUN=function(s) test_func(s))
      }
      bootstrap_sample = simplify2array(bootstrap_sample)

      test_var = sum(apply(bootstrap_sample,1,
                           function(bsrow) apply(bsrow,1,
                                                 function(bscol)
                                                   stats::var(bscol[2,])*
                                                   num_observations)))
      test_cov = sum(apply(bootstrap_sample,1,
                           function(bsrow) apply(bsrow,1,
                                                 function(bscol) stats::cov(
                                                   bscol[1,],
                                                   bscol[2,])*
                                                   num_observations)))
      test_mse=
        sum((apply(bootstrap_sample,1,
                   function(bsrow) apply(bsrow,1,
                                         function(bscol) mean(bscol[1,])))-
                      pearson_mat)^2)

      approx_pi =   test_mse
      approx_mse = sum((SCE_mat- pearson_mat )^2)
      (lambda = (test_var-test_cov)/(num_observations*test_mse))
      try_again = FALSE
    } else{
      # the covariance matrix is given by the inverse of the Fisher information
      Fisher_mat = Fisher_information(adj_positions, parm, matList,
                                      #link_der_beta,
                                      interaction_effects=interaction_effects)
      Fisher_mat = try(chol2inv(chol(Fisher_mat)), silent = TRUE)
      if (inherits(Fisher_mat, "try-error")) {
        use_bootstrap = TRUE
        try_again = TRUE
      } else{
        Sigma_der = GradLogLikParm_02(adj_positions, parm, matList,
                                      dataset=matrix(0,ncol=dimension,
                                                     nrow=dimension),
                                      interaction_effects=interaction_effects,
                                      #link_der_beta=link_der_beta,
                                      return_Sigma_der=TRUE)
        dimension_small = dimension
        total_der = matrix(0, ncol=dimension_small*dimension_small,
                           nrow=length(parm))
        for(i in seq_along(parm)){
          total_der[i,] = c(as.matrix(Sigma_der[[i]]))
        }

        #compute approximate covariance matrix of the SCE by the delta method
        approx_var_SCE = sum(diag((t(total_der)%*%Fisher_mat%*%total_der)))
        approx_pi = sum((1-pearson_mat^2)^2)
        approx_mse = sum((SCE_mat- pearson_mat )^2)
        approx_var_SCE = sum(diag(t(total_der)%*%Fisher_mat%*%total_der))
        (lambda = (sqrt(approx_pi)*(sqrt(approx_pi)-sqrt(approx_var_SCE)))/
            (num_observations*approx_mse))
        try_again = FALSE
      }
    }
  }

  lambda = min(lambda,1)
  lambda = max(lambda,0)
  #lambda has to be between 0 and 1

  bic = -2*true_LogLikParm_02(adj_positions, sce_init, matList, dataset,
                              interaction_effects=interaction_effects) +
    length(init)*log(dim(dataset)[1])

  average_effects = avg_effect(parm=sce_init,
                               matList=matList,
                               adj_positions=adj_positions,
                               interaction_effects=interaction_effects)
  comb_mat = combined_matList(matList,interaction_effects=interaction_effects)
  effect_names = names(comb_mat$matList_full)
  names(average_effects) = effect_names
  if(!is.null(adj_matrix)){
    names(sce_init) = c(effect_names,"beta")
  } else{
    names(sce_init) = c(effect_names)
  }

  corrmat_estim = lambda*as.matrix(SCE_mat) +
    (1-lambda)*as.matrix(pearson_mat)

  return(list(parm=sce_init,
              average_effects = average_effects,
              corrmat_estim = corrmat_estim,
              covmat_estim = diag(sd_estim)%*%corrmat_estim%*%diag(sd_estim),
              bic=bic,
              lambda=lambda))
}
