#' Computes the initial value estimator (IVE)
#'
#' This function computes the IVE estimator for large covariances in the
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
#' @importFrom stats          sd
#' @importFrom parallel       mclapply
#' @importFrom parallel       detectCores
#'
#' @returns Returns a named list with the following elements:
#'
#'          parm,             estimated parameters of pairwise, spatial effects
#'          average_effects,  average effects of the covariates
#'          corrmat_estim,    estimator of the correlation matrix
#'          covmat_estim,     estimator of the covariance matrix
#'
#'
#' @references Metodiev, M., Perrot-Dockès, M., Ouadah, S., Fosdick, B. K.,
#' Robin, S., Latouche, P., & Raftery, A. E. (2024). A Structured Estimator for
#' large Covariance Matrices in the Presence of Pairwise and Spatial Covariates.
#' arXiv preprint arXiv:2411.04520.
#'
#' @keywords internal
ive = function(pairwise_covariate_matrices=NULL, adj_matrix=NULL,
               dataset, mean_estim = NULL, sd_estim = NULL,
               grid_size = 100, parallelize=FALSE, ncores=8,
               adj_positions=1:nrow(adj_matrix), interaction_effects = list()){

  if(is.null(pairwise_covariate_matrices)&is.null(adj_matrix)){
    print("ERROR: (pairwise_covariate_matrices, adj_matrix) both NULL!")
    return(-1)
  }

  if(!is.null(pairwise_covariate_matrices)){
    if("spatial" %in% names(pairwise_covariate_matrices)){
      print("ERROR: Please use another name than - spatial - for the covariates.")
      return(-1)
    }

    if("beta" %in% names(pairwise_covariate_matrices)){
      print("ERROR: Please use another name than - beta - for the covariates.")
      return(-1)
    }

    # names the pairwise_covariate_matrices if they are not named
    if(is.null(names(pairwise_covariate_matrices))){
      names(pairwise_covariate_matrices) =
        LETTERS[1:length(pairwise_covariate_matrices)]
    }

    # if not positive semidefinite,
    # the matrices are mapped to positive semidefinite matrix
    pairwise_covariate_matrices=to_positive_definite(pairwise_covariate_matrices)
  }

  # uses standard mean and sd estimators if they are not already given
  if(sum(is.na(dataset))==0){
    if(is.null(mean_estim)){
      mean_estim = colMeans(dataset)
    }
    if(is.null(sd_estim)){
      sd_estim = apply(dataset,2,stats::sd)
    }
  } else{
    if(is.null(mean_estim)){
      mean_estim = apply(dataset,2,function(x) mean(x[!is.na(x)]))
    }
    if(is.null(sd_estim)){
      sd_estim = apply(dataset,2,function(x) stats::sd(x[!is.na(x)]))
    }
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
                                         beta=0.5)[adj_positions, adj_positions]),
                   Fk = pairwise_covariate_matrices,
                   vois = adj_matrix)
  } else{
    matList = list(Fk = pairwise_covariate_matrices)
  }

  comb_mat = combined_matList(matList,
                              interaction_effects=interaction_effects,
                              check_redundancy=TRUE)

  if(is.atomic(comb_mat)){
    return(-1)
  }

  correlation_matrix = correlation_matrix(dataset, mean_estim, sd_estim)
  dimension = dim(matList$Fk[[1]])[1]
  s = dim(matList$Al[[1]])[1]

  # grid-search not necessary if there is no spatial effect
  if(!is.null(adj_matrix)){
    # Need grid-search for beta because the norm is not quadratic w.r.t. beta
    beta=1.5
    xi = (1:(grid_size+1))/(grid_size+1)
    # tan-hyperbolic-spaced grid because beta approaches 1
    beta_vec = (1-tanh(beta*(1+xi))/tanh(beta))/
      (min((1-tanh(beta*(1+xi))/tanh(beta))))
    beta_vec = beta_vec[-length(beta_vec)]

    is_on_edge = TRUE # solution can lie on the edge of the parameter space
    edge_constraints = list()
    # params which lie on the edge will be adjusted in constraints

    counter = 0
    while(is_on_edge){
      counter = counter + 1
      grid_search = function(beta){
        matList$Gl[[1]]=
          tilde_G_inv(M=matList$Ml[[1]],A=matList$Al[[1]],
                      beta=beta, U_full=matList$U_full,
                      solve_U_full=matList$solve_U_full,
                      solve_M_no_islands=matList$solve_M_no_islands,
                      eigen_real=matList$eigen_real)[adj_positions,
                                                     adj_positions]
        matList_full = combined_matList(matList,
                                        interaction_effects)$matList_full
        res = calc_Sigma_opt_frob(matList_full, correlation_matrix,
                                  edge_constraints=edge_constraints)
        return(list(value=res$value,init=res$init))
      }

      test0=grid_search(beta_vec[1])
      test1=grid_search(beta_vec[grid_size])

      if(parallelize){

        # parallelize process
        cores=parallel::detectCores()
        this = simplify2array(
          parallel::mclapply((1: grid_size),
                             function(s) grid_search(beta_vec[s]),
                             mc.cores = min(cores[1]-1,ncores)))
      } else{
        this = sapply((1: grid_size), function(s) grid_search(beta_vec[s]))
      }
      res = sapply((1: grid_size), function(s) this[,s]$value)

      init <- c(this[,which.min(res)]$init,log(beta_vec[which.min(res)]))

      null_vec = which(round(c(1-sum(exp(init)[1:(length(exp(init))-1)]),
                               exp(init)[1:(length(exp(init))-1)]),15)<=0)
      if(length(null_vec)>0){
        matList$Gl[[1]] = tilde_G_inv(matList$Ml[[1]],matList$Al[[1]],
                                      exp(init)[length(init)])[adj_positions,
                                                               adj_positions]
        matList_full = combined_matList(matList,
                                        interaction_effects)$matList_full
        matList_full_extended = c(list(matrix(0,dimension,dimension)),
                                  matList_full)
        for(mat in matList_full_extended){
          diag(mat)=0
        }

        one_vec = (1:length(init))[-unique(c(1,null_vec))]
        # the null matrix can never be a target since its correlation is always 0

        ## choose vector pair with smallest distance in supports ##
        dist_matrix =
          sapply(null_vec, function(s) sapply(one_vec,
                                              function(t) mat_support_distance(
                                                matList_full_extended[[s]],
                                                matList_full_extended[[t]])))

        if(is.vector(dist_matrix)){
          # if there is only one option, choose the one
          min_vec = sapply(seq_along(null_vec),
                           function(s) which.min(dist_matrix[s]))
          arg_min1 = which.min(sapply(seq_along(null_vec),
                                      function(s) dist_matrix[s]))
          arg_min2 = 1
        } else{
          min_vec = sapply(seq_along(null_vec),
                           function(s) which.min(dist_matrix[,s]))
          arg_min1 = which.min(sapply(seq_along(null_vec),
                                      function(s) dist_matrix[min_vec[s],s]))
          arg_min2 = min_vec[arg_min1]
        }
        r_min = null_vec[arg_min1]-1
        # b chosen for the constraint of the form b>a/K
        s_min = one_vec[arg_min2]-1
        # a chosen for the constraint of the form b>a/K
        constraint_digit = ((exp(init)[-length(init)])[s_min])*
          mean(
            matList_full_extended[[s_min+
                                     1]][matList_full_extended[[s_min+1]]>0])/
          (length(matList_full)+1)
        # K chosen for the constraint of the form b>a/K

        # In the case that ALL of the values are 0
        if(length(one_vec)==0){
          s_min=0
          r_min=1
          constraint_digit = 1e-15
        }

        edge_constraints[[counter]] =
          list(r_min=r_min, s_min=s_min,
               constraint_digit=max(constraint_digit,1e-15))
        ## End: choose vector pair with smallest distance in supports ##
      } else{
        is_on_edge = FALSE
      }
    }
  } else{
    is_on_edge = TRUE # solution can lie on the edge of the parameter space
    edge_constraints = list()
    # params which lie on the edge will be adjusted in constraints

    counter = 0
    while(is_on_edge){
      counter = counter + 1
      matList_full = combined_matList(matList, interaction_effects)$matList_full

      res = calc_Sigma_opt_frob(matList_full, correlation_matrix,
                                edge_constraints=edge_constraints)



      init = res$init
      null_vec = which(round(c(1-sum(exp(init)[1:(length(exp(init))-1)]),
                               exp(init)[1:(length(exp(init))-1)]),15)<=0)
      if(length(null_vec)>0){
        matList_full = combined_matList(matList,
                                        interaction_effects)$matList_full
        matList_full_extended = c(list(matrix(0,dimension,dimension)),
                                  matList_full)
        for(mat in matList_full_extended){
          diag(mat)=0
        }

        one_vec = (1:length(init))[-unique(c(1,null_vec))]
        # the null matrix can never be a target since its correlation is 0

        ## choose vector pair with smallest distance in supports ##
        dist_matrix =
          sapply(null_vec,
                 function(s) sapply(one_vec,
                                    function(t) mat_support_distance(
                                      matList_full_extended[[s]],
                                      matList_full_extended[[t]])))

        if(is.vector(dist_matrix)){
          # if there is only one option, choose the one
          min_vec = sapply(seq_along(null_vec),
                           function(s) which.min(dist_matrix[s]))
          arg_min1 = which.min(sapply(seq_along(null_vec),
                                      function(s) dist_matrix[s]))
          arg_min2 = 1
        } else{
          min_vec = sapply(seq_along(null_vec),
                           function(s) which.min(dist_matrix[,s]))
          arg_min1 = which.min(sapply(seq_along(null_vec),
                                      function(s) dist_matrix[min_vec[s],s]))
          arg_min2 = min_vec[arg_min1]
        }
        r_min = null_vec[arg_min1]-1
        # b chosen for the constraint of the form b>a/K
        s_min = one_vec[arg_min2]-1
        # a chosen for the constraint of the form b>a/K
        constraint_digit = ((exp(init)[-length(init)])[s_min])*
          mean(
            matList_full_extended[[s_min+
                                     1]][matList_full_extended[[s_min+1]]>0])/
          (length(matList_full)+1)
        # K chosen for the constraint of the form b>a/K

        # In the case that ALL of the values are 0
        if(length(one_vec)==0){
          s_min=0
          r_min=1
          constraint_digit = 1e-15
        }

        edge_constraints[[counter]] =
          list(r_min=r_min, s_min=s_min,
               constraint_digit=max(constraint_digit,1e-15))
        ## End: choose vector pair with smallest distance in supports ##
      } else{
        is_on_edge = FALSE
      }
    }
  }

  init = exp(init)
  average_effects = avg_effect(parm=init,
                               matList=matList,
                               adj_positions=adj_positions,
                               interaction_effects=interaction_effects)
  effect_names=names(comb_mat$matList_full)
  names(average_effects) = effect_names
  if(!is.null(matList$Al)){
    names(init) = c(effect_names,"beta")
  } else{
    names(init) = c(effect_names)
  }

  corrmat_estim = CovMat_03(parm=init,
                            matList=matList,
                            adj_positions = adj_positions,
                            interaction_effects = interaction_effects)$Sigma
  return(list(parm = init,
              average_effects = average_effects,
              corrmat_estim = corrmat_estim,
              covmat_estim = diag(sd_estim) %*%
                corrmat_estim %*%
                diag(sd_estim)))
}
