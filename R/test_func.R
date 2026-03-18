#' Bootstrap simulations for the WSCE
#'
#' This function is computed multiple times as a bootstrap sample
#'
#' @param s                             label of the sample
#' @param num_bootstrap_iters           number of iterations
#' @param seed                          the seed
#' @param dimension                     the dimension of the data
#' @param num_observations              size of the data
#' @param pearson_mat                   the pearson correlation matrix
#' @param dataset                       the dataset
#' @param parm                          the parameter estimated from the sce
#' @param pairwise_covariate_matrices   the list of covariate matrices
#' @param adj_matrix                    adjacency matrix of the spatial effect
#' @param grid_size                     size of the grid for the spatial effect
#' @param adj_positions                 positions within the adjacency matrix
#' @param interaction_effects           list of interaction effects
#' @param verbose                       prints progress if TRUE
#' @param joint_estimation              estimates everything jointly if TRUE,
#'                                      uses a 2 step procedure if FALSE
#'
#' @returns Returns an array containing the simulated value of the SCE and
#'          the Pearson correlation matrix
#'
#' @keywords internal
#' @importFrom stats          cor
#' @importFrom missMDA        imputePCA
#' @importFrom mvtnorm        rmvnorm
#' @importFrom withr          with_seed
#' @importFrom purrr          quietly
test_func = function(s,
                     num_bootstrap_iters,
                     seed,
                     dimension,
                     num_observations,
                     pearson_mat,
                     dataset,
                     parm,
                     pairwise_covariate_matrices,
                     adj_matrix,
                     grid_size,
                     adj_positions,
                     interaction_effects,
                     verbose,
                     joint_estimation){
  perc = round((s / num_bootstrap_iters) * 100)

  # make sure it is divisable by 10
  if(verbose){
    num_rounded = num_bootstrap_iters - (num_bootstrap_iters %% 10)
    if(num_rounded == 0){
      system(sprintf('echo "%s"', paste0("+10 percent", collapse="")))
    } else{
      if(s==1){
        system(sprintf('echo "%s"', paste0("+10 percent", collapse="")))
      }
      if(((s / num_bootstrap_iters) * 10) ==
         round((s / num_bootstrap_iters) * 10)){
        system(sprintf('echo "%s"', paste0("+10 percent", collapse="")))
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
                   joint_estimation=joint_estimation,
                   init=parm)$corrmat_estim
  compute_quietly = purrr::quietly(compute_noisy)
  SCE_test = compute_quietly()$result

  return(array(c(SCE_test, pearson_test), dim = c(dimension,
                                                  dimension,
                                                  2)))
}
