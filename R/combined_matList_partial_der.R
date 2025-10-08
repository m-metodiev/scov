#' Computes the derivative of the correlation matrix w.r.t. beta
#'
#' @param matList                   a list of all matrices
#' @param link_matList              not in use (could be removed)
#' @param tilde_G_inv_partial_beta  derivative of G w.r.t. beta
#' @param interaction_effects       list of vector pairs indicating interactions
#' @keywords internal
#'
#' @returns the derivative of the correlation matrix w.r.t. beta
combined_matList_partial_der = function(matList, link_matList,
                                        tilde_G_inv_partial_beta,
                                        interaction_effects=list()){

  comb_mat = combined_matList(matList,interaction_effects=interaction_effects)
  matList_full =  comb_mat$matList_full
  dimension = dim(comb_mat$matList_full[[1]])[1]

  spatial_index = length(matList$Fk) + 1 # the index of the spatial effects
  indices = c(lapply(seq_along(matList$Fk),function(x) x),
              list(spatial_index),comb_mat$matrix_pairs)

  # use the product rule to compute the derivative
  for(s in seq_along(indices)){
    if(length(indices[[s]])==1){
      if(indices[[s]] == spatial_index){
        matList_full[[s]] = tilde_G_inv_partial_beta
      } else{
        matList_full[[s]] = matrix(0,dimension,dimension)
      }
    } else{
      if(indices[[s]][1] == spatial_index){
        matList_full[[s]] = matList$Fk[[indices[[s]][2]]] *
          tilde_G_inv_partial_beta
      } else if(indices[[s]][2] == spatial_index){
        matList_full[[s]] = matList$Fk[[indices[[s]][1]]] *
          tilde_G_inv_partial_beta
      } else{
        matList_full[[s]] = matrix(0,dimension,dimension)
      }
    }
  }
  return(matList_full)
}
