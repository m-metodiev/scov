#' Calculates the matrix D used in the quadratic minimization problem
#'
#' @inheritParams calc_Sigma_opt_frob
#' @keywords internal
#'
#' @returns the matrix D
calc_Dmat = function(matList_full){
  Dmat = matrix(nrow=length(matList_full),ncol=length(matList_full))
  for(i in seq_along(matList_full)){
    for(j in seq_along(matList_full)){
      Dmat[i,j] = frob_scalar_prod(matList_full[[i]],matList_full[[j]])
    }
  }
  return(2*Dmat)
}
