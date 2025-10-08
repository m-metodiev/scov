#' Calculates the vector d used in the quadratic minimization problem
#' @param S            the Pearson correlation matrix
#' @param matList_full combined list of pairwise and spatial matrices
#' @keywords internal
#'
#' @returns the vector d
calc_dvec = function(S, matList_full){
  dvec = numeric(length(matList_full))
  for(i in seq_along(matList_full)){
    dvec[i] = frob_scalar_prod(S, matList_full[[i]])
  }
  return(2*dvec)
}
