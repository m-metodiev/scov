#' Computes the matrices needed for the spatial effect
#'
#' @param adj_matrix the adjacency matrix of the spatial effect
#' @keywords internal
#' @returns a list containing Ml and Al, which are used to compute the
#'          correlation matrix for the spatial effect
get_M_A = function(adj_matrix){
  M = diag(rowSums(adj_matrix))
  A = adj_matrix
  A = (diag(1/rowSums(A))%*%A)
  Ml = list(M)
  Al = list(A)
  return(list(Ml=Ml,Al=Al))
}
