#' Calculates the Frobenius inner product between to square matrices
#'
#' @param A a square matrix
#' @param B a square matrix
#'
#' @returns the Frobenius inner product
#' @keywords internal
frob_scalar_prod = function(A,B){
  sum(diag(t(A)%*%B))
}
