#' Computes a measure of distance between the support of two matrices
#'
#' @param mat1 a square matrix
#' @param mat2 a square matrix
#'
#' @returns a measure of distance between the support of two matrices
#' @importFrom pracma          ceil
#' @keywords internal
mat_support_distance = function(mat1, mat2){
  return(mean(abs(pracma::ceil(mat1)-ceil(mat2))))
}
