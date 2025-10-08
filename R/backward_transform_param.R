#' Calculates the backward transformation of the parameter
#'
#' @inheritParams avg_effect
#'
#' @returns the parameter after applying the backward transformation
#' @keywords internal
#' @importFrom ohenery  smax
#' @importFrom pracma   sigmoid
backward_transform_param = function(param){
  if(names(param)[length(param)]=="beta"){
    transformed_parm = param
    transformed_parm[1:(length(param)-1)] =
      ohenery::smax(c(-sum(transformed_parm[1:(length(param)-1)]),
                      transformed_parm[1:(length(param)-1)]))[-1]
    transformed_parm[length(param)] =
      pracma::sigmoid(transformed_parm[length(param)])
  } else{
    transformed_parm = ohenery::smax(c(-sum(param), param))[-1]
  }

  return(transformed_parm)
}
