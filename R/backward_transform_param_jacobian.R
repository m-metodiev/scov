#' Calculates the Jacobian of the backwards transformation
#'
#' @inheritParams avg_effect
#'
#' @returns the Jacobian of the backwards transformation
#' @keywords internal
#' @importFrom pracma          sigmoid
backward_transform_param_jacobian = function(param){
  # Jacobian of the softmax
  soft_max_der = function(par) diag(par)-
    matrix(rep(par,length(par)),
           ncol=length(par),
           nrow=length(par))*t(matrix(rep(par,length(par)),
                                      ncol=length(par),nrow=length(par)))

  # multivariate chain rule component
  sum_der = function(par) rbind(-1,diag(length(par)))

  if(names(param)[length(param)]=="beta"){
    jacobian_parm_full = param
    jacobian_parm = param[1:(length(param)-1)]
    jacobian_parm_matrix = cbind(0,diag(length(jacobian_parm)))%*%
      soft_max_der(smax(c(-sum(jacobian_parm),jacobian_parm)))%*%
      sum_der(jacobian_parm)

    jacobian_parm_full[length(param)] =
      pracma::sigmoid(jacobian_parm_full[length(param)])*
      (1-pracma::sigmoid(jacobian_parm_full[length(param)]))
    # derivative of the sigmoid function

    jacobian_parm_matrix = rbind(cbind(jacobian_parm_matrix,
                                       rep(0,length(param)-1)),
                                 c(rep(0,length(param)-1),
                                   jacobian_parm_full[length(param)]))
  } else{
    jacobian_parm = param
    jacobian_parm_matrix = cbind(0,diag(length(jacobian_parm)))%*%
      soft_max_der(smax(c(-sum(jacobian_parm),jacobian_parm)))%*%
      sum_der(jacobian_parm)
  }
  return(jacobian_parm_matrix)
}
