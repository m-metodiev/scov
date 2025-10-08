#' Transforms the parameter using a logit and inverse softmax
#'
#' @param param the parameter
#'
#' @returns the transformed parameter
#' @keywords internal
#' @importFrom pracma           logit
#' @importFrom ohenery          inv_smax
forward_transform_param = function(param){

  if(names(param)[length(param)]=="beta"){
    # "Round down" if parameters are too close to the edge
    param[1:(length(param)-1)][param[1:(length(param)-1)]>=(1-1e-8)]=.99
    param[1:(length(param)-1)][param[1:(length(param)-1)]<1e-8]=
      rep(0.001/(length(param)-1),sum(param[1:(length(param)-1)]<1e-8))
    if(sum(param[1:(length(param)-1)])>=(1-1e-8)){
      param[1:(length(param)-1)] =
        param[1:(length(param)-1)]/sum(param[1:(length(param)-1)])*(1-1e-8)
    }
    param[length(param)] =
      param[length(param)]*(param[length(param)]<1-1e-4) +
      1e-4*(param[length(param)]>=1-1e-4)

    transformed_init = param
    transformed_init[1:(length(param)-1)] =
      ohenery::inv_smax(c(1-sum(transformed_init[1:(length(param)-1)]),
                          transformed_init[1:(length(param)-1)]))[-1]
    transformed_init[length(param)] =
      pracma::logit(transformed_init[length(param)])
  } else{
    # "Round down" if parameters are too close to the edge
    param[param>=(1-1e-8)]=.99
    param[param<1e-8]=
      rep(0.001/length(param),sum(param<1e-8))
    if(sum(param)>=(1-1e-8)){
      param = param/sum(param)*(1-1e-8)
    }
    transformed_init = ohenery::inv_smax(c(1-sum(param), param))[-1]
  }

  return(transformed_init)
}
