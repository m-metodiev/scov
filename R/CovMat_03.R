#' Computes the correlation matrix corresponding to the SCE model
#'
#' @param parm                  the parameter
#' @param matList               a list of the matrices (pairwise+spatial)
#' @param adj_positions         position of the spatial effect (if embedded)
#' @param combined_effects      if yes, computes the Fosdick-Raftery model
#' @param interaction_effects   lists the interaction effects
#'
#' @keywords internal
#' @returns a named list including the following elements:
#'            Phi:                sum of scaled pairwise effect matrices
#'            Gamma:              the null-matrix (no longer used)
#'            matList_combined    list of all matrices
#'                                (pairwise+spatial+interaction)
#'            Sigma               the SCE-model correlation matrix
#'            alpha_delta         linear effects (pairwise+spatial)
#'
#' the correlation matrix corresponding to the SCE model
CovMat_03 <- function(parm, matList, adj_positions, combined_effects=FALSE,
                      interaction_effects=list()){
  if(!is.null(matList$Fk)){
    dimension <- ncol(matList$Fk[[1]])
    alpha <- parm[1:length(matList$Fk)]
  } else{
    dimension <- length(adj_positions)
  }


  if(combined_effects=="FosdickRaftery"){
    for(k in 1:length(matList$Fk)){Phi <- Phi + alpha[k]*matList$Fk[[k]]}
    adj_matrix = matList$Al[[1]]
    adj_matrix[is.na(adj_matrix)]=0
    adj_matrix = (adj_matrix != 0) + 0

    Gamma <- delta[1] * adj_matrix[adj_positions,adj_positions]
    Sigma <- .5*(Phi + Gamma + t(Phi+Gamma))
    diag(Sigma) = 1
    return(list(Phi=Phi, Gamma=Gamma, Sigma=Sigma))
  }

  if(!is.null(matList$Al)){
    delta <- parm[length(matList$Fk)+(1:length(matList$Ml))]
    betal = parm[length(matList$Fk)+length(matList$Ml)+(1:length(matList$Ml))]
    s = dim(matList$Ml[[1]])[1]
    matList$Gl[[1]] = tilde_G_inv(matList$Ml[[1]],
                                  matList$Al[[1]],
                                  parm[length(parm)])[adj_positions,
                                                      adj_positions]

  }
  matList_combined = combined_matList(matList,
                                      interaction_effects)$matList_full

  alpha_delta = parm[1:length(matList_combined)]
  Phi <- Gamma <- matrix(0, dimension, dimension)
  for(k in 1:length(matList_combined)){
    Phi <- Phi + alpha_delta[k]*matList_combined[[k]]
  }
  # Sigma is just a convex combination of all matrices
  Sigma <- .5*(Phi + t(Phi));
  # assure symmetry (due to numerical issues)

  diag(Sigma) = 1
  return(list(Phi=Phi,
              Gamma=Gamma,
              Sigma=Sigma,
              matList_combined=matList_combined,
              alpha_delta=alpha_delta))
}
