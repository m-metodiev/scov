#' Adds combined effects to the matList via the Hadamard product
#'
#' @inheritParams combined_matList_partial_der
#'
#' @returns named list with all matList combinations and their positions
#' @keywords internal
#' @importFrom Matrix          rankMatrix
combined_matList = function(matList,
                            interaction_effects=NULL,
                            check_redundancy=FALSE){
  matList_full = c(matList$Fk,matList$Gl)

  for(i in seq_along(matList_full)){
    matList_full[[i]] = as.matrix(matList_full[[i]])
  }

  if(!is.null(matList$Gl)){
    matList_full_names = c(names(matList$Fk),"spatial")
  } else{
    matList_full_names = names(matList$Fk)
  }


  # add names of the matrices
  if(!is.null(matList$Gl)){
    base_effects = c(names(matList$Fk),"spatial")
  } else{
    base_effects = names(matList$Fk)
  }

  if(check_redundancy){
    not_redundant = (Matrix::rankMatrix(sapply(matList_full,
                                               function(x) c(x)))[1] -
                       length(matList_full)) == 0
    if(!not_redundant){
      qr_decomp = qr(sapply(matList_full, function(x) c(x)))
      warning("Redundant effects have been found.")
      warning("The following effects are not redundant:")
      warning(paste0(base_effects[qr_decomp$pivot[1:qr_decomp$rank]],
                     collapse=","))
      warning("Please use them and start over.")
    }
  } else{
    not_redundant = TRUE
  }

  counter = length(matList_full)
  sequence = seq_along(matList_full)
  matrix_pairs = list()
  found_redundant_pairs = FALSE

  if(!is.null(matList$Fk)){
    #calculate all possible Hadamard-products; Exclude global effect matrix
    for(i in sequence){
      for(j in sequence[-(1:i)]){
        matprod = matList_full[[i]] * matList_full[[j]]
        if(check_redundancy){
          not_redundant = (Matrix::rankMatrix(sapply(c(matList_full,
                                                       list(matprod)),
                                                     function(x) c(x)))[1] -
                             (length(matList_full) + 1)) == 0
        } else{
          not_redundant = TRUE
        }
        if(!is.null(interaction_effects)){
          # only loop over those
          in_left = matList_full_names[i] ==
            sapply(interaction_effects,function(interaction) interaction[1])
          in_both = sum(matList_full_names[j] ==
                          sapply(interaction_effects[in_left],
                                 function(interaction) interaction[2])) >= 1
          in_right = matList_full_names[j] ==
            sapply(interaction_effects,function(interaction) interaction[1])
          in_both = in_both |
            (sum(matList_full_names[j] ==
                   sapply(interaction_effects[in_right],
                          function(interaction) interaction[2]))>=1)
          include_matrix_pair = not_redundant & in_both
          if((!not_redundant) & in_both){
            found_redundant_pairs = TRUE
          }
        } else{
          include_matrix_pair = not_redundant
        }
        if(include_matrix_pair){
          counter = counter + 1
          matList_full[[counter]] = matprod
          matrix_pairs = c(matrix_pairs,list(c(i,j)))
        }
      }
    }
  }

  effects = base_effects
  for(matrix_pair in matrix_pairs){
    effects[length(effects)+1] = paste0(base_effects[matrix_pair[1]],
                                        " and ",
                                        base_effects[matrix_pair[2]])
  }
  names(matList_full) = effects
  if(found_redundant_pairs){
    if(length(effects) == length(base_effects)){
      warning("The interaction effects are redundant.")
      warning("Please remove them and restart.")
    } else{
      warning("Redundant effect pairs have been found.")
      warning("The following pairs are not redundant:")
      warning(paste0(effects[(length(base_effects)+1):length(effects)]),
            quote=FALSE)
      warning("Please use these effect pairs as input and start over.")
    }
    return(-1)
  } else{
    return(list(matList_full=matList_full,matrix_pairs=matrix_pairs))
  }
}
