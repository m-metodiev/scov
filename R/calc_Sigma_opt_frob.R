#' Minimizes the Frobenius norm via quadratic optimization
#'
#' @param matList_full      combined list of pairwise and spatial matrices
#' @param corY              a correlation matrix estimate
#' @param edge_constraints  a vector indicating which effect is constrained
#'
#' @returns a list giving the optimal parameter and corresponding matrix
#' @keywords internal
#' @importFrom quadprog          solve.QP
#'
calc_Sigma_opt_frob = function(matList_full, corY, edge_constraints=c()){
  if(length(matList_full)==0){
    print("ERROR: no matrix found")
    return(-1)
  }

  dimension = dim(matList_full[[1]])[1]
  diag(corY)=0
  for(i in (1:length(matList_full))){
    diag(matList_full[[i]]) = 0
  } # we are only comparing the covariates, not the diagonals

  # The Frobenius inner product is included in each element
  Dmat = calc_Dmat(matList_full)
  dvec =  calc_dvec(corY,matList_full)
  num_param = length(dvec)

  #browser()
  #add constraint that sum has to be smaller than 1
  Dmat = rbind(cbind(Dmat,rep(0,dim(Dmat)[1])),c(rep(0,dim(Dmat)[1]),1))
  dvec = c(dvec,0)
  Amat = rbind(cbind(rep(-1,length(matList_full)),
                     diag(1,length(matList_full))),
               c(rep(0,length(matList_full)+1)))
  bvec = -t(t(c(1,numeric(length(matList_full)))))

  num_extra_var = 1 # some variables contribute nothing and only add constraints
  Sigma_0_opt = quadprog::solve.QP(Dmat,dvec,Amat,bvec)
  init <- c(log(abs(Sigma_0_opt$solution)))[1:num_param]

  # need to stop parameters from being on the edge (messes up init)
  if(length(edge_constraints)>0){
    for(edge_constraint in edge_constraints){
      #browser()
      num_extra_var = num_extra_var + 1
      Dmat = rbind(cbind(Dmat,rep(0,dim(Dmat)[1])),c(rep(0,dim(Dmat)[1]),1))
      dvec = c(dvec,0)

      if(edge_constraint$r_min==0){

        constraint_vec = rep(-1,length(matList_full))
        constraint_digit = -1+edge_constraint$constraint_digit
      } else{
        if(edge_constraint$s_min==0){
          constraint_vec = rep(0,length(matList_full))
          constraint_vec[edge_constraint$r_min] = 1
          constraint_digit = edge_constraint$constraint_digit
        } else{
          constraint_vec = rep(0,length(matList_full))
          constraint_vec[edge_constraint$r_min] = 1
          constraint_digit = edge_constraint$constraint_digit
        }
      }
      Amat = rbind(cbind(c(constraint_vec,
                           rep(0,dim(Amat)[1]-length(matList_full))),
                         Amat),
                   c(rep(0,dim(Amat)[1]+1)))
      bvec = t(t(c(constraint_digit,bvec)))
      Sigma_0_opt = quadprog::solve.QP(Dmat,dvec,Amat,bvec)
      init <- c(log(abs(Sigma_0_opt$solution)))[1:num_param]
    }
  }
  value = Sigma_0_opt$value

  return(list(init=init,value=Sigma_0_opt$value))
}
