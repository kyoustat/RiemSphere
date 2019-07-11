#' Mixture of
#' 
#' @export
mix.spnorm <- function(x, k=2, n.start=20, hard.assign=FALSE, maxiter=496, eps=1e-5){
  #-----------------------------------------------------------------------------------------
  ## Preprocessing 
  check_datamat(x)  # checking the datamatrix
  if (is.vector(x)||(nrow(x)==1)){
    stop("* mix.spnorm : input data should consist of at least two points on the unitsphere.")
  }
  n = nrow(x)
  if (k >= (n-1)){
    stop("* mix.spnorm : the number of clusters is too big. Try with a smaller number.")
  }
  if (k<2){
    stop("* mix.spnorm : the number of clusters should be at least 2.")
  }
  d = ncol(x) # S^{d-1) in R^p
  old.cluster <- as.integer(as.vector(stats::kmeans(x, k, nstart=n.start)$cluster)) # initialize the label
  
  #-----------------------------------------------------------------------------------------
  ## Iteration
  for (it in 1:maxiter){
    # E-STEP
    
    # Hard Assignment if necessary
    # M-STEP
  }
}



# Auxiliary Functions for MIX.SPNORM --------------------------------------
# 1. evaluating mixture density and return (log)-likelihood for the whole
#' @keywords internal
#' @noRd
mix.spnorm.eval <- function(x, centers, concentrations, weights, log=FALSE){
  n = nrow(x)
  k = nrow(centers)
  if (length(concentrations)==1){
    concentrations = rep(concentrations,k)
  }
  
  densities = array(0,c(n,k))
  for (i in 1:k){
    densities[,i] = dspnorm(x, as.vector(centers[i,]), concentrations[i])*weights[i]
  }
  density.final = rowSums(densities)
  if (log){
    return(sum(log(density.final)))
  } else {
    return(prod(density.final))
  }
}
