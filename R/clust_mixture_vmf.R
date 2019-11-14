#' Clustering with Mixture of Distributions on Hypersphere
#' 
#' 
#' @export
sp.mixclust <- function(x, k=2, n.start=20, kernel=c("vmf","spnorm")){
  #-----------------------------------------------------------------------------------------
  ## Preprocessing 
  check_datamat(x)  # checking the datamatrix
  if (is.vector(x)||(nrow(x)==1)){
    stop("* sp.mixclust : input data should consist of at least two points on the unitsphere.")
  }
  p = ncol(x)
  
  #-----------------------------------------------------------------------------------------
  # Apply the method
  myk = round(k)
  nst = round(n.start)
  tmp = Directional::mix.vmf(x, g=myk, n.start=nst)
  
  #-----------------------------------------------------------------------------------------
  # Separate temporary outputs
  
}