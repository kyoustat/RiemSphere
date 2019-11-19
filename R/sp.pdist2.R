#' Pairwise Distance of Two Sets of Data on Hypersphere
#' 
#' @param x1 an \eqn{(m\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.
#' @param x2 an \eqn{(n\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' 
#' @return an \eqn{(m \times n)} matrix of cross distances between \code{x1} and \code{x2}.
#' @export
sp.pdist2 <- function(x1, x2, type=c("intrinsic","extrinsic")){
  ############################################################
  # Preprocessing
  # 1. check the data matrix
  if (!check_datamat(x1)){
    stop("* sp.pdist2 : an input 'x1' is not a row-stacked matrix of unit-norm vectors.")
  }
  if (!check_datamat(x2)){
    stop("* sp.pdist2 : an input 'x2' is not a row-stacked matrix of unit-norm vectors.")
  }
  if (ncol(x1)!=ncol(x2)){
    stop("* sp.pdist2 : two inputs 'x1' and 'x2' should have same number of columns.")
  }
  # 2. check the mode
  type = match.arg(type)
  # 3. parameters
  m = base::nrow(x1)
  n = base::nrow(x2)
  
  ############################################################
  # Computation
  if (all(type=="intrinsic")){
    outmat = aux_dist_MtoN(x1, x2)
  } else if (all(type=="extrinsic")){
    outmat = array(0,c(m,n))
    for (i in 1:m){
      tgti = as.vector(x1[i,])
      for (j in 1:n){
        tgtj = as.vector(x2[j,])
        outmat[i,j] = sqrt(sum((tgti-tgtj)^2))
      }
    }
  }
  
  ############################################################
  # Report
  return(outmat)
}