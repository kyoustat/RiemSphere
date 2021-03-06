#' Pairwise Distance of Data on Hypersphere
#' 
#' @param x an \eqn{(n\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param as.dist a logical; \code{TRUE} to return an object of class \code{dist} or \code{FALSE} a symmetric matrix.
#' 
#' @return an \eqn{(n\times n)} matrix of pairwise distances or \code{dist} object.
#' 
#' @export
sp.pdist <- function(x, type=c("intrinsic","extrinsic"), as.dist=FALSE){
  ############################################################
  # Preprocessing
  # 1. check the data matrix
  if (!check_datamat(x)){
    stop("* sp.pdist : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  # 2. check the mode
  type = match.arg(type)
  # 3. parameters
  n = nrow(x)
  
  ############################################################
  # Computation
  if (all(type=="intrinsic")){
    output = cppdist_pair_int(x)
  } else {
    output = cppdist_pair_ext(x)
  }
  
  ############################################################
  # Report
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}


# internal ----------------------------------------------------------------
#' @keywords internal
#' @noRd
sp.pdist.internal <- function(x, type=c("intrinsic","extrinsic"), as.dist=FALSE){
  ############################################################
  # Preprocessing
  type = match.arg(type)
  n = nrow(x)
  
  ############################################################
  # Computation
  if (all(type=="intrinsic")){
    output = cppdist_pair_int(x)
  } else {
    output = cppdist_pair_ext(x)
  }
  
  ############################################################
  # Report
  if (as.dist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}


# # Personalized Example
# mykap = 10 # larger value, higher concentration
# mymu1 = c(0,0,0,1)
# mymu2 = c(-1,0,0,0)
# x1 = rvmf(50,mymu1,kappa=mykap)
# x2 = rvmf(50,mymu2,kappa=mykap)
# XX = rbind(x1,x2)
# 
# out.int = RiemSphere::sp.pdist(XX, type="intrinsic")
# out.ext = RiemSphere::sp.pdist(XX, type="extrinsic")
# par(mfrow=c(1,2), pty="s")
# image(out.int[,nrow(out.int):1], main="intrinsic")
# image(out.ext[,nrow(out.ext):1], main="extrinsic")