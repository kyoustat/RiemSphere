#' Multidimensional Scaling on Sphere
#' 
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)  # center of class 1
#' mymu2 = c(-1,0,0,0) # center of class 2
#' 
#' x1 = rvmf(50, mymu1, kappa=15)
#' x2 = rvmf(50, mymu2, kappa=15)
#' xx = rbind(x1,x2)
#' 
#' ## compute 2d embedding and visualization
#' mds2d <- sp.mds(xx, ndim=2)
#' plot(mds2d$embed[,1], mds2d$embed[,2], pch=19)
#' 
#' @references 
#' \insertRef{torgerson_multidimensional_1952}{DAS}
#' 
#' @export
sp.mds <- function(x, ndim=2, type=c("intrinsic","extrinsic")){
  ############################################################
  # Preprocessing
  if (!check_datamat(x)){
    stop("* sp.mds : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  myndim = round(ndim)
  mytype = match.arg(type)
  
  ############################################################
  # Compute Pairwise distance and Use DAS package
  dmat   = sp.pdist.internal(x, type=mytype, as.dist=TRUE)
  output = DAS::cmds(dmat, ndim=myndim)
  
  ############################################################
  # Return
  return(output)
}
  