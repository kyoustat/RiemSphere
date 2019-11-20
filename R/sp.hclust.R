#' Hierarchical Clustering for Data on Hypersphere
#' 
#' @param x an \eqn{(n\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param method the agglomeration method to be used. This must be (an unambiguous abbreviation of) one of \code{"single"},
#' \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details.
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)  # center of class 1
#' mymu2 = c(-1,0,0,0) # center of class 2
#' 
#' x1 = rvmf(20, mymu1, kappa=5)
#' x2 = rvmf(20, mymu2, kappa=5)
#' xx = rbind(x1,x2)
#' 
#' ## apply hierarchical clustering with different methods
#' hc1 <- sp.hclust(xx, method="single")
#' hc2 <- sp.hclust(xx, method="complete")
#' hc3 <- sp.hclust(xx, method="average")
#' 
#' ## visualize
#' \dontrun{
#' opar <- par(mfrow=c(1,3), pty="s")
#' plot(hc1, main="'single'")
#' plot(hc2, main="'complete'")
#' plot(hc3, main="'average'")
#' par(opar)
#' }
#' 
#' @export
sp.hclust <- function(x, type=c("intrinsic","extrinsic"), 
                      method=c("single","complete","average","mcquitty",
                               "ward.D","ward.D2","centroid","median"), members=NULL){
  ############################################################
  # Preprocessing
  if (!check_datamat(x)){
    stop("* sp.hclust : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  mytype   = match.arg(type)
  mymethod = match.arg(method)
  
  ############################################################
  # Compute Pairwise Distance
  dmat   = sp.pdist.internal(x, type=mytype, as.dist=TRUE)
  
  ############################################################
  # HCLUST via RiemBaseExt and Return
  output = RiemBaseExt::rclust.hclust(dmat)
}