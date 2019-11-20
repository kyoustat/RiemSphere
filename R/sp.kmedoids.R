#' k-Medoids Clustering for Data on Hypersphere
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)  # center of class 1
#' mymu2 = c(-1,0,0,0) # center of class 2
#' 
#' x1 = rvmf(50, mymu1, kappa=5)
#' x2 = rvmf(50, mymu2, kappa=10)
#' xx = rbind(x1,x2)
#' 
#' ## apply clustering with different k values
#' cl2 <- sp.kmedoids(xx, k=2)
#' cl3 <- sp.kmedoids(xx, k=3)
#' cl4 <- sp.kmedoids(xx, k=4)
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- sp.mds(xx, ndim=2)
#' mdsx  <- mds2d$embed[,1]
#' mdsy  <- mds2d$embed[,2]
#' 
#' ## compare via visualization
#' opar  <- par(mfrow=c(1,3), pty="s")
#' plot(mdsx, mdsy, col=cl2$cluster, main="k=2 medoids", pch=19)
#' plot(mdsx, mdsy, col=cl3$cluster, main="k=3 medoids", pch=19)
#' plot(mdsx, mdsy, col=cl4$cluster, main="k=4 medoids", pch=19)
#' par(opar)
#' 
#' 
#' @export
sp.kmedoids <- function(x, k=2, type=c("intrinsic","extrinsic")){
  ############################################################
  # Preprocessing
  if (!check_datamat(x)){
    stop("* sp.kmedoids : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  myk    = round(k)
  mytype = match.arg(type)
  
  ############################################################
  # Compute Pairwise Distances and do PAM
  dmat   = sp.pdist.internal(x, type=mytype, as.dist=TRUE)
  tmprun = RiemBaseExt::rclust.kmedoids(dmat, k=myk)
  
  ############################################################
  # Recap
  output = list()
  output$cluster = tmprun$clustering
  output$medoids = x[tmprun$id.med,]
  return(output)
}