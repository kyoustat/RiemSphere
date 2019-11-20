#' k-Means Clustering for Data on Hypersphere
#' 
#' @param x an \eqn{(n\times p)} row-stacked matrix for \eqn{\mathbb{S}^{p-1}}.
#' @param k the number of clusters to be found.
#' @param n.start the number of random starting-point configurations.
#' @param maxiter maximum number of iterations to be run.
#' @param type type of distance, either \code{"intrinsic"} or \code{"extrinsic"}.
#' 
#' @return a named list containing
#' \describe{
#' \item{cluster}{length-\eqn{n} vector of class membership indices.}
#' }
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)  # center of class 1
#' mymu2 = c(-1,0,0,0) # center of class 2
#' 
#' x1 = rvmf(50, mymu1, kappa=10)
#' x2 = rvmf(50, mymu2, kappa=10)
#' xx = rbind(x1,x2)
#' 
#' ## apply clustering with different k values
#' cl2 <- sp.kmeans(xx, k=2)
#' cl3 <- sp.kmeans(xx, k=3)
#' cl4 <- sp.kmeans(xx, k=4)
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- sp.mds(xx, ndim=2)
#' mdsx  <- mds2d$embed[,1]
#' mdsy  <- mds2d$embed[,2]
#' 
#' ## compare via visualization
#' opar  <- par(mfrow=c(1,3), pty="s")
#' plot(mdsx, mdsy, col=cl2$cluster, main="k=2 means", pch=19)
#' plot(mdsx, mdsy, col=cl3$cluster, main="k=3 means", pch=19)
#' plot(mdsx, mdsy, col=cl4$cluster, main="k=4 means", pch=19)
#' par(opar)
#' 
#' @export
sp.kmeans <- function(x, k=2, n.start=5, maxiter = 100, type=c("intrinsic","extrinsic")){
  ############################################################
  # Preprocessing
  if (!check_datamat(x)){
    stop("* sp.kmeans : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  myn     = nrow(x)   # number of observations
  myk     = round(k)  # desired number of clusters
  mytype  = match.arg(type)
  maxiter = round(maxiter)
  
  ############################################################
  # Initialize
  label.old  = stats::kmeans(x, myk, nstart=round(n.start))$cluster # label
  if (aux_strcmp(mytype, "intrinsic")){
    center.old = sp.kmeans.center.int(x, label.old, myk)
  } else {
    center.old = sp.kmeans.center.ext(x, label.old, myk)
  }

  ############################################################
  # Naive Algorithm
  for (it in 1:maxiter){
    # Assignment Step
    # A-1. compute pairwise distance (N x K)
    pdmat = sp.pdist2.internal(x, center.old, type=mytype)
    # A-2. class assignment
    label.new = rep(0,myn)
    for (i in 1:myn){
      idmins = which.min(as.vector(pdmat[i,]))
      if (length(idmins)==1){
        label.new[i] = idmins
      } else {
        label.new[i] = base::sample(idmins, 1)
      }
    }
    # Update Step
    if (aux_strcmp(mytype, "intrinsic")){
      center.new = sp.kmeans.center.int(x, label.new, myk)
    } else {
      center.new = sp.kmeans.center.ext(x, label.new, myk)
    } 
    # Iteration Control
    labeldel   = base::norm(as.matrix(label.old-label.new),"f")
    label.old  = label.new
    center.old = center.new
    if ((labeldel<1e-6)&&(it>=5)){
      break
    }
  }
  
  ############################################################
  # Return
  output = list()
  output$cluster = label.old
  output$centers = center.old
  return(output)
}


#   -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
sp.kmeans.center.int <- function(x, label, k){
  n = nrow(x)
  p = ncol(x)
  
  label = round(label)
  k     = round(k)
  
  centers = array(0,c(k,p)) # each row is label
  for (i in 1:k){
    idnow = which(label==i)
    if (length(idnow)==1){
      centers[i,] = x[idnow,]
    } else {
      centers[i,] = aux_intmean(x[idnow,])
    }
  }
  return(centers)
}
#' @keywords internal
#' @noRd
sp.kmeans.center.ext <- function(x, label, k){
  n = nrow(x)
  p = ncol(x)
  
  label = round(label)
  k     = round(k)
  
  centers = array(0,c(k,p)) # each row is label
  for (i in 1:k){
    idnow = which(label==i)
    if (length(idnow)==1){
      centers[i,] = x[idnow,]
    } else {
      cnow = base::colMeans(x[idnow,])
      centers[i,] = cnow/sqrt(sum(cnow^2))
    }
  }
  return(centers)
}

