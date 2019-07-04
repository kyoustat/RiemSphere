#' Geodesic Spherical k-Means
#' 
#' 
#' @export
gskmeans <- function(x, k=2, maxiter=100){
  ## Preprocessing
  check_datamat(x)  # checking the datamatrix
  if (is.vector(x)||(nrow(x)==1)){
    stop("* gskemans : input data should consist of at least two points on the unitsphere.")
  }
  n = nrow(x)
  if (k >= (n-1)){
    stop("* gskmeans : the number of clusters is too big. Try with a smaller number.")
  }
  
  ## Initialization for both centers (ctd) and labels (label)
  naivek  = stats::kmeans(x, k)
  
  old.ctd = naivek$centers
  old.ctd = old.ctd/sqrt(rowSums(old.ctd^2))
  old.label = as.vector(as.integer(naivek$cluster))
  
  ## Main Iteration
  for (it in 1:maxiter){
    # 1. compute pairwise distances
    distmat = aux_dist_MtoN(x, old.ctd)

    # 2. update label information using minimal distance matching
    new.label = aux_assignment(distmat, max.type=FALSE) 

    # 3. update cluster means
    new.ctd = aux_clustermean(x, new.label)

    # 4. simple distance calculation and update olds
    labeldist = sum((old.label-new.label)^2)

    old.label = new.label
    old.ctd   = new.ctd
    if ((it >= 10)&&(labeldist < sqrt(.Machine$double.eps))){
      break;
    }
  }
  
  ## Return Results
  output = list()
  output$cluster = old.label
  output$centers = old.ctd
  return(output)
}




# trial -------------------------------------------------------------------
# (1) gskmeans, gskmedoids work all well compared to mixture !
#
# k <- runif(3, 10, 20)
# prob <- rep(1/3, 3)
# mu <- matrix(rnorm(3*2), ncol = 2)
# mu <- mu / sqrt( rowSums(mu^2) )
# myx <- rmixvmf(200, prob, mu, k)$x
# myd <- data.frame(x=myx[,1], y=myx[,2])
# 
# k1 <- kmeans(myx, 3)$cluster
# k2 <- gskmeans(myx, 3)$cluster
# k3 <- mix.vmf(myx, 3)$pred
# k4 <- gskmedoids(myx, 3)$cluster
# 
# myd$k1 = k1
# myd$k2 = k2
# myd$k3 = k3
# myd$k4 = k4
# 
# library(cowplot)
# f1 = ggplot(myd, aes(x=x,y=y, col=factor(myd$k1))) + geom_point(size=2.5)
# f2 = ggplot(myd, aes(x=x,y=y, col=factor(myd$k2))) + geom_point(size=2.5)
# f3 = ggplot(myd, aes(x=x,y=y, col=factor(myd$k3))) + geom_point(size=2.5)
# f4 = ggplot(myd, aes(x=x,y=y, col=factor(myd$k4))) + geom_point(size=2.5)
# plot_grid(f1,f2,f3,f4, labels=c("k-means","gskmeans","mixture","gskmedoids"), nrow=2, align="h")
