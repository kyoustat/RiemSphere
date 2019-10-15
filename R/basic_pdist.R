#' Pairwise Distance
#' 
#' @param x (n x p) data matrix on S^{p-1}
#' 
#' @export
pdist <- function(x, mode=c("intrinsic","extrinsic"), as.dist=FALSE){
  ############################################################
  # Preprocessing
  # 1. check the data matrix
  if (!check_datamat(x)){
    stop("* pdist : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  # 2. check the mode
  mode = match.arg(mode)
  # 3. parameters
  n = nrow(x)
  
  ############################################################
  # Computation
  if (all(mode=="intrinsic")){
    output = array(0,c(n,n))
    for (i in 1:(n-1)){
      tgt1 = as.vector(x[i,])
      for (j in (i+1):n){
        tgt2  = as.vector(x[j,])
        logxy = aux_log(tgt1, tgt2)
        output[i,j] <- output[j,i] <- sqrt(sum((logxy)^2))
      }
    }
  } else {
    output = as.matrix(stats::dist(x))
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
# out.int = RiemSphere::pdist(XX, mode="intrinsic")
# out.ext = RiemSphere::pdist(XX, mode="extrinsic")
# par(mfrow=c(1,2), pty="s")
# image(out.int[,nrow(out.int):1], main="intrinsic")
# image(out.ext[,nrow(out.ext):1], main="extrinsic")