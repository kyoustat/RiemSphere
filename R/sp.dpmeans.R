#' DP-Means Clustering for Data on Hypersphere
#' 
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
#' ## apply clustering with different threshold values
#' case1 <- sp.dpmeans(xx, lambda=0.5)
#' case2 <- sp.dpmeans(xx, lambda=1.0)
#' case3 <- sp.dpmeans(xx, lambda=1.5)
#' case4 <- sp.dpmeans(xx, lambda=2.0)
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- sp.mds(xx, ndim=2)
#' mdsx  <- mds2d$embed[,1]
#' mdsy  <- mds2d$embed[,2]
#' 
#' ## compare via visualization
#' opar  <- par(mfrow=c(2,2), pty="s")
#' plot(mdsx, mdsy, col=case1$cluster, main="DP lambda=0.5", pch=19)
#' plot(mdsx, mdsy, col=case2$cluster, main="DP lambda=1.0", pch=19)
#' plot(mdsx, mdsy, col=case3$cluster, main="DP lambda=1.5", pch=19)
#' plot(mdsx, mdsy, col=case4$cluster, main="DP lambda=2.0", pch=19)
#' par(opar)
#' 
#' @export
sp.dpmeans <- function(x, lambda=1, type=c("intrinsic","extrinsic"), 
                       maxiter=1234, abstol=1e-6, permute.order=FALSE){
  ############################################################
  # Preprocessing
  if (!check_datamat(x)){
    stop("* sp.dpmeans : an input 'x' is not a row-stacked matrix of unit-norm vectors.")
  }
  myn     = nrow(x)   # number of observations
  mytype  = match.arg(type)
  maxiter = round(maxiter)

  ############################################################
  # Initialize
  lambda  = as.double(lambda)
  # if (lambda > (pi)){
  #   warning(" sp.dpmeans : if threshold 'lambda' is greater than 3.141592..., it gives a single large component.")
  # }
  mu      = as.vector(aux_intmean(x)) # global mean
  labels  = rep(1,myn)                # labels={1,2,...,n}
  k       = 1
  p       = ncol(x)
  
  ############################################################
  # Main Iteration
  ss.old = sp.compute.ss(x, labels, mu, mytype)
  ss.new = 0
  for (iter in 1:maxiter){
    # 0. updating order of observations
    if (permute.order){
      idseq = sample(1:myn)
    } else {
      idseq = 1:myn
    }
    # 1. update the class membership per each class
    for (i in idseq){
      # 1-1. compute distances to the centers
      # dic = rep(0, k); for (j in 1:k){dic[j] = sum((as.vector(data[i,])-as.vector(mu[j,]))^2)}
      if ((is.vector(mu))||(nrow(mu)==1)){
        if (aux_strcmp(mytype, "extrinsic")){
          dic = as.vector(cppdist_ext_1toN(as.vector(x[i,]), matrix(mu,nrow=1)))  
        } else {
          dic = as.vector(cppdist_int_1toN(as.vector(x[i,]), matrix(mu,nrow=1)))  
        }
      } else {
        if (aux_strcmp(mytype, "extrinsic")){
          dic = as.vector(cppdist_ext_1toN(as.vector(x[i,]), mu))  
        } else {
          dic = as.vector(cppdist_int_1toN(as.vector(x[i,]), mu))  
        }
      }
      # 1-2. assign new or stay
      if (min(dic) > lambda){
        k = k+1
        labels[i] = k
        mu = rbind(mu, x[i,])
      } else {
        labels[i] = aux_whichmin(dic)
        # idmins = which(dic==min(dic))
        # if (length(idmins)>1){
        #   labels[i] = sample(idmins, 1)
        # } else {
        #   labels[i] = idmins 
        # }
      }
    }
    # 2. rearrange the label (remove empty ones)
    labels = as.factor(labels)
    ulabel = sort(unique(labels))
    labnew = rep(0,myn)
    for (i in 1:length(ulabel)){
      labnew[(labels==ulabel[i])] = i
    }
    labels = labnew
    k      = round(max(labels))
    # 3. compute per-class means
    uassign = sort(unique(labels))
    mu = array(0,c(k,p))
    for (i in 1:k){
      idmean = which(labels==uassign[i])
      if (length(idmean)==1){
        mu[i,] = as.vector(x[idmean,])
      } else {
        mu[i,] = as.vector(aux_intmean(x[idmean,]))
      }
    }
    # 4. compute updated sum of squared errors
    ss.new   = sp.compute.ss(x, labels, mu, mytype)
    ss.delta = ss.old-ss.new
    ss.old   = ss.new
    # 5. stop if updating is not significant
    if ((ss.delta < abstol)&&(iter >= min(10,maxiter))){
      break
    }
    # print(paste("iteration ",iter,"/",maxiter," complete with ss=",ss.new,sep=""))
  }
  
  ############################################################
  # Return the results
  output = list()
  output$cluster = labels
  output$centers = mu
  return(output)
}


# auxiliary functions -----------------------------------------------------
#' @keywords internal
#' @noRd
sp.compute.ss <- function(data, label, centers, dtype){
  p = ncol(data)
  if (is.vector(centers)){
    centers = matrix(centers, nrow=1)
  }
  ulabel = sort(unique(label))
  output = 0
  for (i in 1:length(ulabel)){
    subdata = data[which(label==ulabel[i]),]
    if (is.vector(subdata)){
      output = output + sum(as.vector(sp.pdist2.internal(matrix(subdata,nrow=1),matrix(centers[i,],nrow=1),type=dtype))^2)
    } else {
      output = output + sum(as.vector(sp.pdist2.internal(subdata,matrix(centers[i,],nrow=1),type=dtype))^2)
    }
  }
  return(output)
}

  
