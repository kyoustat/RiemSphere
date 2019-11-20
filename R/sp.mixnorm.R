#' Mixture of Spherical Normal Distributions
#' 
#' @examples 
#' ## generate two-cluster data
#' mymu1 = c(0,0,0,1)        # center of class 1
#' mymu2 = c(-1,0,0,0)       # center of class 2
#' 
#' x1 = rvmf(50, mymu1, kappa=10)
#' x2 = rvmf(50, mymu2, kappa=10)
#' xx = rbind(x1,x2)
#' 
#' ## apply clustering with different k values
#' mix2 <- sp.mixnorm(xx, k=2)
#' mix3 <- sp.mixnorm(xx, k=3)
#' mix4 <- sp.mixnorm(xx, k=4)
#' 
#' ## compute 2-dimensional embedding for visualization
#' mds2d <- sp.mds(xx, ndim=2)
#' mdsx  <- mds2d$embed[,1]
#' mdsy  <- mds2d$embed[,2]
#' 
#' ## compare via visualization
#' opar  <- par(mfrow=c(1,3), pty="s")
#' plot(mdsx, mdsy, col=mix2$cluster, main="k=2 mixture", pch=19)
#' plot(mdsx, mdsy, col=mix3$cluster, main="k=3 mixture", pch=19)
#' plot(mdsx, mdsy, col=mix4$cluster, main="k=4 mixture", pch=19)
#' par(opar)
#' 
#' @export
sp.mixnorm <- function(x, k=2, init=c("kmeans","random"), maxiter=496, same.lambda=TRUE){
  ###################################################################
  ## Preprocessing 
  if (!check_datamat(x)){
    stop("* sp.mixnorm : input data should consist of at least two points on the unitsphere.")
  }

  myn = nrow(x)
  myp = ncol(x)
  myk = round(k)
  myiter = round(maxiter)
  myinit = match.arg(init)
  
  if (myk >= (myn-1)){
    stop("* sp.mixnorm : the number of clusters is too big. Try with a smaller number.")
  }
  if (myk<2){
    stop("* sp.mixnorm : the number of clusters should be at least 2.")
  }
  
  ###################################################################
  ## Initialize all the parameters
  if (all(myinit=="random")){
    initlabel = sample(c(1:myk, sample(1:myk, myn-myk, replace = TRUE)))
  } else {
    initlabel = stats::kmeans(x, myk, nstart=round(5))$cluster # label  
  }
  par.eta <- array(0,c(myn,myk))
  for (i in 1:myn){
    par.eta[i,initlabel[i]] = 1
  }
  par.mu = mixnorm.frechet(x, par.eta)
  d2mat  = (sp.pdist2.internal(x, par.mu, type="intrinsic")^2)
  if (same.lambda){
    par.lambda = mixnorm.lambda.homo(d2mat, par.eta, myp)
  } else {
    par.lambda = mixnorm.lambda.hetero(d2mat, par.eta, myp)
  }
  par.pi = base::colSums(par.eta)/myn
  loglkd.old = mixnorm.loglkd(x, par.mu, par.lambda, par.pi)

  ###################################################################
  ## Let's Iterate
  for (it in 1:myiter){
    # E-Step
    par.eta = mixnorm.eta(x, par.mu, par.lambda, par.pi)
    
    # M-Step Parameter Update
    #   Preliminary. d2mat
    d2mat = (sp.pdist2.internal(x, par.mu, type="intrinsic")^2)
    #   M1. lambda / concentration
    if (same.lambda){
      par.lambda = mixnorm.lambda.homo(d2mat, par.eta, myp)
    } else {
      par.lambda = mixnorm.lambda.hetero(d2mat, par.eta, myp)
    }
    #   M2. proportion
    par.pi = base::colSums(par.eta)/myn
    #   M3. mu / centers
    par.mu = mixnorm.frechet(x, par.eta)
    #   Mextra for computing log-likelihood
    loglkd.new  = mixnorm.loglkd(x, par.mu, par.lambda, par.pi)
    loglkd.diff = loglkd.new - loglkd.old
    # iteration over
    loglkd.old = loglkd.new
    if ((it > 5)&&(loglkd.diff <= 0)){
      break
    }
  }
  
  ###################################################################
  # Information Criterion
  if (same.lambda){
    par.k = ((myp-1)*myk) + myk + (myk-1)  
  } else {
    par.k = ((myp-1)*myk) + 1   + (myk-1)  
  }
  AIC = -2*loglkd.old + 2*par.k
  BIC = -2*loglkd.old + par.k*log(myn)
  HQIC = -2*loglkd.old + 2*par.k*log(log(myn))
  AICc = AIC + (2*(par.k^2) + 2*par.k)/(myn-par.k-1)
  
  infov = matrix(c(AIC, AICc, BIC, HQIC), nrow=1)
  colnames(infov) = c("AIC","AICc","BIC","HQIC")
  
  
  ###################################################################
  # Return
  output = list()
  output$cluster  = mixnorm.extract.label(par.eta) # 1. clutering label
  output$loglkd   = loglkd.old  # max loglkd
  output$criteria = infov       # min AIC/AICc/BIC/HQIC
  output$parameters = list(pro=par.pi, centers=par.mu, concentration=par.lambda)
  output$membership = par.eta
  return(output)
}


#   -----------------------------------------------------------------------
#' 1. compute weighted frechet mean
#' @keywords internal
#' @noRd
mixnorm.frechet <- function(dat, eta){
  n = nrow(dat)
  p = ncol(dat)
  k = ncol(eta)
  output = array(0,c(k,p))

    for (kk in 1:k){
    if (sum(eta[,kk]==1)==1){
      output[kk,] = dat[which(eta[,kk]==1),]
    } else {
      output[kk,] = as.vector(aux_wfrechet(dat, weight=as.vector(eta[,kk])))  
    }
  }
  return(output)
}
#' 2. compute lambda for homogeneous/heterogeneous
#' @keywords internal
#' @noRd
mixnorm.lambda.homo <- function(d2mat, eta, myp){
  k = ncol(eta)
  term1 = base::sum(d2mat*eta)
  term2 = base::sum(eta)
  
  myfun <- function(lambda){ # function to be minimized
    return((lambda*term1/2) + (term2*(log(dspnorm.constant(lambda,myp)))))
  }
  
  output = stats::optimize(myfun, interval=c(10*.Machine$double.eps,12345), maximum = FALSE)$minimum
  return(rep(output, k))
}
#' @keywords internal
#' @noRd
mixnorm.lambda.hetero <- function(d2mat, eta, myp){
  kk = ncol(eta)
  outvec = rep(0,kk)
  d2meta = d2mat*eta
  
  for (k in 1:kk){
    term1 = as.vector(sum(d2meta[,k]))
    term2 = as.vector(sum(eta[,k]))
    
    myfun <- function(lambda){
      return(((lambda*term1/2) + (term2*(log(dspnorm.constant(lambda,myp))))))
    }
    
    outvec[k] = stats::optimize(myfun, interval=c(10*.Machine$double.eps,12345), maximum = FALSE)$minimum
  }
  return(outvec)
}
#' 3. compute eta parameter
#' @keywords internal
#' @noRd
mixnorm.eta <- function(x, par.mu, par.lambda, par.pi){
  kk = length(par.pi)
  nn = nrow(x)
  
  par.eta = array(0,c(nn,kk))
  for (k in 1:kk){
    par.eta[,k] = as.vector(dspnorm(x, as.vector(par.mu[k,]), lambda=par.lambda[k], log = FALSE))
  }
  output = par.eta/base::rowSums(par.eta)
  return(output)
}
#' 4. compute log-likelihood
#' @keywords internal
#' @noRd
mixnorm.loglkd <- function(x, par.mu, par.lambda, par.pi){
  kk = length(par.pi)
  nn = nrow(x)
  
  output = 0
  for (k in 1:kk){
    output = output + sum(as.vector(dspnorm(x, as.vector(par.mu[k,]), lambda=par.lambda[k], log = TRUE)))
  }
  return(output)
}
#' 5. extract cluster label
#' @keywords internal
#' @noRd
mixnorm.extract.label = function(eta){
  nn = nrow(eta)
  kk = ncol(eta)
  
  lvec = rep(0,nn)
  for (n in 1:nn){
    idmins = which.min(as.vector(eta[n,]))
    if (length(idmins)==1){
      lvec[n] = idmins
    } else {
      lvec[n] = base::sample(idmins, 1)
    }
  }
  return(lvec)
}