# Auxiliary Computation Functions -----------------------------------------
# 01. aux_log             : mu (point) and x (point)
# 02. aux_exp             : mu (point) and v (tangent)
# 03. aux_intmean         : compute intrinsic mean returned as a vector
# 04. aux_dist_1toN       : compute intrinsic distance for a vector vs matrix
# 05. aux_rotation        : compute a rotation matrix R from 'a' to 'b'; R%*%a = b
# 06. aux_dist_MtoN       : compute intrinsic distance from M-row to N-row matrix.
# 07. aux_dist_MtoM       : copmute pairwise distance matrix
# 08. aux_assignment      : for each row, either MIN or MAX, decide an assignment
# 09. aux_clustermean     : compute cluster means per class
# 10. aux_stack3d         : stack riemdata as 3d array
# 11. aux_wfrechet        : compute weighted frechet mean
# 12. aux_latent2hard     : for each row, assign 1 to the largest and 0 others.
# 13. aux_vmf_Apk            Apk from von-Mises Fisher
#     aux_vmf_dApk           1st derivate
#     aux_vmf_d2Apk          2nd derivate
#     aux_vmf_Rbar       
# 13. aux_vmf_apk            from the book
# 14. aux_strcmp          : strcmp of MATLAB
# 15. aux_vmf_ReML_kappa  : compute kappa given mu0 for vMF model 
# 16. aux_besselI         : Song (2012)'s truncated power-series method
# 17. aux_whichmax        : replicate 'which.is.max' function
#     aux_whichmin        : replicate 'which.is.min' function


# 01. aux_log -------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_log <- function(mu, x){
  # theta = base::acos(sum(x*mu))
  theta = tryCatch({base::acos(sum(x*mu))},
           warning=function(w){
             0
           },error=function(e){
             0
           })
  if (abs(theta)<10*(.Machine$double.eps)){
    output = x-mu*(sum(x*mu))
  } else {
    output = (x-mu*(sum(x*mu)))*theta/sin(theta)
  }
  return(output)
}


# 02. aux_exp -------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_exp <- function(x, d){
  nrm_td = norm(matrix(d),"f")
  if (nrm_td < sqrt(.Machine$double.eps)){
    output = x;
  } else {
    output = cos(nrm_td)*x + (sin(nrm_td)/nrm_td)*d; 
  }
  return(output)
}



# 03. aux_intmean ---------------------------------------------------------
#' @keywords internal
#' @noRd
aux_intmean <- function(myx){
  output = aux_wfrechet(myx) 
  if (is.list(output)){
    return(as.vector(output$x))
  } else {
    return(as.vector(output))
  }
}



# 04. aux_dist_1toN -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_dist_1toN <- function(x, maty){
  dist_one <- function(y){
    logxy = aux_log(x, y)
    return(sqrt(sum((logxy)^2)))
  }
  return(as.vector(apply(maty, 1, dist_one)))
}


# 05. aux_rotation --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_rotation <- function(a, b) {
  p <- length(a)
  ab <- sum(a * b)
  ca <- a - b * ab
  ca <- ca/sqrt(sum(ca^2))
  A <- b %*% t(ca)
  A <- A - t(A)
  theta <- acos(ab)
  diag(p) + sin(theta) * A + (cos(theta) - 1) * (b %*% 
                                                   t(b) + ca %*% t(ca))
}

# 06. aux_dist_MtoN -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_dist_MtoN <- function(matM, matN){
  M = nrow(matM)
  N = nrow(matN)
  output = array(0,c(M,N))
  for (i in 1:M){
    x = as.vector(matM[i,])
    output[i,] = aux_dist_1toN(x, matN)
  }
  return(output)
}


# 07. aux_dist_MtoM -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_dist_MtoM <- function(mat){
  n = nrow(mat)
  output = array(0,c(n,n))
  for (i in 1:(n-1)){
    tgtx = as.vector(mat[i,])
    tgtM = mat[(i+1):n,]
    if (i==(n-1)){
      tgtM = matrix(tgtM,nrow=1)
    }
    tgtd = aux_dist_1toN(tgtx,tgtM)
    
    output[i,(i+1):n] = tgtd
    output[(i+1):n,i] = tgtd
  }
  return(output)
}

# 08. aux_assignment ------------------------------------------------------
#' @keywords internal
#' @noRd
aux_assignment <- function(mat, max.type=TRUE){
  # for each row, use 'which.max' or 'which.max'
  if (isTRUE(max.type)){
    label = as.vector(apply(mat, 1, aux_whichmax))  
  } else {
    label = as.vector(apply(mat, 1, aux_whichmin))
  } 
  return(label)
}

# 09. aux_clustermean -----------------------------------------------------
#' @keywords internal
#' @noRd
aux_clustermean <- function(dat, assignment){
  # 1. make group labels
  ulabel = sort(unique(assignment), decreasing = FALSE)
  nlabel = length(ulabel)
  groups = list()
  for (i in 1:nlabel){
    groups[[i]] = which(assignment==ulabel[i])
  }
  
  # 2. compute means
  centers = array(0,c(nlabel,ncol(dat)))
  for (i in 1:nlabel){
    tgtlabel = groups[[i]]
    if (length(tgtlabel)==1){
      centers[i,] = as.vector(dat[tgtlabel,])
    } else {
      centers[i,] = aux_intmean(dat[tgtlabel,])
    }
  }
  
  # 3. return the result
  return(centers)
}

# 10. aux_stack3d ---------------------------------------------------------
#' @keywords internal
#' @noRd
aux_stack3d <- function(riemdata){
  msize = riemdata$size
  ndata = length(riemdata$data)
  
  matdata = array(0,c(msize[1], msize[2], ndata))
  for (i in 1:ndata){
    matdata[,,i] = (riemdata$data[[i]])
  }
  return(matdata)
}


# 11. aux_wfrechet --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_wfrechet <- function(dat, weight=rep(1,nrow(dat))/nrow(dat)){
  eps = 1e-6;
  maxiter = 496;
  
  rdat    = RiemBase::riemfactory(t(dat), name="sphere")
  newdata = aux_stack3d(rdat)  # arg 1. 
  mfdname = tolower(rdat$name) # arg 2.
  
  output  = as.vector(engine_wmean(newdata, mfdname, as.integer(maxiter), as.double(eps), weight)$x)
  return(output)
}
# TEST 1. rbase.mean may be incorrect
# library(Directional)
# mymu = rnorm(4); mymu = mymu/sqrt(sum(mymu^2));
# mykk = 1.0;
# 
# diff1 = rep(0,100)
# diff2 = rep(0,100)
# seqnn = round(seq(from=10,to=500,length.out=100))
# for (i in 1:100){
#   x = rspnorm(seqnn[i], mymu, mykk)
# 
#   out1 = as.vector(aux_wfrechet(x)$x)
#   out2 = as.vector(RiemBase::rbase.mean(riemfactory(t(x), name="sphere"))$x)
# 
#   diff1[i] = sum((mymu-out1)^2)
#   diff2[i] = sum((mymu-out2)^2)
#   print(paste("iteration ",i,"/100 complete..",sep=""))
# }
# par(mfrow=c(1,2))
# plot(seqnn,diff1,"l")
# plot(seqnn,diff2,"l")

# # TEST 2. engine_wmean's average performance
# mymu = rnorm(4); mymu = mymu/sqrt(sum(mymu^2));
# mykk = 1.0;
# diff  = rep(0,100)
# seqnn = round(seq(from=10,to=400,length.out=100))
# for (i in 1:100){
#   tmp = rep(0,10)
#   for (j in 1:10){
#     x   = rspnorm(seqnn[i], mymu, mykk)
#     out = as.vector(aux_wfrechet(x)$x)
#     tmp[j] = sqrt(sum((mymu-out)^2))
#   }
#   diff[i] = mean(tmp)
#   print(paste("iteration ",i," complete..",sep=""))
# }
# graphics.off()
# plot(seqnn,diff,"l")


# 12. aux_latent2hard -----------------------------------------------------
#' @keywords internal
#' @noRd
aux_latent2hard <- function(xx){
  n = nrow(xx)
  p = ncol(xx)
  output = array(0,c(n,p))
  for (i in 1:n){
    ids = aux_whichmax(as.vector(xx[i,]))
    output[i,ids] = 1
  }
  return(output)
}



# 13. aux_ for von-Mises  -------------------------------------------------
#' @keywords internal
#' @noRd
aux_vmf_Apk <- function(p, kappa){
  if (kappa > 100){
    output = tryCatch({
      exp(besselIasym(kappa, p/2, log=TRUE) - besselIasym(kappa, (p/2)-1, log=TRUE))
    }, error = function(e){
      exp(aux_besselI(kappa, p/2)-aux_besselI(kappa, (p/2)-1))
    }, warning = function(w){
      exp(aux_besselI(kappa, p/2)-aux_besselI(kappa, (p/2)-1))
    })
  } else if (p > 100){
    output = tryCatch({
      exp(besselI.nuAsym(kappa, p/2, 5, log=TRUE)-besselI.nuAsym(kappa, (p/2)-1, 5, log=TRUE))
    }, error = function(cond){
      exp(aux_besselI(kappa, p/2)-aux_besselI(kappa, (p/2)-1))
    }, warning = function(){
      exp(aux_besselI(kappa, p/2)-aux_besselI(kappa, (p/2)-1))
    })
  } else {
    output = exp(log(besselI(kappa, p/2)) - log(besselI(kappa,(p/2)-1)))
  }
  return(output)
}
#' @keywords internal
#' @noRd
aux_vmf_dApk <- function(p, kappa){  # first derivative
  Apk = aux_vmf_Apk(p, kappa)
  return(1 - (Apk^2) - ((p-1)*Apk/kappa))
}
#' @keywords internal
#' @noRd
aux_vmf_d2Apk <- function(p, kappa){ # 2nd derivative
  Apk = aux_vmf_Apk(p, kappa)
  term1 = 2*(Apk^3)
  term2 = 3*((p-1)/kappa)*(Apk^2)
  term3 = ((p*(p-1) - 2*(kappa^2))/(kappa^2))*Apk
  term4 = -(p-1)/kappa
  return(term1+term2+term3+term4)
}
#' @keywords internal
#' @noRd
aux_vmf_Rbar <- function(dat){ # (n x p) convention
  xbar = colMeans(dat)
  return(sqrt(sum(xbar^2)))
}
#' @keywords internal
#' @noRd
aux_vmf_apk <- function(p, kappa){
  term1 = (1-(p/2))*(log(kappa)-log(2))
  term2 = lgamma(p/2)
  term3 = aux_besselI(kappa, (p/2)-1)

  output = term1+term2+term3
  return(output)
}


# 14. aux_strcmp ----------------------------------------------------------
#' @keywords internal
#' @noRd
aux_strcmp <- function(s1, s2) {
  if (!is.vector(s1, mode="character") || !is.vector(s1, mode="character"))
    stop("Arguments 's1' and 's2' must be character vectors.")
  
  if (length(s1) == length(s2)){
    return(all(s1 == s2))
  } else {
    return(FALSE)
  }
}

# 15. aux_vmf_ReML_kappa --------------------------------------------------
#' @keywords internal
#' @noRd
aux_vmf_ReML_kappa <- function(dat, mu0){
  # parameters
  n = nrow(dat)
  p = ncol(dat)
  xbar = as.vector(colMeans(dat))
  C    = sum(xbar*as.vector(mu0))
  
  min.fun <- function(kappa){
    return(aux_vmf_apk(p, kappa) - kappa*C)
  }
  
  kap.old = vmf_2005banerjee(dat)
  maxiter = 123
  epsthr  = 1e-6
  epserr.ReML  = 100
  citer   = 1
  while (epserr.ReML > epsthr){
    h = min(abs(kap.old),1e-3)/2
    fr = min.fun(kap.old + h)
    fc = min.fun(kap.old)
    fl = min.fun(kap.old - h)
    
    kap.new = kap.old - (h/2)*(fr-fl)/(fr-(2*fc)+fl)
    epserr.ReML  = abs((kap.old - kap.new))
    kap.old = kap.new
    citer   = citer + 1
    if (citer > maxiter){
      break
    }
  }
  return(kap.old)
}


# 16. aux_besselI ---------------------------------------------------------
#' @export
aux_besselI <- function(x, s, log=TRUE){
  c0 = 1.000000000190015
  c1 = 76.18009172947146
  c2 = -86.50532032941677
  c3 = 24.01409824083091
  c4 = 1.231739572450155
  c5 = 1.208650973866179*(1e-3)
  c6 = -5.395239384953*(1e-6)
  
  t1 = c0 + c1/(s+1) + c2/(s+2) + c3/(s+3) + c4/(s+4) + c5/(s+5) + c6/(s+6)
  t2 = (s+0.5)*log(s+5.5) - (s+5.5)
  t1 = s*log(0.5*x) - 0.5*log(2*pi) - t2 - log(t1)
  R = 1.0
  M = 1.0
  k = 1
  
  epsthr  = 1e-6
  epserr.bessel  = 1000
  maxiter = 500
  citer   = 1
  x2      = (x^2)
  while (epserr.bessel > epsthr){
    citer = citer + 1
    R = R*(0.25*(x2))/((s+k)*k)
    M = M + R
    epserr.bessel = R/M
    k = (k+1)
    citer = citer + 1
    if (citer > maxiter){
      break
    }
  }
  output = t1 + log(M)
  if (log==FALSE){
    output = exp(output)
  }
  return(output)
}


# 17. aux_whichmax --------------------------------------------------------
#' @keywords internal
#' @noRd
aux_whichmax <- function(vec){
  mval = base::max(vec)
  idlarge = which(vec>=mval)
  if (length(idlarge)==1){
    return(idlarge)
  } else {
    return(base::sample(idlarge, 1))
  }
}
#' @keywords internal
#' @noRd
aux_whichmin <- function(vec){
  mval = base::min(vec)
  idlarge = which(vec<=mval)
  if (length(idlarge)==1){
    return(idlarge)
  } else {
    return(base::sample(idlarge, 1))
  }
}