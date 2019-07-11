# Auxiliary Computation Functions -----------------------------------------
# 01. aux_log         : mu (point) and x (point)
# 02. aux_exp         : mu (point) and v (tangent)
# 03. aux_intmean     : compute intrinsic mean returned as a vector
# 04. aux_dist_1toN   : compute intrinsic distance for a vector vs matrix
# 05. aux_rotation    : compute a rotation matrix R from 'a' to 'b'; R%*%a = b
# 06. aux_dist_MtoN   : compute intrinsic distance from M-row to N-row matrix.
# 07. aux_dist_MtoM   : copmute pairwise distance matrix
# 08. aux_assignment  : for each row, either MIN or MAX, decide an assignment
# 09. aux_clustermean : compute cluster means per class
# 10. aux_stack3d     : stack riemdata as 3d array
# 11. aux_wfrechet    : compute weighted frechet mean
# 12. aux_latent2hard : for each row, assign 1 to the largest and 0 others.

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
  return(as.vector(aux_wfrechet(myx)$x))
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
    label = as.vector(apply(mat, 1, which.max))  
  } else {
    label = as.vector(apply(mat, 1, which.min))
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
  
  output  = engine_wmean(newdata, mfdname, as.integer(maxiter), as.double(eps), weight)
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
    ids = which.max(as.vector(xx[i,]))
    if (length(ids)>1){
      ids = sample(ids,1)
    }
    output[i,ids] = 1
  }
  return(output)
}