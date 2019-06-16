#' MLE for SPNORM
#' 
#' @export
mle.spnorm <- function(x, method=1){
  ## PREPROCESSING : check data
  check_datamat(x)
  if ((is.vector(x))||(nrow(x)==1)){
    stop("* mle.spnorm : computing MLE requires more than one observations.")
  }
  
  ## STEP 1. intrinsic mean
  opt.mean = aux_intmean(x)
  
  ## STEP 2. optimal lambda can be computed ?
  opt.lambda = switch (method,
    "1" = lambda_method1(x, opt.mean), # DEoptim not working well
    "2" = lambda_method2(x, opt.mean),
    "3" = lambda_method3(x, opt.mean)
  )
  
  ## RETURN
  output = list()
  output$mu = opt.mean
  output$lambda = opt.lambda
  return(output)
}


# several methods ---------------------------------------------------------
# METHOD 1. DIFFERENTIAL EVOLUTION
#' @keywords internal
#' @noRd
lambda_method1 <- function(data, mean){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  C = sum(aux_dist_1toN(mean, data)^2)
  
  # 3. compute a of log-likelihood function
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = -(lambda*C)/2
    term2 = -n*log(dspnorm.constant(lambda,D))
    return(-(term1+term2))
  }
  
  # 4. optimize a log-likelihood function with DEoptim; careful with negative sign; take the last one as the optimum
  output = as.double(tail(DEoptim(opt.fun, 10*.Machine$double.eps, 12345678, control=DEoptim.control(trace=FALSE))$member$bestmemit, n=1L))
  return(output)
}

# METHOD 2. STATS::OPTIMIZE
#' @keywords internal
#' @noRd
lambda_method2 <- function(data, mean){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  C = sum(aux_dist_1toN(mean, data)^2)
  
  # 3. compute a log-likelihood function
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = -(lambda*C)/2
    term2 = -n*log(dspnorm.constant(lambda,D))
    return(term1+term2)
  }
  
  # 4. optimize a log-likelihood function with DEoptim
  output = stats::optimize(opt.fun, interval=c(10*.Machine$double.eps,12345678), maximum=TRUE)$maximum
  return(output)
}

# METHOD 3. KISUNG'S IMPLEMENTATION OF NEWTON'S METHOD, NUMERICAL
#' @keywords internal
#' @noRd
lambda_method3 <- function(data, mean){
  # 1. parameters
  D = length(mean)
  n = nrow(data)
  
  # 2. compute a constant
  C = sum(aux_dist_1toN(mean, data)^2)
  
  # 3. compute a negative log-likelihood function to be minimized
  opt.fun <- function(lambda){
    dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
      myfunc <- function(r){
        return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
      }
      t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
      t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
      return(t1*t2)
    }
    myfun <- function(r){
      return(exp(-lambda*(r^2)/2)*((sin(r))^(D-2)))
    }
    term1 = (lambda*C)/2
    term2 = n*log(dspnorm.constant(lambda,D))
    return(term1+term2)
  }
  
  # 4. run Newton's iteration
  # 4-1. try several inputs and rough start over a grid
  candidates = sort(stats::runif(10, min=0, max=12345))
  canvals    = apply(matrix(candidates), 1, opt.fun)
  xold       = candidates[which.max(canvals)]
  # xold = 1
  
  # 4-2. run iterations
  maxiter = 1000
  for (i in 1:maxiter){
    # print(paste("iteration for Newton : ",i," initiated..", sep=""))
    h = min(abs(xold)/2, 1e-4)
    g.right = opt.fun(xold+h)
    g.mid   = opt.fun(xold)
    g.left  = opt.fun(xold-h)
    xnew    = xold - (h/2)*(g.right-g.left)/(g.right-(2*g.mid)+g.left)
    xinc    = abs(xnew-xold)
    xold    = xnew
    
    if (xinc < .Machine$double.eps^0.25){
      break
    }
  }
  
  # 5. return an optimal solution
  return(xold)
}





  
# # 
# # TESTER FOR MLE ESTIMATION -----------------------------------------------
# myp   = 5
# mylbd = stats::runif(1, min=0.0001, max=15)
# myn   = 2000
# 
# aux_log <- function(mu, x){
#   theta = base::acos(sum(x*mu))
#   if (abs(theta)<sqrt(.Machine$double.eps)){
#     output = x-mu*(sum(x*mu))
#   } else {
#     output = (x-mu*(sum(x*mu)))*theta/sin(theta)
#   }
#   return(output)
# }
# check_single <- function(x){
#   return(sum(x^2))
# }
# aux_dist_1toN <- function(x, maty){
#   dist_one <- function(y){
#     logxy = aux_log(x, y)
#     return(sqrt(sum(logxy^2)))
#   }
#   return(apply(maty, 1, dist_one))
# }
# library(RiemBase)
# mymu  = rnorm(myp)
# mymu  = mymu/sqrt(sum(mymu^2))
# myx   = RiemSphere::rspnorm(myn, mymu, lambda=mylbd)
# mymean = as.vector(rbase.mean(riemfactory(t(myx), name="sphere"))$x)
# dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
#   myfunc <- function(r){
#     return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
#   }
#   t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
#   t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
#   return(t1*t2)
# }
# 
# C = sum(aux_dist_1toN(mymean, myx)^2)
# D = ncol(myx)
# n = nrow(myx)
# 
# 
# # test 1. shape of log-likelihood function --------------------------------
# vec.lambda = seq(from=0,to=20,length.out=200)
# vec.cost   = rep(0,length(vec.lambda))
# for (i in 1:length(vec.lambda)){
#   tl = vec.lambda[i]
#   term1 = -(tl/2)*C
#   term2 = -n*log(dspnorm.constant(tl,D))
#   vec.cost[i] = term1+term2
# }
# lopt = vec.lambda[which.max(vec.cost)]
# 
# hey = mle.spnorm(myx)
# 
# plot(vec.lambda, vec.cost, main="red-MLE, blue-TRUE")
# abline(v=hey$method3, lwd=2, col="red")
# abline(v=mylbd, lwd=2, col="blue")
# 
# 
# # test 2. time comparison for concentration estimation algorithms ---------
# x11()
# library(ggplot2)
# library(microbenchmark)  # time comparison of multiple methods
# lbdtime <- microbenchmark(
#   deoptim = mle.spnorm(myx, method=1),
#   statopt = mle.spnorm(myx, method=2),
#   newtons = mle.spnorm(myx, method=3), times=20L
# )
# autoplot(lbdtime)
# 
# 
# # test 3. estimation over iterations --------------------------------------
# rec.dir <- rep(0,myn-1)
# rec.lbd <- rep(0,myn-1)
# for (i in 1:(myn-1)){
#   tgtmle = mle.spnorm(myx[1:(i+1),])
#   tgt.mean = tgtmle$mu
#   tgt.lbd  = tgtmle$lambda
# 
#   rec.dir[i] = sqrt(sum((mymu-tgt.mean)^2))
#   rec.lbd[i] = tgt.lbd
#   print(paste("iteration ",i," complete..",sep=""))
# }
# 
# selectid = round(seq(from=2,to=myn,length.out=100))-1
# xid      = 2:myn
# x11()
# par(mfrow=c(1,2))
# plot(xid[selectid], rec.dir[selectid], "b", cex=0.5, main="evolution : mean direction")
# abline(h=0, col="blue", lwd=2)
# plot(xid[selectid], rec.lbd[selectid], "b", cex=0.5, main="evolution : concentration")
# abline(h=mylbd, col="red", lwd=1.5)
