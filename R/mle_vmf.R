#' MLE for von Mises-Fisher
#' 
#' 'Marginal Approximation (MA)' ; 1,2,3
#' 'Marginal Exact         (ME)' ; 1,2,3
#' 
#' @export
mle.vmf <- function(x, method=c("Banerjee","Christie","nmarg1","nmarg2","nmle1","nmle2",
                                "Song","Sra","Tanabe","Uniroot")){
  ########################################################################
  ## PREPROCESSING : check data
  check_datamat(x)
  if ((is.vector(x))||(nrow(x)==1)){
    stop("* mle.vmf : computing MLE requires more than one observations.")
  }
  if (missing(method)){
    method = "Banerjee"
  }
  method = tolower(method)
  alldip = c("banerjee","sra","tanabe","uniroot","song","nmarg1","nmarg2","me1","me2",
             "nmle1","nmle2","christie") # me1 and me2 are just saved
  method = match.arg(method, alldip)
  
  
  ########################################################################
  ## Compute 1 : extrinsic mean
  opt.mean = as.vector(colMeans(x))
  opt.mean = opt.mean/sqrt(sum(opt.mean^2))
  
  ########################################################################
  ## Compute 2 : kappa concentration
  opt.kappa <- switch (method,
    "banerjee" = vmf_2005banerjee(x),
    "sra"      = vmf_2012sra(x),
    "tanabe"   = vmf_2007tanabe(x),
    "uniroot"  = vmf_uniroot(x),
    "song"     = vmf_2012song(x),
    "nmarg1"   = vmf_ma1(x), # NMARG<-MA's : finite difference approximation
    "nmarg2"   = vmf_ma2(x),
    "me1"      = vmf_me1(x), # ME: marginal exact
    "me2"      = vmf_me2(x),
    "nmle1"    = vmf_nmle1(x),
    "nmle2"    = vmf_nmle2(x),
    "christie" = vmf_christie(x)
  )
  
  ########################################################################
  ## Report Output
  output = list()
  output$mu    = opt.mean
  output$kappa = opt.kappa
  return(output)
}






######## SEVERAL KAPPA APPROXIMATION METHODS ##############################
# 1. 2005 Banerjee --------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2005banerjee <- function(x){
  rbar <- aux_vmf_Rbar(x)
  p    <- ncol(x)
  
  term1 = rbar*(p-(rbar^2))
  term2 = 1-(rbar^2)
  return(term1/term2)
}


# 2. 2007 Tanabe ----------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2007tanabe <- function(x){
  p       = ncol(x)
  rbar    = aux_vmf_Rbar(x)
  kap.old = vmf_2005banerjee(x)  # initialization
  
  maxiter = 123
  reltol  = 1e-7
  valtol  = 12345
  citer   = 1
  
  while (valtol > reltol){
    kap.new = rbar*kap.old/aux_vmf_Apk(p, kap.old)
    valtol  = abs(kap.new-kap.old)/kap.old
    citer = citer + 1
    
    kap.old = kap.new
    if (citer >= maxiter){
      break
    }
  }
  return(kap.old)
}


# 3. 2012 Sra -------------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2012sra <- function(x){
  p    = ncol(x)
  rbar = aux_vmf_Rbar(x)
  kap0 = vmf_2005banerjee(x)     # initialization
  kap0Apk = aux_vmf_Apk(p, kap0) # A_p(kap0)
  
  # first iteration
  kap1 = kap0 - ((kap0Apk-rbar)/(1-(kap0Apk^2)-(((p-1)/kap0)*kap0Apk)))
  kap1Apk = aux_vmf_Apk(p, kap1)
  
  # second iteration
  kap2 = kap1 - ((kap1Apk-rbar)/(1-(kap1Apk^2)-(((p-1)/kap1)*kap1Apk)))
  return(kap2)
  
}


# 4. uniroot --------------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_uniroot <- function(x){
  rbar = aux_vmf_Rbar(x)
  p    = ncol(x)
  
  myfun <- function(kappa){
    return(aux_vmf_Apk(p, kappa)-rbar)
  }
  tgt = vmf_2005banerjee(x)
  return(stats::uniroot(myfun, lower=tgt/2, upper=tgt*2)$root)
}


# 5. 2012 Song ------------------------------------------------------------
#' @keywords internal
#' @noRd
vmf_2012song <- function(x){
  p = ncol(x)
  kap0 = vmf_2005banerjee(x)
  rbar = aux_vmf_Rbar(x)
  
  # iteration 1.
  f0   = aux_vmf_Apk(p, kap0)-rbar
  df0  = aux_vmf_dApk(p, kap0)
  ddf0 = aux_vmf_d2Apk(p, kap0)
  kap1 = kap0 - (2*f0*df0)/(2*(df0^2) - (f0*ddf0))
  
  # iteration 2.
  f1   = aux_vmf_Apk(p, kap1)-rbar
  df1  = aux_vmf_dApk(p, kap1)
  ddf1 = aux_vmf_d2Apk(p, kap1)
  kap2 = kap1 - (2*f1*df1)/(2*(df1^2) - (f1*ddf1))
  return(kap2)
}


# 6. numerical marginal mle -----------------------------------------------
#' @keywords internal
#' @noRd
vmf_ma1 <- function(x){
  n = nrow(x)
  p = ncol(x)
  kap.old = vmf_2005banerjee(x)
  Rbar = aux_vmf_Rbar(x)
  
  myfun <- function(kappa){
    return(aux_vmf_Apk(p,kappa) - Rbar*aux_vmf_Apk(p, n*Rbar*kappa))
  }
  
  maxiter = 123
  epsthr  = 1e-8
  sqrteps = sqrt(.Machine$double.eps)
  if (Rbar < 1/sqrt(n)){
    return(0)
  } else {
    incval = 100000
    citer  = 1
    while (incval > epsthr){
      hh      = min(abs(kap.old),sqrteps)/2
      
      ## Critical Part : updating rules
      fx0 = myfun(kap.old) # evaluation
      fx1 = (myfun(kap.old+hh)-myfun(kap.old-hh))/(2*hh)
      kap.new = kap.old-(fx0/fx1)
      
      # update
      incval  = abs(kap.new-kap.old)
      citer   = citer + 1  # update the iteration
      kap.old = kap.new
      if (citer > maxiter){
        break
      }
    }
    return(kap.old)
  }
}

#' @keywords internal
#' @noRd
vmf_ma2 <- function(x){
  n = nrow(x)
  p = ncol(x)
  kap.old = vmf_2005banerjee(x)
  Rbar = aux_vmf_Rbar(x)
  
  myfun <- function(kappa){
    return(aux_vmf_Apk(p,kappa) - Rbar*aux_vmf_Apk(p, n*Rbar*kappa))
  }
  
  maxiter = 123
  epsthr  = 1e-8
  sqrteps = sqrt(.Machine$double.eps)
  if (Rbar < 1/sqrt(n)){
    return(0)
  } else {
    incval = 100000
    citer  = 1
    while (incval > epsthr){
      hh      = min(abs(kap.old),sqrteps)/2
      
      ## Critical Part : updating rules
      fxx = myfun(kap.old)     # center point
      fxl = myfun(kap.old-hh)  # right 
      fxr = myfun(kap.old+hh)  # left
      
      fx0 = fxx                    # f(x)
      fx1 = (fxr-fxl)/(2*hh)       # f'(x)
      fx2 = (fxr-2*fxx+fxl)/(hh^2)  # f''(x)

      term1 = 2*fx0*fx1
      term2 = (2*(fx1^2) - fx0*fx2)
      kap.new = kap.old - term1/term2
      
      # update
      incval  = abs(kap.new-kap.old)
      citer   = citer + 1  # update the iteration
      kap.old = kap.new
      if (citer > maxiter){
        break
      }
    }
    return(kap.old)
  }
}

# 7. exact marginal mle ---------------------------------------------------
#' @keywords internal
#' @noRd
vmf_Apakdiff <- function(Apak, a, p, kappa, order=1){
  if(order==1){
    return(a - ((p-1)/kappa)*Apak - a*(Apak^2))
  } else if (order==2){
    return(-(a*(p-1)/kappa) + (((p^2)-p)/(kappa^2) - 2*(a^2))*Apak + 3*a*((p-1)/kappa)*(Apak^2) + 2*(a^2)*(Apak^3))
  } else if (order==3){
    term1 = a*(p^2-p)/(kappa^2) - 2*(a^3)
    term2 = ((p-p^3)/(kappa^3) + 8*(a^2)*(p-1)/kappa)*Apak
    term3 = (8*a^3 - (a*(7*p^2 - 13*p + 6))/(kappa^2))*(Apak^2)
    term4 = -(12*(a^2)*(p-1)/kappa)*(Apak^3)
    term5 = -6*(a^3)*(Apak^4)
    return(term1+term2+term3+term4+term5)
  }
}
#' @keywords internal
#' @noRd
vmf_mexact <- function(kappa, n, p, Rbar){
  Apk  = aux_vmf_Apk(p, kappa)
  Apak = aux_vmf_Apk(p, n*Rbar*kappa)
  
  output = list()
  output$ord0 = Apk-Rbar*Apak
  output$ord1 = vmf_Apakdiff(Apk,1,p,kappa,order=1)-Rbar*vmf_Apakdiff(Apak,n*Rbar,p,kappa,order=1)
  output$ord2 = vmf_Apakdiff(Apk,1,p,kappa,order=2)-Rbar*vmf_Apakdiff(Apak,n*Rbar,p,kappa,order=2)
  output$ord3 = vmf_Apakdiff(Apk,1,p,kappa,order=3)-Rbar*vmf_Apakdiff(Apak,n*Rbar,p,kappa,order=3)
  return(output)
}

#' @keywords internal
#' @noRd
vmf_me1 <- function(x){
  n = nrow(x)
  p = ncol(x)
  kap.old = vmf_2005banerjee(x)
  Rbar = aux_vmf_Rbar(x)
  
  maxiter = 123
  epsthr  = 1e-8
  sqrteps = sqrt(.Machine$double.eps)
  if (Rbar < 1/sqrt(n)){
    return(0)
  } else {
    incval = 100000
    citer  = 1
    while (incval > epsthr){
      hh      = min(abs(kap.old),sqrteps)/2
      
      ## compute something
      mexact <- vmf_mexact(kap.old, n, p, Rbar)
      f0 <- mexact$ord0
      f1 <- mexact$ord1
      f2 <- mexact$ord2
      f3 <- mexact$ord3
      
      ## update the kappa
      term1 = f0
      term2 = f1
      
      kap.new = kap.old - term1/term2
      
      # update
      incval  = abs(kap.new-kap.old)
      citer   = citer + 1  # update the iteration
      kap.old = kap.new
      if (citer > maxiter){
        break
      }
    }
    return(kap.old)
  }
}

#' @keywords internal
#' @noRd
vmf_me2 <- function(x){
  n = nrow(x)
  p = ncol(x)
  kap.old = vmf_2005banerjee(x)
  Rbar = aux_vmf_Rbar(x)
  
  maxiter = 123
  epsthr  = 1e-8
  sqrteps = sqrt(.Machine$double.eps)
  if (Rbar < 1/sqrt(n)){
    return(0)
  } else {
    incval = 100000
    citer  = 1
    while (incval > epsthr){
      hh      = min(abs(kap.old),sqrteps)/2
      
      ## compute something
      mexact <- vmf_mexact(kap.old, n, p, Rbar)
      f0 <- mexact$ord0
      f1 <- mexact$ord1
      f2 <- mexact$ord2
      f3 <- mexact$ord3
      
      ## update the kappa
      term1 = 2*f0*f1
      term2 = 2*(f1^2) - f0*f2
      
      kap.new = kap.old - term1/term2
      
      # update
      incval  = abs(kap.new-kap.old)
      citer   = citer + 1  # update the iteration
      kap.old = kap.new
      if (citer > maxiter){
        break
      }
    }
    return(kap.old)
  }
}


# 8. numerical mle of orders 1 and 2 --------------------------------------
#' @keywords internal
#' @noRd
vmf_nmle1 <- function(dat){
  p       = ncol(x)
  Rbar    = aux_vmf_Rbar(x)
  kap.old = vmf_2005banerjee(x)  # initialization
  
  maxiter = 123
  reltol  = 1e-7
  valtol  = 12345
  citer   = 1
  sqrteps = (.Machine$double.eps)^(1/6)
  
  my.fun <- function(kappa){
    return(aux_vmf_Apk(p,kappa)-Rbar)
  }
  
  while (valtol > reltol){
    hh = min(abs(kap.old),sqrteps)/2
    
    fc = my.fun(kap.old)
    fr = my.fun(kap.old+hh)
    fl = my.fun(kap.old-hh)
    
    f0 = fc
    f1 = (fr-fl)/(2*hh)
    
    kap.new = kap.old - (f0/f1)
    valtol  = abs(kap.new-kap.old)/kap.old
    citer = citer + 1
    
    kap.old = kap.new
    if (citer >= maxiter){
      break
    }
  }
  return(kap.old)
}
#' @keywords internal
#' @noRd
vmf_nmle2 <- function(dat){
  p       = ncol(x)
  Rbar    = aux_vmf_Rbar(x)
  kap.old = vmf_2005banerjee(x)  # initialization
  
  maxiter = 123
  reltol  = 1e-7
  valtol  = 12345
  citer   = 1
  sqrteps = (.Machine$double.eps)^(1/6)
  
  my.fun <- function(kappa){
    return(aux_vmf_Apk(p,kappa)-Rbar)
  }
  
  while (valtol > reltol){
    hh = min(abs(kap.old),sqrteps)/2
    
    fc = my.fun(kap.old)
    fr = my.fun(kap.old+hh)
    fl = my.fun(kap.old-hh)
    
    f0 = fc
    f1 = (fr-fl)/(2*hh)
    f2 = (fr-(2*fc)+fl)/(hh^2)
    
    kap.new = kap.old - (2*f0*f1)/(2*(f1^2) - f0*f2)
    valtol  = abs(kap.new-kap.old)/kap.old
    citer = citer + 1
    
    kap.old = kap.new
    if (citer >= maxiter){
      break
    }
  }
  return(kap.old)
}


# 9. Christie (2014) Taylor Series ----------------------------------------
#  I'll go with 4-th order approximation.
#' @keywords internal
#' @noRd
vmf_christie <- function(x){
  ## preliminary computations
  p = ncol(x)
  Rbar = aux_vmf_Rbar(x)
  kap0 = vmf_2005banerjee(x)
  R0   = aux_vmf_Apk(p, kap0)
  
  ## order 0
  cr0 = p-kap0*((1-(R0^2))/R0)
  ## order 1
  Phi = 1/(1-cr0) - 2*(R0^2)/(1-(R0^2)) - 1
  cr1 = (1/R0)*(cr0-p)*Phi
  ## order 2
  dPhi = cr1/(1-(cr0^2)) - 4*R0/((1-(R0^2))^2)
  cr2  = (cr1*(Phi-1) + (cr0-p)*dPhi)/R0
  ## order 3
  ddPhi = 2*(cr1^2)/((1-cr0)^3) + cr2/((1-cr0)^2) - 4*(1+(3*(R0^2)))/((1-R0)^3)
  cr3   = (cr2*(Phi-2)+(2*cr1*dPhi)+((cr0-p)*ddPhi))/R0
  ## order 4
  dddPhi = (6*(cr1^3))/((1-cr0)^4) + 6*cr1*cr2/((1-cr0)^3) + cr3/((1-cr0)^2) - 48*R0*(1+(R0)^2)/((1-R0)^4)
  cr4    = (cr3*(Phi-3) + 3*cr2*dPhi + 3*cr1*ddPhi + (cr0-p)*dddPhi)/R0 # 3*cr1 part is unclear
  
  ## compute the output
  cRbar  = cr0 + cr1*(Rbar-R0) + (cr2/2)*((Rbar-R0)^2) + (cr3/6)*((Rbar-R0)^3) + (cr4/24)*((Rbar-R0)^3)
  output = Rbar*(p-cRbar)/(1-(Rbar^2))
  return(output)
}

# ## simply test
# x = rvmf(500, c(1,rep(0,7)), 10)
# mle.vmf(x, method="banerjee")
# mle.vmf(x, method="song")
# mle.vmf(x, method="sra")
# mle.vmf(x, method="tanabe")
# mle.vmf(x, method="uniroot")
# mle.vmf(x, method="nmarg1")
# mle.vmf(x, method="nmarg2")
# mle.vmf(x, method="me1")
# mle.vmf(x, method="me2")
# mle.vmf(x, method="nmle1")
# mle.vmf(x, method="nmle2")
# mle.vmf(x, method="christie")
