#' Group of functions page title
#' 
#' Group of functions Description section
#' 
#' Group of functions Details paragraph.
#'
#' @section After Arguments and Value sections:
#' Despite its location, this actually comes after the Arguments and Value sections.
#' Also, don't need to use null, could annotate first function, and then
#' using function name as the groupBy name is more intuitive.
#' 
#' @param x a param for toBar and notToBar
#' @param y a param just for notToBar
#' @return Hard to have one return section for all functions,
#' might want to have a manual list here.
#' @name splap
NULL

#' @rdname splap
#' @export
dsplap <- function(x, mu, sigma=1, log = FALSE){ # sigma : dispersion/scale parameter
  ## PREPROCESSING
  # 1. mu : unit vector
  check_unitvec(mu)
  p = length(mu)
  # 2. lambda 
  check_num_nonneg(lambda)
  # 3. check data
  check_datamat(x)
  
  ## EVALUATION
  #   1. normalizing constant
  nconstant = dsplap.constant(sigma, p)
  #   2. case branching
  if (is.vector(x)){
    logmux = aux_log(mu, x)
    output = exp(-sqrt(sum(logmux*logmux))/sigma)
  } else {
    nx     = nrow(x)
    output = rep(0,nx)
    for (i in 1:nx){
      logmux = aux_log(mu, as.vector(x[i,]))
      output[i] = exp(-sqrt(sum(logmux*logmux))/sigma)
    }
  }
  #   3. scale by normalizing constant and RETURN
  if (log){
    return(log(output)-log(nconstant))
  } else {
    return(output/nconstant)  
  }
}

#' @keywords internal
#' @noRd
dsplap.constant <- function(sigma, D){
  myfunc <- function(r){
    return(exp(-r/sigma)*((sin(r))^(D-2)))
  }
  t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
  t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
  return(t1*t2)
}

