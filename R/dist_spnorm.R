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
#' @name spnorm
NULL

#' @rdname spnorm
#' @export
dspnorm <- function(x, mu, lambda=1, log = FALSE){
  ## PREPROCESSING
  # 1. mu : unit vector
  check_unitvec(mu)
  p = length(mu)
  # 2. lambda 
  check_num_nonneg(lambda)
  # 3. check data
  if (!check_datamat(x)){
    stop("* asdf")
  }
  
  ## EVALUATION
  #   1. normalizing constant
  nconstant = dspnorm.constant(lambda, p)
  #   2. case branching
  if (is.vector(x)){
    logmux = aux_log(mu, x)
    output = exp(-(lambda/2)*sum(logmux*logmux))
  } else {
    dvec   = as.vector(cppdist_int_1toN(mu, x));
    output = exp((-lambda/2)*(dvec^2))
    
    # nx     = nrow(x)
    # output = rep(0,nx)
    # for (i in 1:nx){
    #   logmux = aux_log(mu, as.vector(x[i,]))
    #   output[i] = exp(-(lambda/2)*sum(logmux*logmux))
    # }
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
dspnorm.constant <- function(lbd, D){ # lbd : lambda / D : dimension
  myfunc <- function(r){
    return(exp(-lbd*(r^2)/2)*((sin(r))^(D-2)))
  }
  t1 = 2*(pi^((D-1)/2))/gamma((D-1)/2) # one possible source of error
  t2 = stats::integrate(myfunc, lower=0, upper=pi, rel.tol=sqrt(.Machine$double.eps))$value
  return(t1*t2)
}




#   -----------------------------------------------------------------------

#' @rdname spnorm
#' @param z a param just for theQuestion.
#'
#' Could put here, but easiest if all params are
#' described in the same place, with page documentation, so none get repeated.
#' @export
rspnorm <- function(n, mu, lambda=1){
  ## PREPROCESSING
  # 1. mu : unit vector
  check_unitvec(mu)
  D = length(mu)
  # 2. lambda 
  check_num_nonneg(lambda)
  
  ## ITERATE or RANDOM
  if (lambda==0){
    output = rvmf_uniform(n,mu,k=0)
  } else {
    #   1. tangent vectors
    vectors = array(0,c(n,D))
    for (i in 1:n){
      vectors[i,] = rspnorm.single(mu, lambda, D)
    }
    #   2. exponential mapping
    output = array(0,c(n,D))
    for (i in 1:n){
      output[i,] = aux_exp(mu, as.vector(vectors[i,]))
    }  
  }
  
  ## RETURN
  if (n==1){
    return(as.vector(output))
  } else {
    return(output)  
  }
}


# . -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
rspnorm.single <- function(mu, lambda, D){
  status = FALSE
  sqrts  = sqrt(1/(lambda + ((D-2)/pi))) # from the box : This is the correct one.
  #sqrts  = sqrt((lambda + ((D-2)/pi)))  # from the text
  while (status==FALSE){
    v = stats::rnorm(D,mean = 0, sd=sqrts)
    v = v-mu*sum(mu*v)       ## this part is something missed from the paper.
    v.norm = sqrt(sum(v^2))
    if (v.norm <= pi){
      if (v.norm <= sqrt(.Machine$double.eps)){
        r1 = exp(-(lambda/2)*(v.norm^2))  
      } else {
        r1 = exp(-(lambda/2)*(v.norm^2))*((sin(v.norm)/v.norm)^(D-2))
      }
      r2 = exp(-((v.norm^2)/2)*(lambda+((D-2)/pi)))
      r  = r1/r2
      u  = stats::runif(1)
      if (u <= r){
        status = TRUE
      }
    }
  }
  return(v)
}

# multiple function documentation
# https://gist.github.com/jefferys/b79fe87314b0dc72fec9
