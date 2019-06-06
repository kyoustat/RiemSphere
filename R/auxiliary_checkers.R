# Auxiliary Checker Functions ---------------------------------------------
# 01. check_unitvec    : whether input is a unit vector
# 02. check_num_nonneg : whether it's a nonnegative real number



# 01. check_unitvec -------------------------------------------------------
#' @keywords internal
#' @noRd
check_unitvec <- function(x){
  cond1 = is.vector(x)
  cond2 = (abs(sum(x^2)-1) < sqrt(.Machine$double.eps))
  if (cond1&&cond2){
    return(TRUE)
  } else {
    stop(paste("'",deparse(substitute(x)),"' is not a unit vector.",sep=""))
  }
}

# 02. check_num_nonneg ----------------------------------------------------
#' @keywords internal
#' @noRd
check_num_nonneg <- function(x){
  cond1 = (length(x)==1)
  cond2 = ((is.finite(x))&&(!is.na(x))&&(x>0))
  if (cond1&&cond2){
    return(TRUE)
  } else {
    stop(paste("'",deparse(substitute(x)),"' is not a nonnegative number.",sep=""))
  }
}
