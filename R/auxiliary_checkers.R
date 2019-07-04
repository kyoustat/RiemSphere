# Auxiliary Checker Functions ---------------------------------------------
# 01. check_unitvec    : whether input is a unit vector
# 02. check_num_nonneg : whether it's a nonnegative real number
# 03. check_datamat    : whether a given input is a vector or matrix of data on sphere



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
  cond2 = ((all(is.finite(x)))&&(!any(is.na(x)))&&(all(x>=0)))
  if (cond1&&cond2){
    return(TRUE)
  } else {
    stop(paste("'",deparse(substitute(x)),"' is not a nonnegative number.",sep=""))
  }
}


# 03. check_datamat -------------------------------------------------------
#' @keywords internal
#' @noRd
check_datamat <- function(x){
  single_checker <- function(vec){
    if (abs(sum(vec^2)-1)>sqrt(.Machine$double.eps)){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
  if (is.vector(x)){
    return(single_checker(x))
  } else if (is.matrix(x)){
    allvec = apply(x, 1, single_checker)
    if (all(allvec==TRUE)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}

# 03. check_datamat    : whether a given input is a vector or matrix of data on sphere