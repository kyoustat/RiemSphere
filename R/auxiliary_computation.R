# Auxiliary Computation Functions -----------------------------------------
# 01. aux_log       : mu (point) and x (point)
# 02. aux_exp       : mu (point) and v (tangent)
# 03. aux_intmean   : compute intrinsic mean returned as a vector
# 04. aux_dist_1toN : compute intrinsic distance for a vector vs matrix
# 05. aux_rotation  : compute a rotation matrix R from 'a' to 'b'; R%*%a = b

# 01. aux_log -------------------------------------------------------------
#' @keywords internal
#' @noRd
aux_log <- function(mu, x){
  theta = base::acos(sum(x*mu))
  if (abs(theta)<sqrt(.Machine$double.eps)){
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
  return(as.vector(rbase.mean(riemfactory(t(myx), name="sphere"))$x))
}



# 04. aux_dist_1toN -------------------------------------------------------
#' @keywords internal
#' @noRd
aux_dist_1toN <- function(x, maty){
  dist_one <- function(y){
    logxy = aux_log(x, y)
    return(sqrt(sum(logxy^2)))
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
# 05. aux_rotation  : compute a rotation matrix from 'a' to 'b'