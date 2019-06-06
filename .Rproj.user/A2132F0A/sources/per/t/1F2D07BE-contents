# Auxiliary Computation Functions -----------------------------------------
# 01. aux_log : mu (point) and x (point)
# 02. aux_exp : mu (point) and v (tangent)



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
aux_exp <- function(mu, v){
  norm.v = sqrt(sum(v^2))
  if (norm.v < sqrt(.Machine$double.eps)){
    output = mu*cos(norm.v) + v
  } else {
    output = mu*cos(norm.v) + (sin(norm.v)/norm.v)*v
  }
  return(output)
}


