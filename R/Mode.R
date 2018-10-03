#' Computes the mode (most frequent value) of a vector or array.
#' @title Mode: Most Frequent Value in Vector
#' 
#' @author Brain Dynamics
#' 
#' @param x  A vector, array or matrix whose most frequent value is computed. Needs not be numeric.
#' 
#' @return The mode of the input vector.
#' 
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
