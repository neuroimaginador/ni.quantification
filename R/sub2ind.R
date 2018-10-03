#' Convert multi-dimensional indexing to linear indexing
#' @title Linear Index from Multi-dimensional Subscripts
#'
#' @author Brain Dynamics
#'
#' @param dims     The dimensions of the array.
#' @param subs     Multi-dimensional subscript to convert to linear indices.
#' @param offset   (optional, default = 1) The offset for the linear indices. In R, linear indices start at 1, whileas in other languages, the usual offset is 0.
#'
#' @return Linear indices corresponding to the multi-dimensional subscripts for an array of given dimensions.
#'
sub2ind <- function(dims, subs, offset = 1L) {

  p <- cumprod(c(1, dims))
  nd <- length(dims)
  idx <- vector(mode = "integer", length = nrow(subs))
  for (i in 1:nd) {

    idx <- idx + as.integer((subs[ , i] - 1) * p[i])

  }

  idx <- idx + offset
  return(idx)

}
