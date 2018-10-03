#' The ind2sub command determines the equivalent subscript values corresponding to a single index into an array.
#' @title Subscripts from Linear Index
#'
#' @author Brain Dynamics
#'
#' @param dims   The dimensions (in vector form) of an array.
#' @param idx    Linear indices to be converted to subscript (n-dimensional).
#'
#' @return This function returns the n-dimensional indices that can be used as subscript in an array.
#'
ind2sub <- function(dims, idx) {

  p <- cumprod(c(1, dims[1:(length(dims) - 1)]))
  nout <- length(dims)
  res <- array(0, c(length(idx), nout))

  for (i in nout:1) {

    v <- (idx - 1) %% p[i] + 1
    vj <- (idx - v)/p[i] + 1
    idx <- v
    res[ , i] <- vj

  }

  return(res)
}
