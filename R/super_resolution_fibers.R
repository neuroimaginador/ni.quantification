#' Generates more detailed and smoother versions of the fibers
#'
#' @title Better Visualization of Fibers
#'
#' @author Brain Dynamics
#'
#' @param fibers A \code{fibers} object as given by package \pkg{dti}.
#'
#' @return A \code{fibers} object as the original, but each fiber has been smoothed and more detail is given.
#'
super.resolution.fibers <- function(fibers){

  require(npbd.algorithms.generic)

  super.resolution <- function(f) {

    fibra <- f[ , 1:3]
    newfibra <- bezierCurve(fibra)
    newdirect <- local.direction(newfibra)

    return(cbind(newfibra, newdirect))

  }

  Nfibers <- length(fibers@startind)
  fiberstart <- fibers@startind
  fibers_ <- fibers@fibers
  fiberlength <- diff(c(fiberstart, dim(fibers_)[1] + 1))

  inda <- fiberstart
  inde <- c(fiberstart, dim(fibers_)[1] + 1)[2:(Nfibers + 1)] - 1

  new.fiberlength <- 2 * fiberlength - 1 # Duplicate the number of points
  new.fiberstart <- c(0, cumsum(new.fiberlength))[1:length(new.fiberlength)] + 1

  indices <- c(mapply(":", inda, inde))

  new.fibers <- sapply(indices, function(idx) super.resolution(fibers_[idx, ]))
  new.f2 <- do.call(rbind, lapply(new.fibers, matrix, ncol = 6))

  obj2 <- fibers
  obj2@fibers <- new.f2
  obj2@startind <- as.integer(new.fiberstart)

  return(obj2)
}
