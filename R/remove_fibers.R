#' Removes fibers given by indices
#'
#' @title Remove fibers
#'
#' @author Brain Dynamics
#'
#' @param fibers A \code{fibers} object as given by package \pkg{dti}.
#' @param idx    A vector indicating which fibers to remove.
#'
#' @return A \code{fibers} object as the original, but after removal of the indicated fibers.
#'
remove.fibers <- function(fibers, idx) {

  fibers_ <- fibers@fibers
  fiberstart <- fibers@startind
  fiberlength <- diff(c(fiberstart, dim(fibers_)[1] + 1))

  remove <- (1:length(fiberstart))[idx]
  if (length(remove) > 0) {

    inda <- fiberstart[remove]
    inde <- c(fiberstart, dim(fibers_)[1] + 1)[remove + 1] - 1
    fibers_ <- fibers_[-c(mapply(":", inda, inde), recursive = TRUE), ]
    fiberlength <- fiberlength[-remove]
    fiberstart <- c(0, cumsum(fiberlength))[1:length(fiberlength)] + 1
    
  }

  res <- fibers
  res@fibers <- fibers_
  res@startind <- as.integer(fiberstart)

  return(res)

}
