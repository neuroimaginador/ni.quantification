#' First and last points in each Fiber
#' 
#' @title Fiber Endpoints
#' 
#' @author Brain Dynamics
#' 
#' @param fibers A \code{fiber} object, as given by the \pkg{dti} package
#' 
#' @return A \code{2*nfibers x 6} matrix where each pair of consecutive rows (1-2, 3-4, etc.) represents the first and last points of a corresponding fiber.
#' 
fiber.endpoints <- function(fibers){
  
  Nfibers <- length(fibers@startind)
  Npoints <- dim(fibers@fibers)[1]
  
  startind <- c(fibers@startind, Npoints + 1)
  
  ST <- cbind(startind[1:Nfibers], startind[2:(Nfibers + 1)] - 1)
  
  return(fibers@fibers[ST, ])
  
}