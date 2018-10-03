#' Computes local (pointwise-) direction of fibers
#'
#' @title Local Directionality
#'
#' @author Brain Dynamics
#'
#' @param coords A \code{npoints x 3} matrix of (x, y, z) coordinates of fibers.
#'
#' @return A \code{npoints x 3} matrix of (dx, dy, dz) local directions.
#'
local.direction <- function(coords){

  Npoints <- dim(coords)[1]
  difs <- diff(coords)
  difs1 <- rbind(difs[1, ], difs)
  difsend <- rbind(difs, difs[Npoints - 1, ])
  newdifs <- 0.5 * (difs1 + difsend)
  M <- apply(abs(newdifs), 1, max) + 0.001
  d <- dim(difs)[2]
  for (i in 1:d) {

    newdifs[ , i] <- newdifs[ , i] / M

  }

  return(newdifs)

}
