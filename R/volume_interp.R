#' Interpolates a volume at given points.
#' @title Volume Interpolation
#'
#' @author Brain Dynamics
#'
#' @param V            Volume to be interpolated.
#' @param points       \code{N x 3} matrix of points where to interpolate.
#' @param voxel.dims   Vector of voxel dimensions (in mm.). Defaults to \code{(1, 1, 1)}.
#' @param interp.mode  Character string specifying the type of interpolation. If \code{"nearest"} (the default), the interpolation uses the nearest neighbour technique. If not, uses trilinear interpolation.
#'
#' @return The values of V interpolated at the specified points.
#'
volume.interp <- function(V, points, voxel.dims = c(1,1,1), interp.mode = "nearest") {

  # Requirements
  require(RANN)

  d <- dim(V)

  subs.orig <- ind2sub(dims = d, idx = 1:length(V))

  subs.orig[ , 1] <- (subs.orig[ , 1] - 1) * voxel.dims[1]
  subs.orig[ , 2] <- (subs.orig[ , 2] - 1) * voxel.dims[2]
  subs.orig[ , 3] <- (subs.orig[ , 3] - 1) * voxel.dims[3]

  nd <- length(d)


  if (interp.mode == "nearest") {

    res <- nn2(query = points, data = subs.orig, k = 2 ^ nd,
               treetype = "kd", searchtype = "standard")

    V.interp <- as.vector(V[res$nn.idx[ , 1]])

  } else {

    require(oce)
    V.interp <- approx3d(x = (0:(d[1] - 1)) * voxel.dims[1],
                         y = (0:(d[2] - 1)) * voxel.dims[2],
                         z = (0:(d[3] - 1)) * voxel.dims[3],
                         f = V,
                         xout = points[ , 1],
                         yout = points[ , 2],
                         zout = points[ , 3])

  }

  return(V.interp)
}
