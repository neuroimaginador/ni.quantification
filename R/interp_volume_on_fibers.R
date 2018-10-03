#' Interpolate a Volume in points given by a set of fibers
#'
#' @param fibers       (a dwiFiber object) Fibers where to interpolate
#' @param s0           (a nifti object) volume in the same sapce as the fibers
#' @param V            (array) Volume to interpolate
#' @param interp.mode  (character string) "nearest" for 1-NN interpolation. Other value will give trilinear interpolation.
#' @param with_mask    (boolean) if TRUE, use only positive values in V to interpolate. Otherwise, use the whole volume.
#'
#' @return A vector with the value of volume V interpolated at each of the points of the set of fibers, in the same order.
#'
interp_volume_on_fibers <- function(fibers, s0, V, interp.mode = "nearest", with_mask = TRUE) {

  fiber2 <- fiber_to_voxel(fibers, s0)
  points <- fiber2@fibers[, 1:3]

  # Dependencies
  require(RANN)
  require(npbd.algorithms.generic)

  d <- dim(V)

  if (with_mask) {

    idx <- which(V > 0)

  } else {

    idx <- 1:length(V)

  }
  subs.orig <- ind2sub(dims = d, idx = idx)

  nd <- length(d)


  if (interp.mode == "nearest") {

    res <- nn2(query = points, data = subs.orig, k = 2 ^ nd,
               treetype = "kd", searchtype = "standard")

    V.interp <- as.vector(V[idx][res$nn.idx[ , 1]])

  } else {

    require(oce)
    V.interp <- approx3d(x = 1:d[1], y = 1:d[2], z = 1:d[3],
                         f = V,
                         xout = points[ , 1],
                         yout = points[ , 2],
                         zout = points[ , 3])

  }

  return(V.interp)

}
