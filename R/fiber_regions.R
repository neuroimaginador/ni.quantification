#' Computes the indices of brain regions intersected by group of fibers
#'
#' @title Brain Regions per fiber.
#'
#' @author Brain Dynamics
#'
#' @param fibers A \code{fibers} object as given by package \pkg{dti}.
#' @param labels A \code{nifti} object (package \pkg{oro.nifti}) with the anatomical parcellation.
#' @param s0     (optional) The image in whose space the points are given, an object of class \code{nifti}.
#'
#' @note Both the anatomical parcellation and fibers must be in the same coordinate space, given by \code{s0}. If \code{s0} is not present, its default value is using the cooridnate space given by the \code{nifti} object \code{labels}.
#'
#' @return This function returns a list with one component for each fiber, listing the indices (as teaken from the anatomical parcellation) of the regions crossed by the given fiber.
#'
fiber.regions <- function(fibers, labels, s0 = labels){

  elapsed <- system.time(rounded <- round(fastTransformWorldToVoxel(points = fibers@fibers[, 1:3],
                                         pixdim(s0)[2:4])))[3]

  rounded[which(rounded[, 1] <= 0), 1] <- 1
  rounded[which(rounded[, 2] <= 0), 2] <- 1
  rounded[which(rounded[, 3] <= 0), 3] <- 1

  elapsed <- system.time(sub2inded <- sub2ind(dims = dim(labels),
                                              subs = rounded))[3]

  elapsed <- system.time(WMpoints <- labels[sub2inded])[3]

  Nfibers <- length(fibers@startind)
  fiberstart <- fibers@startind
  fibers_ <- fibers@fibers
  fiberlength <- diff(c(fiberstart, dim(fibers_)[1] + 1))

  inda <- fiberstart
  inde <- c(fiberstart, dim(fibers_)[1] + 1)[2:(Nfibers + 1)] - 1

  indices <- c(mapply(":", inda, inde))

  elapsed <- system.time(WM_ROIs <- sapply(indices, function(idx) unique(WMpoints[idx])))[3]
  elapsed <- system.time(WM_ROIs <- sapply(WM_ROIs, function(x) x[x > 0]))[3]

  return(WM_ROIs)
}