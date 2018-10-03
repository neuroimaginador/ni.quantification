#' Extract fibers crossing a given region
#'
#' @title Extract fibers passing through a brain region
#'
#' @author Brain Dynamics
#'
#' @param fibers    A \code{fibers} object as given by package \pkg{dti}.
#' @param labels    A \code{nifti} object with the anatomical parcellation.
#' @param roi.index The index of the region to filter the fibers.
#'
#' @return A \code{fibers} object with just fibers crossing the specific region.
#'
select.fibers.ROI <- function(fibers, labels, roi.index){

  Nfibers <- length(fibers@startind)
  idx <- which(labels != roi.index)
  idx[idx > Nfibers] <- idx[idx > Nfibers] - Nfibers
  idx <- unique(idx)

  f2 <- remove.fibers(fibers, idx)

  return(f2)

}
