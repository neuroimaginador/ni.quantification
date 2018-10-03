#' Extract regions connected by each fiber
#'
#' @title Regions connected by each fibers
#'
#' @author Brain Dynamics
#'
#' @param fibers A \code{fibers} object (as given by the \pkg{dti} package)
#' @param labels A 3D array with the anatomical parcellation
#' @param voxel.dims A vector of length 3 representing the extent of a voxel in \code{labels}
#' @param add.ends A boolean (default \code{TRUE}) indicating whether to expand original \code{fibers} with new endpoints reaching each connected region (to prevent roundup errors). For visualization purpose mainly.
#'
#' @note Some fibers may be removed in the process, those that do not lie within a specific distance of a known labelled region.
#'
#' @return A list with 3 components: \code{fibers}, a new \code{fibers} object with the new fibers (expanded if needed, and after removal or distant-to-regions fibers); \code{rois}, a \code{nfibers x 2} matrix giving, for each row, the indices of the two regions connected by the corresponding fiber; and \code{removed}, a vector of the (original) indices of removed fibers.
#'
label.fibers <- function(fibers, labels, voxel.dims, add.ends = TRUE){

  # Requirements
  require(RANN)

  query <- fiber.endpoints(fibers)[ , 1:3]
  Nfibers <- dim(query)[1]/2

  idx <- which(labels >= 1 & labels <= 116)
  data <- ind2sub(dims = dim(labels), idx = idx) - 1
  for (i in 1:3) {

    data[ , i] <- data[ , i] * voxel.dims[i]

  }

  bound <- 1.5 * max(voxel.dims)

  res <- nn2(data = data, query = query, k = 1, treetype = "kd", searchtype = "standard")
  idx.orig <- which(res$nn.dists[ , 1] > bound)
  idx2 <- idx.orig
  idx2[idx2 > Nfibers] <- idx2[idx2 > Nfibers] - Nfibers
  idx2 <- unique(idx2)
  fibers2 <- remove.fibers(fibers, idx2)

  query <- fiber.endpoints(fibers2)[ , 1:3]
  Nfibers <- dim(query)[1]/2
  res <- nn2(data = data, query = query, k = 10, treetype = "kd", searchtype = "standard")
  print(max(res$nn.dists[ , 1]))
  v <- matrix(labels[idx[res$nn.idx]], nrow = nrow(res$nn.idx))
  rois <- matrix(apply(v, 1, Mode), nrow = Nfibers)

  if (add.ends) { # For visualization purposes, add new endpoints to every fiber

    fiberstart <- fibers2@startind
    fibers_ <- fibers2@fibers
    fiberlength <- diff(c(fiberstart, dim(fibers_)[1] + 1))

    inda <- fiberstart
    inde <- c(fiberstart, dim(fibers_)[1] + 1)[2:(Nfibers + 1)] - 1

    new.fiberlength <- fiberlength + 2 # Add two points for every tract
    new.fiberstart <- c(0, cumsum(new.fiberlength))[1:length(new.fiberlength)] + 1

    new.points <- cbind(data[res$nn.idx[ , 1], ], 1, 1, 1)


    ids <- c(1:dim(fibers_)[1], inda - 0.4, inde + 0.4)
    all.points <- rbind(fibers_, new.points)
    fibers_ <- all.points[order(ids), ]
    fibers2@fibers <- fibers_
    fibers2@startind <- as.integer(new.fiberstart)

  }

  return(list(fibers = fibers2, rois = rois, removed = idx2))

}
