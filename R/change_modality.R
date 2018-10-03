#' Changes the modality of a T2 image (from a DWI acquisition) to match the usual grayscale of a T1.
#' @title Simulated Modality Transformation
#'
#' @author Brain Dynamics
#'
#' @param V      3D array with the T2 image.
#' @param mask   3D array with the same dimensions as \code{V}, indicating (\code{mask > 0}) which voxels are used to make the transformation.
#' @param trans  (string) Type of transformation to be done. Currently only \code{"t2t1"} (that is, from a T2 to a T1) is implemented. Any other value will result in the original image.
#'
#' @return If \code{trans = "t2t1"}, then the result is the transformed T2 to a T1. If not, the result is the original image.
#'
change.modality <- function(V, mask, trans = "t2t1"){

  res <- V

  try(
    switch(as.character(trans),
           t2t1 = {

             V2 <- histeq1d(V) * array(mask, dim = dim(mask))
             res <- V2 / V

           },
           t1t2 = {

             res <- V
             res <- max(res) + min(res) - res
             res[!mask] <- 0

           }
    )
  )

  return(res)
}
