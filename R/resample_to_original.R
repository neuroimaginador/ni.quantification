#' Resamples an image to its original voxel size
#'
#' @param img          (nifti object) Image to be resampled
#' @param original     (nifti object) Template image with final voxel dimensions.
#' @param interp       Type of interpolation. 0 for linear, 1 for nearest-neighbor - useful for resampling anatomical parcellations.
#'
#' @return A nifti object with the same dimensions and voxel size as the \code{original}.
#'
resample_to_original <- function(img, original, interp = 0) {
  
  require(ANTsR)
  require(oro.nifti)
  
  if (is.nifti(img)) {
    
    spacing <- pixdim(img)[2:4]
    
  } else {
    
    spacing <- c(1, 1, 1)
    
  }
  
  new_img <- ANTsR::as.array(resampleImage(as.antsImage(img, spacing = spacing), 
                                           resampleParams = pixdim(original)[2:4], 
                                           0, interpType = interp))
  
  new_img <- as.nifti(from = new_img, value = original)
  
  return(new_img)
  
}
