#' Resample image to isotropic 1mm3 voxels
#'
#' @param img       (nifti object) Image to be resampled.
#' @param interp    Type of interpolation (0 = linear, 1 = nearest-neighbour - useful for resampling parcellations)
#'
#' @return A 3D array with the resampled image
#'
resample_to_1mm3 <- function(img, interp = 0) {
  
  require(ANTsR)
  
  if (is.nifti(img))
    
    return(ANTsR::as.array(resampleImage(as.antsImage(img, spacing = pixdim(img)[2:4]), c(1, 1, 1), 0, interpType = interp)))
  
  else
    
    return(img)
  
}