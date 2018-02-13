resample <- function(image) {

  file <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(image, file = file)

  return(ANTsR::as.array(resampleImage(image = antsImageRead(file), resampleParams = c(1, 1, 1), interpType = 0)))

}
