corregister <- function(ref_image, mov_image) {

  ref_file <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(ref_image, file = ref_file)

  mov_file <- tempfile(fileext = ".nii.gz")
  RNifti::writeNifti(mov_image, file = mov_file)

  require(ANTsR)

  transformed <- antsRegistration(fixed = antsImageRead(ref_file),
                                  moving = antsImageRead(mov_file),
                                  typeofTransform = "AffineFast",
                                  verbose = TRUE)

  return(ANTsR::as.array(transformed$warpedmovout))

}