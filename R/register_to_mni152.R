register_to_mni152 <- function(image) {

  require(ni.datasets)

  filename <- get_atlas("MNI152_T1")
  tmp_file <- tempfile(fileext = ".nii.gz")

  RNifti::writeNifti(image, file = tmp_file)

  require(ANTsR)

  transformed <- antsRegistration(fixed = antsImageRead(filename),
                                  moving = antsImageRead(tmp_file),
                                  typeofTransform = "AffineFast",
                                  verbose = TRUE)

  return(ANTsR::as.array(transformed$warpedmovout))

}