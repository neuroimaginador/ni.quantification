corregister_intensive <- function(ref_image, mov_image) {

  # ref_file <- tempfile(fileext = ".nii.gz")
  # RNifti::writeNifti(ref_image, file = ref_file)
  #
  # mov_file <- tempfile(fileext = ".nii.gz")
  # RNifti::writeNifti(mov_image, file = mov_file)
  #
  require(ANTsR)

  # transformed <- antsRegistration(fixed = antsImageRead(ref_file),
  #                                 moving = antsImageRead(mov_file),
  #                                 typeofTransform = "SyN",
  #                                 verbose = TRUE)
  #
  if (is.character(ref_image)) ref_image <- antsImageRead(ref_image) %>% as.array()

  transformed <- antsRegistration(fixed = as.antsImage(ref_image),
                                  moving = as.antsImage(mov_image),
                                  typeofTransform = "AffineFast",
                                  verbose = TRUE)

  return(as.array(transformed$warpedmovout))

}