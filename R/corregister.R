corregister <- function(ref_image, mov_image) {

  # ref_file <- tempfile(fileext = ".nii.gz")
  # RNifti::writeNifti(ref_image, file = ref_file)
  #
  # mov_file <- tempfile(fileext = ".nii.gz")
  # RNifti::writeNifti(mov_image, file = mov_file)
  #
  # require(ANTsR)
  #
  # transformed <- antsRegistration(fixed = antsImageRead(ref_file),
  #                                 moving = antsImageRead(mov_file),
  #                                 typeofTransform = "AffineFast",
  #                                 verbose = TRUE)
  # transformed <- antsRegistration(fixed = as.antsImage(ref_image),
  #                                 moving = as.antsImage(mov_image),
  #                                 typeofTransform = "SyNCC",
  #                                 verbose = TRUE)

  require(RNiftyReg)

  if (is.character(ref_image)) ref_image <- utils4ni::read_nifti(ref_image)

  ref_image <- ref_image %>% as.array() %>% oro.nifti::as.nifti()
  mov_image <- mov_image %>% as.array() %>% oro.nifti::as.nifti()

  transformed <- niftyreg.linear(source = mov_image,
                                 target = ref_image,
                                 verbose = TRUE)

  # return(as.array(transformed$warpedmovout))
  return(as.array(transformed$image))

}