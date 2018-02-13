bias_field_correction <- function(image) {

  require(ANTsR)
  ANTsR::as.array(n4BiasFieldCorrection(img = as.antsImage(image)))

}
