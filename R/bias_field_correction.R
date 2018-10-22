bias_field_correction <- function(image) {

  image <- image %>% as.array()

  require(ANTsR)
  img <- abpN4(img = as.antsImage(image))

  return(img %>% as.array())

}
