make_segmentation <- function(betted_image, sigma = 0) {

  require(ni.datasets)
  require(utils4ni)

  # Recursive k-means to obtain class means
  rk_seg <- rkmeans(imag = betted_image, mask = betted_image > 0)
  params <- base::sort(rk_seg$means)

  # Partial volume estimation
  PVE <- segmentation(betted_image, params)

  if (sigma > 0) {

    kernel <- gaussian_kernel(sigma = sigma, dim = 3, size = 3)
    PVE <- regularize(PVE, kernel)

  }

  # Actual segmentation
  segmentation <- (defuzzify(PVE)) * (betted_image > 0)

  return(segmentation)

}