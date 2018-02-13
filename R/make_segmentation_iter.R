make_segmentation_iter <- function(betted_image, sigma = 0) {

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
  segmentation_init <- (defuzzify(PVE)) * (betted_image > 0)
  PVE_init <- PVE
  params_init <- params

  diff <- 1.e10
  previous_seg <- segmentation_init

  print(params)

  while (diff > 1) {

    mean_by_ROI_df <- mean_by_ROI(labelled = previous_seg, values = betted_image)
    PVE <- segmentation(betted_image, mean_by_ROI_df)

    if (sigma > 0) {

      PVE <- regularize(PVE, kernel)

    }

    previous_seg <- defuzzify(PVE) * (betted_image > 0)
    diff <- max(abs(mean_by_ROI_df - params))
    params <- mean_by_ROI_df
    print(diff)
    print(params)

  }

  return(previous_seg)

}

make_segmentation_iter_sigma1 <- function(betted_image) make_segmentation_iter(betted_image, sigma = 1)
