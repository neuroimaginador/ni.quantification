make_segmentation_malf <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("segmentation")
  info <- get_problem_info(path)

  # Templates
  templates <- read_nifti_batch_4d(info$inputs[[1]])

  labels <- read_nifti_batch_4d(info$outputs)
  original_image <- image
  mask <- sum_4d(labels)

  reference <- build_average_image(templates)

  image_to_reference <- find_transform(source_image = image, target_image = reference)
  reference_to_image <- find_transform(source_image = image_to_reference$image_norm, target_image = image)

  new_image <- map(source = image_to_reference$image_norm, target = reference, nbins = 128)

  # Apply atlas fusion algorithm to the image
  label_ids <- c(0, info$values)
  segmentation <- malf(input_image = new_image,
                       mask = mask,
                       template4D = templates,
                       labels4D = labels,
                       label_ids = label_ids,
                       patch_size = 13,
                       search_size = 5,
                       lambda = 0.7,
                       max_random_neighbours = 3,
                       kernel_width = 0,
                       return_memberships = TRUE)

  segmentation2 <- 0 * segmentation
  for (i in seq_along(label_ids)) {

    segmentation2[, , , i] <- deform_volume(V = segmentation[, , , i],
                                     Dx = reference_to_image$transform$Dx,
                                     Dy = reference_to_image$transform$Dy,
                                     Dz = reference_to_image$transform$Dz,
                                     target_dims = dim(reference), method = 3)

  }

  kernel <- gaussian_kernel(sigma = 1, size = 5)
  segmentation2 <- regularize(image = segmentation2, kernel = kernel, ncores = parallel::detectCores() - 1)
  segmentation <- defuzzify(segmentation2) - 1

  segmentation <- segmentation * (original_image > 0)

  return(segmentation)

}