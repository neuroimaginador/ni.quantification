make_brain_extraction <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("brain_extraction")
  info <- get_problem_info(path, num_subjects = 20)

  # Templates
  templates <- read_nifti_batch_4d(info$inputs[[1]])
  labels <- read_nifti_batch_4d(info$outputs)

  reference <- build_average_image(templates)
  mask <- build_average_image(labels)

  image_to_reference <- find_transform(source_image = image, target_image = reference,
                                       inner_iter = 1,
                                       max_search = 9,
                                       tol = 1.e-7,
                                       resolution = 25,
                                       scheme = 0)
  reference_to_image <- find_transform(source_image = image_to_reference$image, target_image = image)

  new_image <- map(source = image_to_reference$image, target = reference, nbins = 128)

  # Apply atlas fusion algorithm to the image
  label_ids <- c(0, info$values)
  brain <- malf(input_image = new_image,
                # mask = mask,
                template4D = templates,
                labels4D = labels,
                label_ids = label_ids,
                patch_size = 11,
                search_size = 7,
                lambda = 0.7,
                max_random_neighbours = 3,
                kernel_width = 0,
                return_memberships = TRUE)

  brain2 <- 0 * brain
  for (i in seq_along(label_ids)) {

    brain2[, , , i] <- deform_volume(V = brain[, , , i],
                                     Dx = reference_to_image$transform$Dx,
                                     Dy = reference_to_image$transform$Dy,
                                     Dz = reference_to_image$transform$Dz,
                                     target_dims = dim(reference), method = 3)

  }

  kernel <- gaussian_kernel(sigma = 1, size = 5)
  brain2 <- regularize(image = brain2, kernel = kernel, ncores = parallel::detectCores() - 1)
  brain <- defuzzify(brain2) - 1

  brain_cc <- connected_components(brain)
  counts <- count_by_ROI(brain_cc)
  brain_cc[brain_cc != which.max(counts)] <- 0

  return(brain_cc > 0)

}
