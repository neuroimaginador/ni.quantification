make_brain_extraction_syn <- function(image) {

  require(ni.datasets)
  require(utils4ni)
  require(ANTsR)

  # Prepare priors
  path <- get_dataset("brain_extraction")
  info <- get_problem_info(path, num_subjects = 20)

  # Templates
  templates <- read_nifti_batch_4d(info$inputs[[1]])
  labels <- read_nifti_batch_4d(info$outputs)

  reference <- build_average_image(templates)
  mask <- build_average_image(labels)

  xtr <- antsRegistration(fixed = as.antsImage(reference),
                          moving = as.antsImage(image),
                          typeofTransform = "SyN",
                          verbose = TRUE)

  new_image <- utils4ni::map(source = as.array(xtr$warpedmovout),
                             target = reference,
                             nbins = 128)

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

  kernel <- gaussian_kernel(sigma = 1, size = 5)
  brain2 <- regularize(image = brain, kernel = kernel, ncores = parallel::detectCores() - 1)
  brain <- defuzzify(brain2) - 1

  brain2 <- antsApplyTransforms(fixed = as.antsImage(image),
                                moving = as.antsImage(brain),
                                transformlist = xtr$invtransforms,
                                interpolator = "genericLabel") %>% as.array()

  brain_cc <- connected_components(brain2)
  counts <- count_by_ROI(brain_cc)
  brain_cc[brain_cc != which.max(counts)] <- 0

  return(brain_cc > 0)

}
