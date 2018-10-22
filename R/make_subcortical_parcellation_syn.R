make_subcortical_parcellation_syn <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("parcellation")
  # info <- get_problem_info(path, num_subjects = 10, interactive = FALSE)
  info <- get_problem_info(problem_path = "",
                           input_path = file.path(path, "inputs"),
                           output_path = file.path(path, "outputs", "transformed"),
                           num_subjects = 10,
                           interactive = FALSE)

  # Detect subcortical GM
  scgm_labels <- c(10, 11, 12, 13, 17, 18, 49:54)
  spinal_cord_labels <- c(16)
  ventricles_labels <- c(4, 5, 14, 15, 24, 43, 44, 72)
  all_labels <- c(0, scgm_labels, spinal_cord_labels, ventricles_labels)

  info %>% subset_problem(subset_classes = c(scgm_labels, spinal_cord_labels),
                          unify_classes = ventricles_labels,
                          use_all = TRUE)

  # Templates
  templates <- read_nifti_batch_4d(info$inputs[[1]])

  labels <- read_nifti_batch_4d(info$outputs)

  labels2 <- labels %>% map_ids_cpp(remap_classes = info$remap_classes)

  original_image <- image
  mask <- sum_4d(labels2)

  reference <- build_average_image(templates)

  xtr <- antsRegistration(fixed = as.antsImage(reference),
                          moving = as.antsImage(image),
                          typeofTransform = "SyN",
                          verbose = TRUE)

  new_image <- utils4ni::map(source = as.array(xtr$warpedmovout),
                             target = reference,
                             nbins = 128)


  # Apply atlas fusion algorithm to the image
  label_ids <- c(0, unique(info$remap_classes$target, info$remap_classes$remaining))
  scgm_probs <- malf(input_image = new_image,
                     mask = mask,
                     template4D = templates,
                     labels4D = labels2,
                     label_ids = label_ids,
                     kernel_width = 0,
                     return_memberships = TRUE)

  kernel <- gaussian_kernel(sigma = 1, size = 3)
  scgm2 <- regularize(image = scgm_probs, kernel = kernel, ncores = parallel::detectCores() - 1)
  scgm <- defuzzify(scgm2) - 1

  scgm <- antsApplyTransforms(fixed = as.antsImage(image),
                              moving = as.antsImage(scgm),
                              transformlist = xtr$invtransforms,
                              interpolator = "genericLabel") %>% as.array()

  ventricles_idx <- info$remap_classes$target[info$remap_classes$source %in% ventricles_labels][1]
  scgm_idx <- info$remap_classes$target[info$remap_classes$source %in% scgm_labels]
  spinal_cord_idx <- info$remap_classes$target[info$remap_classes$source %in% spinal_cord_labels]
  remaining <- info$remap_classes$remaining

  ventricles_aal <- 700
  scgm_aal <- c(77, 71, 73, 75, 37, 41, 78, 72, 74, 76, 38, 42)
  spinal_cord_aal <- 300

  remap2aal116 <- list(source = c(ventricles_idx, scgm_idx, spinal_cord_idx, remaining),
                       target = c(ventricles_aal, scgm_aal, spinal_cord_aal, 0))

  parcellation <- scgm %>%
    map_ids_cpp(remap_classes = remap2aal116) %>%
    as.integer() %>%
    array(dim = dim(image))

  return(parcellation)

}