make_subcortical_parcellation <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("parcellation")
  # info <- get_problem_info(path, num_subjects = 10, interactive = FALSE)
  info <- get_problem_info(problem_path = "",
                           input_path = file.path(path, "inputs"),
                           output_path = file.path(path, "outputs", "transformed"),
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

  image_to_reference <- find_transform(source_image = image, target_image = reference, max_search = 11, resolution = 25, inner_iter = 5, max_iter = 10)
  reference_to_image <- find_transform(source_image = image_to_reference$image_norm, target_image = image, max_search = 11, resolution = 25, inner_iter = 5, max_iter = 10)

  new_image <- map(source = image_to_reference$image_norm, target = reference)

  # original_image <- image
  # image <- map(source = image, target = templates[, , , 1], nbins = 128)

  # Apply atlas fusion algorithm to the image
  label_ids <- c(0, unique(info$remap_classes$target, info$remap_classes$remaining))
  scgm_probs <- malf(input_image = new_image,
                     mask = mask,
                     template4D = templates,
                     labels4D = labels2,
                     label_ids = label_ids,
                     kernel_width = 0,
                     return_memberships = TRUE)

  scgm2 <- 0 * scgm_probs
  for (i in seq_along(label_ids)) {

    scgm2[, , , i] <- deform_volume(V = scgm_probs[, , , i],
                                    Dx = reference_to_image$transform$Dx,
                                    Dy = reference_to_image$transform$Dy,
                                    Dz = reference_to_image$transform$Dz,
                                    target_dims = dim(reference), method = 3)

  }

  kernel <- gaussian_kernel(sigma = 1, size = 5)
  scgm2 <- regularize(image = scgm2, kernel = kernel, ncores = parallel::detectCores() - 1)
  scgm <- defuzzify(scgm2) - 1

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