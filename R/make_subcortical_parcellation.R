make_subcortical_parcellation <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("parcellation")
  # info <- get_problem_info(path, num_subjects = 10, interactive = FALSE)
  info <- get_problem_info(problem_path = "",
                           input_path = file.path(path, "inputs"),
                           output_path = file.path(path, "outputs", "transformed"),
                           num_subjects = 10, interactive = FALSE)

  # Detect subcortical GM
  scgm_labels <- c(10, 11, 12, 13, 17, 18, 49:54)
  spinal_cord_labels <- c(16)
  ventricles_labels <- c(4, 5, 14, 15, 24, 43, 44, 72)
  all_labels <- c(0, scgm_labels, spinal_cord_labels, ventricles_labels)

  info %>% subset_problem(subset_classes = c(scgm_labels, spinal_cord_labels),
                          unify_classes = ventricles_labels,
                          use_all = TRUE)

  # Templates
  templates <- read_nifti_batch(info$inputs[[1]]) %>%
    unlist() %>%
    array(dim = c(dim(image), length(info$inputs[[1]])))

  labels <- read_nifti_batch(info$outputs) %>%
    unlist() %>% as.integer() %>%
    array(dim = c(dim(image), length(info$outputs)))


  # system.time(labels1 <- labels %>% map_ids(remap_classes = info$remap_classes))
  labels2 <- labels %>% map_ids_cpp(remap_classes = info$remap_classes)

  # labels <- labels %>% map_ids(remap_classes = info$remap_classes)

  # Apply atlas fusion algorithm to the image
  scgm <- apply_atlas_fusion_opal_omp(input_image = image,
                                      template4D = templates,
                                      labels4D = labels2,
                                      label_ids = c(0, unique(info$remap_classes$target, info$remap_classes$remaining)),
                                      k = 7,
                                      kernel_width = 0)

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