make_brain_extraction <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("brain_extraction")
  info <- get_problem_info(path, num_subjects = 10)

  # Templates
  templates <-  read_nifti_batch(info$inputs[[1]]) %>%
    unlist() %>%
    array(dim = c(dim(image), length(info$inputs[[1]])))

  labels <- read_nifti_batch(info$outputs) %>%
    unlist() %>%
    array(dim = c(dim(image), length(info$outputs)))

  # Apply atlas fusion algorithm to the image
  brain <- apply_atlas_fusion_opal_omp(input_image = image,
                                       template4D = templates,
                                       labels4D = labels,
                                       label_ids = c(0, info$values),
                                       patch_size = 11,
                                       search_size = 7,
                                       lambda = 0.8,
                                       k = 10)

  brain_cc <- connected_components(brain)
  counts <- count_by_ROI(brain_cc)
  brain_cc[brain_cc != which.max(counts)] <- 0

  brain <- image * (brain_cc > 0)

  return(brain)

}
