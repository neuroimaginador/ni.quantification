make_segmentation_malf <- function(image) {

  require(ni.datasets)
  require(utils4ni)

  # Prepare priors
  path <- get_dataset("segmentation")
  info <- get_problem_info(path, num_subjects = 10)

  # Templates
  templates <-  read_nifti_batch(info$inputs[[1]]) %>%
    unlist() %>%
    array(dim = c(dim(image), length(info$inputs[[1]])))

  labels <- read_nifti_batch(info$outputs) %>%
    unlist() %>%
    array(dim = c(dim(image), length(info$outputs)))

  # Apply atlas fusion algorithm to the image
  segmentation <- apply_atlas_fusion_opal_omp(input_image = image,
                                              template4D = templates,
                                              labels4D = labels,
                                              label_ids = c(0, info$values),
                                              patch_size = 5,
                                              k = 7,
                                              search_size = 5)

  segmentation <- segmentation * (image > 0)

  return(segmentation)

}