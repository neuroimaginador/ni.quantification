fuse_segmentation <- function(segmentation, subcortical_parcellation) {

  new_segmentation <- segmentation

  new_segmentation[segmentation > 1 & subcortical_parcellation > 0 & subcortical_parcellation <= 116] <- 2
  new_segmentation[segmentation > 1 & subcortical_parcellation == 300] <- 3
  new_segmentation[segmentation > 0 & subcortical_parcellation == 700] <- 1

  return(new_segmentation)

}