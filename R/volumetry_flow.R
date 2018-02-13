preprocessing_flow <- function() {

  require(wf4ni)

  flow <- DLflow$new(name = "preprocessing", inputs = "T1")

  flow$add(what = register_to_mni152,
           inputs = "T1",
           output = "T1_mni")

  flow$add(what = bias_field_correction,
           inputs = "T1_mni",
           output = "T1_bias_field")

  return(flow)
}

basic_volumetry_flow <- function() {

  flow <- preprocessing_flow()

  flow$add(what = make_brain_extraction,
           inputs = c("T1_bias_field"),
           output = "betted_image")

  flow$add(what = make_segmentation_iter_sigma1,
           inputs = "betted_image",
           output = "segmentation")

  flow$add(what = quantify_ROIs,
           inputs = "segmentation",
           output = "basic_volumetry")

  return(flow)

}

subcortical_volumetry_flow <- function() {

  flow <- basic_volumetry_flow()

  flow$add(what = make_subcortical_parcellation,
           inputs = "betted_image",
           output = "subcortical_parcellation")

  flow$add(what = count_by_ROI,
           inputs = "subcortical_parcellation",
           output = "subcortical_volumetry")

  flow$add(what = fuse_segmentation,
           inputs = c("segmentation", "subcortical_parcellation"),
           output = "enhanced_segmentation")

  return(flow)

}

complete_volumetry_flow <- function() {

  flow <- subcortical_volumetry_flow()

  flow$add(what = make_parcellation,
           inputs = c("betted_image", "enhanced_segmentation", "subcortical_parcellation"),
           output = "parcellation")

  flow$add(what = count_by_ROI,
           inputs = "parcellation",
           output = "final_volumetry")

  return(flow)

}