##%######################################################%##
#                                                          #
####        Volumetry and Lesion Detection Flows        ####
#                                                          #
##%######################################################%##

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

basic_volumetry_flow_syn <- function() {

  flow <- preprocessing_flow()

  flow$add(what = make_brain_extraction_syn,
           inputs = c("T1_bias_field"),
           output = "brain_mask")

  flow$add(what = function(A, B) {A * B},
           inputs = c("T1_bias_field", "brain_mask"),
           output = "betted_image")

  flow$add(what = make_segmentation_iter_sigma1,
           inputs = "betted_image",
           output = "segmentation")

  flow$add(what = count_by_ROI,
           inputs = "segmentation",
           output = "basic_volumetry")

  return(flow)

}

subcortical_volumetry_flow_syn <- function() {

  flow <- basic_volumetry_flow_syn()

  flow$add(what = make_subcortical_parcellation_syn,
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

complete_volumetry_flow_syn <- function() {

  flow <- subcortical_volumetry_flow_syn()

  flow$add(what = make_parcellation,
           inputs = c("betted_image", "enhanced_segmentation", "subcortical_parcellation"),
           output = "parcellation")

  flow$add(what = count_by_ROI,
           inputs = "parcellation",
           output = "final_volumetry")

  return(flow)

}

lesion_segmentation_flow_syn <- function() {

  flow <- complete_volumetry_flow_syn()

  flow$add(inputs = "FLAIR")

  flow$add(what = function(i) corregister_intensive(ref_image = get_atlas("FLAIR"), mov_image = i),
           inputs = "FLAIR",
           output = "FLAIR_1mm")

  flow$add(what = bias_field_correction,
           inputs = "FLAIR_1mm",
           output = "FLAIR_bias_field")

  flow$add(what = corregister,
           inputs = c("T1_mni", "FLAIR_bias_field"),
           output = "FLAIR_mni")

  flow$add(what = detect_lesions,
           inputs = c("parcellation", "FLAIR_mni"),
           output = "lesion_map")

  flow$add(what = utils4ni::count_by_ROI,
           input = "lesion_map",
           output = "lesion_volumetry")

  return(flow)

}

