lesion_segmentation_flow <- function() {

  flow <- complete_volumetry_flow()

  flow$add(inputs = "FLAIR")

  flow$add(what = function(i) corregister(ref_image = get_atlas("FLAIR"), mov_image = i),
           inputs = "FLAIR",
           output = "FLAIR_1mm")

  flow$add(what = corregister,
           inputs = c("T1_mni", "FLAIR_1mm"),
           output = "FLAIR_mni")

  flow$add(what = detect_lesions,
           inputs = c("parcellation", "FLAIR_mni"),
           output = "lesion_map")

  flow$add(what = utils4ni::count_by_ROI,
           input = "lesion_map",
           output = "lesion_volumetry")

  return(flow)

}
