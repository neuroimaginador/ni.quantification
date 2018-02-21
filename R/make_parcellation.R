make_parcellation <- function(betted_image, segmentation, subcortical_parcellation) {

  require(ni.datasets)
  require(utils4ni)

  GM_mask <- segmentation
  GM_mask[GM_mask != 2] <- 0
  GM_mask[GM_mask > 0] <- 1

  filename <- get_atlas("MNI152_T1")
  brain.mni <- read_nifti_to_array(filename)
  dim(brain.mni) <- c(dim(brain.mni), 1)

  # Invert registration parameters to label in the subject native space
  filename <- get_atlas("AAL116")
  atlasMNI <- read_nifti_to_array(filename)
  dim(atlasMNI) <- c(dim(atlasMNI), 1)

  original_image <- betted_image
  betted_image <- map(source = betted_image, target = brain.mni[, , , 1], nbins = 128)

  GM_parcellation <- malf(input_image = betted_image,
                          mask = original_image > 0,
                          template4D = brain.mni,
                          labels4D = atlasMNI,
                          label_ids = c(0, 1:116),
                          kernel_width = 0)

  GM_parcellation <- GM_parcellation * GM_mask
  GM_parcellation[subcortical_parcellation > 1 & subcortical_parcellation <= 116] <- subcortical_parcellation[subcortical_parcellation > 1 & subcortical_parcellation <= 116]

  GM_parcellation <- remove_small_cc(GM_parcellation, pctg = 0.25)
  GM_parcellation <- extend_labels(pIn = GM_parcellation, maskImage = GM_mask)

  filename <- get_atlas("atlasWM")
  atlasMNI <- read_nifti_to_array(filename)

  WM_mask <- segmentation
  WM_mask[WM_mask != 3] <- 0
  WM_mask[WM_mask > 0] <- 1

  atlas20 <- atlasMNI
  atlas20[atlas20 > 20] <- 0
  dim(atlas20) <- c(dim(atlas20), 1)

  WM_parcellation <- malf(input_image = betted_image,
                          mask = original_image > 0,
                          template4D = brain.mni,
                          labels4D = atlas20,
                          label_ids = c(0, 1:20),
                          kernel_width = 0)

  WM_parcellation <- WM_parcellation * WM_mask

  extended_image <- extend_labels(pIn = GM_parcellation, maskImage = WM_mask)

  parcellation <- extended_image + 136
  parcellation[GM_parcellation > 0] <- GM_parcellation[GM_parcellation > 0]
  parcellation[WM_parcellation > 0] <- 116 + WM_parcellation[WM_parcellation > 0]
  parcellation[subcortical_parcellation > 0] <- subcortical_parcellation[subcortical_parcellation > 0]
  parcellation[segmentation <= 1] <- 0

  parcellation <- parcellation %>% as.integer() %>% array(dim = dim(betted_image))

  new_parcellation <- remove_small_cc(parcellation, pctg = 0.4)
  parcellation <- extend_labels(pIn = new_parcellation, maskImage = parcellation)

  return(parcellation)

}