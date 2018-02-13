detect_hiperintensities <- function(flair_image, mask = flair_image > 0, pre_parcellation = NULL) {

  require(utils4ni)
  require(ANTsR)

  v <- as.vector(flair_image[mask > 0 & pre_parcellation != 300])

  m <- lts(v)

  w_up <- v[v > m]
  w_down <- v[v < m]

  q1_ <- lts(w_down)

  q3_ <- lts(w_up)

  iqr_ <- q3_ - q1_

  threshold_mild <- q3_ + 1.5 * iqr_
  threshold_severe <- q3_ + 3 * iqr_

  lesions <- flair_image > threshold_mild
  lesions[pre_parcellation == 0] <- 0

  if (!is.null(pre_parcellation)) {

    ventricles <- pre_parcellation >= 700 & pre_parcellation < 1000
    ventricles_ants <- as.antsImage(ventricles)
    ventricles_ants <- iMath(ventricles_ants, "MD") # Dilate ventricles
    ventricles <- ANTsR::as.array(ventricles_ants)

    brain_mask <- mask + ventricles > 0
    brain_mask <- ANTsR::as.array(iMath(as.antsImage(brain_mask), "MC", 7))
    brain_mask <- ANTsR::as.array(iMath(as.antsImage(brain_mask), "FillHoles"))

    lesions[brain_mask == 0] <- 0

  }

  lesions_ants <- as.antsImage(lesions)
  lesions_ants <- iMath(lesions_ants, "FillHoles")
  lesions_ants <- iMath(lesions_ants, "MC", 2)
  lesions_ants <- smoothImage(lesions_ants, sigma = 1)

  lesions <- ANTsR::as.array(lesions_ants) #> 0.5

  if (!is.null(pre_parcellation)) {

    lesions[pre_parcellation == 300] <- 0

  }

  new_threshold <- sum(flair_image * lesions) / sum(lesions)
  new_les <- flair_image > new_threshold

  lesions_thresholded <- 2 * (lesions > 0.5)
  lesions <- extend_labels(pIn = lesions_thresholded,
                        maskImage = new_les %>% as.integer() %>% array(dim = dim(new_les)))

  lesions <- lesions / 2
  # min_threshold <- quantile(flair_image[lesions > 0], probs = c(0.3))
  # lesions[flair_image > min_threshold] <- 1
  # lesions[pre_parcellation == 0] <- 0

  lesions[brain_mask == 0] <- 0
  # GM_mask <- pre_parcellation
  # GM_mask[GM_mask > 116] <- 0
  # GM_mask[GM_mask > 0] <- 1
  # GM_mask <- as_array(dilation_by_distance(convert_to_int(create(GM_mask)), radius = 0.1))
  # lesions[GM_mask > 0] <- 0
  # lesions[pre_parcellation == 300] <- 0
  spinal_cord <- pre_parcellation * (pre_parcellation == 300)
  spinal_cord[spinal_cord > 0] <- 1
  spinal_cord_mask <- ANTsR::as.array(iMath(as.antsImage(spinal_cord, "MD")))
  # spinal_cord_mask <- as_array(dilation_by_distance(mask_values(create(spinal_cord), 299.5, 300.5, 1, 0), radius = 1.5))
  lesions[spinal_cord_mask > 0] <- 0

  lesions <- ANTsR::as.array(iMath(as.antsImage(lesions, "MO")))
  # RNiftiExtension::opening(lesions, radius = 1)
  # lesions <- remove_small_cc_ptr(lesions)
  lesions[flair_image > threshold_severe & lesions > 0] <- 2

  return(lesions)

}

detect_lesions <- function(parcellation, FLAIR) {

  detect_hiperintensities(flair_image = FLAIR,
                          mask = parcellation > 116 & parcellation < 400,
                          pre_parcellation = parcellation)

}