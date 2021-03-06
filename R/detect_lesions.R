detect_hiperintensities <- function(flair_image, mask = flair_image > 0, pre_parcellation = NULL) {

  require(utils4ni)
  require(ANTsR)

  detect_outliers <- function(x) {

    m <- mean(x, trim = 0.05, na.rm = TRUE)

    w_up <- x[x > m]
    w_down <- x[x < m]

    q1_ <- mean(w_down, trim = 0.05, na.rm = TRUE)

    q3_ <- mean(w_up, trim = 0.05, na.rm = TRUE)

    iqr_ <- q3_ - q1_


    return(list(lower = q1_ - 1.5 * iqr_, upper = q3_ + 1.5 * iqr_))

  }

  # Normalize images
  flair_norm <- flair_image / max(as.vector(flair_image), na.rm = TRUE)

  if (!is.null(pre_parcellation)) {

    v <- as.vector(flair_norm[mask > 0 & pre_parcellation != 300])

  } else {

    v <- as.vector(flair_norm[mask > 0])

  }

  # m <- lts(v)
  #
  # w_up <- v[v > m]
  # w_down <- v[v < m]
  #
  # q1_ <- lts(w_down)
  #
  # q3_ <- lts(w_up)
  #
  # iqr_ <- q3_ - q1_

  # threshold_mild <- q3_ + 1.5 * iqr_
  # threshold_severe <- q3_ + 3 * iqr_
  #
  threshold_mild <- detect_outliers(v)$upper

  lesions <- flair_norm > threshold_mild
  lesions <- lesions * mask

  if (!is.null(pre_parcellation)) {

    lesions[pre_parcellation == 0] <- 0

    ventricles <- pre_parcellation >= 700 & pre_parcellation < 1000
    ventricles_ants <- as.antsImage(ventricles)
    ventricles_ants <- iMath(ventricles_ants, "MD") # Dilate ventricles
    ventricles <- as.array(ventricles_ants)

    brain_mask <- mask + ventricles > 0
    brain_mask <- as.array(iMath(as.antsImage(brain_mask), "MC", 7))
    brain_mask <- as.array(iMath(as.antsImage(brain_mask), "FillHoles"))

    lesions[brain_mask == 0] <- 0

  }

  lesions_ants <- as.antsImage(lesions)
  lesions_ants <- iMath(lesions_ants, "FillHoles")
  lesions_ants <- iMath(lesions_ants, "MC", 2)
  lesions_ants <- smoothImage(lesions_ants, sigma = 1)

  lesions <- as.array(lesions_ants) #> 0.5

  if (!is.null(pre_parcellation)) {

    lesions[pre_parcellation == 300] <- 0

  }

  new_threshold <- sum(flair_image * lesions, na.rm = TRUE) / sum(lesions, na.rm = TRUE)
  new_les <- flair_image > new_threshold

  lesions_thresholded <- 2 * (lesions > 0.5)
  lesions <- extend_labels(pIn = lesions_thresholded,
                        maskImage = new_les %>% as.integer() %>% array(dim = dim(new_les)))

  lesions <- lesions / 2
  # min_threshold <- quantile(flair_image[lesions > 0], probs = c(0.3))
  # lesions[flair_image > min_threshold] <- 1
  # lesions[pre_parcellation == 0] <- 0

  if (!is.null(pre_parcellation)) {

    lesions[brain_mask == 0] <- 0
    # GM_mask <- pre_parcellation
    # GM_mask[GM_mask > 116] <- 0
    # GM_mask[GM_mask > 0] <- 1
    # GM_mask <- as.array(iMath(as.antsImage(GM_mask), "MD"))
    # # GM_mask <- as_array(dilation_by_distance(convert_to_int(create(GM_mask)), radius = 0.1))
    # lesions[GM_mask > 0] <- 0
    # # lesions[pre_parcellation == 300] <- 0
    spinal_cord <- pre_parcellation * (pre_parcellation == 300)
    spinal_cord[spinal_cord > 0] <- 1
    spinal_cord_mask <- as.array(iMath(as.antsImage(spinal_cord), "MD"))
    # spinal_cord_mask <- as_array(dilation_by_distance(mask_values(create(spinal_cord), 299.5, 300.5, 1, 0), radius = 1.5))
    lesions[spinal_cord_mask > 0] <- 0

  }

  # lesions <- as.array(iMath(as.antsImage(lesions), "MO"))
  # RNiftiExtension::opening(lesions, radius = 1)
  # lesions <- remove_small_cc_ptr(lesions)
  # lesions[flair_image > threshold_severe & lesions > 0] <- 2

  return(lesions)

}

detect_lesions <- function(parcellation, FLAIR) {

  detect_hiperintensities(flair_image = FLAIR,
                          mask = parcellation > 116 & parcellation < 400,
                          pre_parcellation = parcellation)

}