process_lesion_map <- function(lesion_map) {

  cc <- utils4ni::connected_components(lesion_map > 0)
  volumes <- utils4ni::count_by_ROI(cc)

  idx_to_remove <- which(volumes <= 6)
  cc[cc %in% idx_to_remove] <- 0
  lesion_map[cc %in% idx_to_remove] <- 0

  lesion_count <- length(unique(as.vector(cc[cc > 0])))
  lesion_volume <- sum(volumes[volumes > 0]) - sum(volumes[idx_to_remove])

  return(list(count = lesion_count, volume = lesion_volume, map = lesion_map))

}

lesion_evolution <- function(lesion_map_t0, lesion_map_t1) {

  # Difference (t1 - t0). Negative values are not considered.
  # Must return lesion count and volume in both timestamps.
  # Small lesions < 3mm^3 are not considered.

  `%<-%` <- zeallot::`%<-%`
  c(count_t0, volume_t0, cleaned_t0) %<-% (lesion_map_t0 %>% process_lesion_map())
  c(count_t1, volume_t1, cleaned_t1) %<-% (lesion_map_t1 %>% process_lesion_map())

  diff_lesion_map <- (lesion_map_t1 > 0) - (lesion_map_t0 > 0)
  diff_lesion_map[diff_lesion_map < 0] <- 0
  kernel <- utils4ni::gaussian_kernel()
  diff_lesion_map <- utils4ni::regularize(diff_lesion_map, kernel = kernel)
  diff_lesion_map <- diff_lesion_map > 0.85

  c(count_new_or_enlarged,
    volume_new_or_enlarged,
    cleaned_diff) %<-% (diff_lesion_map %>% process_lesion_map())

  return(list(count_t0 = count_t0,
              volume_t0 = volume_t0,
              new_lesion_map_t0 = cleaned_t0,
              count_t1 = count_t1,
              volume_t1 = volume_t1,
              new_lesion_map_t1 = cleaned_t1,
              count_new_or_enlarged = count_new_or_enlarged,
              volume_new_or_enlarged = volume_new_or_enlarged,
              lesion_map_diff = cleaned_diff))

}

