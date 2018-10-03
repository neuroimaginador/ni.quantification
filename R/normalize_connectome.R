normalize_connectome <- function(connectome) {

  # Normalize with respect to the inputs (rows)
  sums <- rowSums(connectome)

  for (row in 1:nrow(connectome)) {

    if (sums[row] > 0)
      connectome[row, ] <- connectome[row, ] / sums[row]

  }

  return(connectome)

}
