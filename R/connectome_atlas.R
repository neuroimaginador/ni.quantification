clean_matrix_list <- function(matrix_list) {

  matrix_list <- matrix_list[!sapply(matrix_list, is.null)]

  return(matrix_list)

}

#' Checks for Basic Properties of a List of Matrices
#'
#' @param matrix_list      List of matrices
#' @param same_dimensions  Logical expressing whether to check if matrices are all of the same dimensions
#' @param squared          Logical expressing whether to check if matrices are all squared
#'
#' @return A logical \code{TRUE} indicating that the list of matrices passed all checks,
#' and \code{FALSE} otherwise.
#'
check_matrix_list <- function(matrix_list, same_dimensions = TRUE, squared = TRUE) {

  # Initially, let us assume matrices meet requirements
  res <- TRUE

  n_rows <- sapply(matrix_list, nrow) %>% unlist()
  n_cols <- sapply(matrix_list, ncol)

  # Check if all matrices have the same dimensions (not necessarily squared)
  if (same_dimensions) {

    same_rows <- max(n_rows) - min(n_rows) == 0
    same_cols <- max(n_cols) - min(n_cols) == 0

    same_dims <- same_rows & same_cols

  } else {

    same_dims <- TRUE

  }

  # Check if all matrices are squared, not necessarily symmetric
  if (squared) {

    sqrd <- all(n_rows == n_cols)

  } else {

    sqrd <- TRUE

  }

  # Conjunction of all checks
  res <- res & same_dims & sqrd

}

#' Compute Statistical Descriptors for a Vector
#'
#' @param sample   (numeric vector) vector to compute statistics for.
#'
#' @return
#' A list with components: \code{Minimum}, with the minimum value in the \code{sample};
#' \code{Percent.5}, with the 5\% percentile; \code{Mean} and \code{Median}, with
#' mean and median of the vector, respectively; \code{Percent.95}, the 95th percentile;
#' \code{Maximum}, the maximum value in the \code{sample}; and \code{Std.Dev}, the
#' standard deviation of the input vector.
#'
#' @details In all cases, possible NAs are removed before performing the computation.
#'
extract_statistics <- function(sample) {

  q <- quantile(sample, probs = c(0.05, 0.95), na.rm = TRUE)

  res <- list("Minimum"    = min(sample, na.rm = TRUE),
              "Percent.5"  = q[1],
              "Mean"       = mean(sample, na.rm = TRUE),
              "Median"     = median(sample, na.rm = TRUE),
              "Percent.95" = q[2],
              "Maximum"    = max(sample, na.rm = TRUE),
              "Std.Dev"    = sd(sample, na.rm = TRUE))

  return(res)

}


#' Extract the Values for a Pair of Regions from a List of Connectivity Matrices
#'
#' @param connectivity_matrix_list   (list of squared matrices) A list with the
#' connectivity matrices
#' @param region_a                   (integer) First region index, must be lower than
#' or equal to the number of rows of the matrices
#' @param region_b                   (integer) Second region index, must be lower than
#' or equal to the number of rows of the matrices
#'
#' @return
#' A numeric vector with M_j[a,b] for all j (being M_j all connectivity matrices
#' in the list)
#'
sample_for_regions <- function(connectivity_matrix_list, region_a, region_b) {

  res <- sapply(connectivity_matrix_list, function(m) m[region_a, region_b])

  return(res)

}

#' Computation of Connectome Atlas based on DTI
#'
#' @param connectivity_matrix_list (list of squared matrices) List of connectivity
#' matrices to build atlas from.
#'
#' @return
#' From a list of connectivity matrices, build a connectome atlas, that is, a set
#' of statistical descriptors of all possible connections in the brain, given, in
#' this case, by DTI.
#' This also means that the connectivity matrices are assumed to be symmetric.
#' These statistical descriptors are used then to model the normal behaviour of
#' fiber connections in the brain of the studied population.
#'
connectome_atlas <- function(connectivity_matrix_list) {

  # Previous basic checks
  connectivity_matrix_list <- clean_matrix_list(connectivity_matrix_list)
  check_result <- check_matrix_list(connectivity_matrix_list, same_dimensions = TRUE, squared = TRUE)
  stopifnot(check_result)

  res <- list()
  nROIs <- nrow(connectivity_matrix_list[[1]])

  # For each ordered pair of regions, compute its descriptors
  for (roi_a in 1:(nROIs - 1)) {

    # Initialize as empty
    res[[roi_a]] <- list()

    for (roi_b in roi_a:nROIs) {
      res[[roi_b]] <- list()

      # Compute the given component in all possible connectivity matrices.
      sample <- sample_for_regions(connectivity_matrix_list, region_a = roi_a, region_b = roi_b)

      # Compute the statistical descriptors for the sample and store it in a "matrix-wise" list
      res[[roi_a]][[roi_b]] <- extract_statistics(sample = sample)
      res[[roi_b]][[roi_a]] <- res[[roi_a]][[roi_b]]

    }

  }

  return(res)

}