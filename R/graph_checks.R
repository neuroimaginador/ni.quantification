#' Check if a matrix is square
#'
#' @param M   (matrix) The matrix to check
#'
#' @return \code{TRUE} if the matrix is square, \code{FALSE} otherwise.
#'
#' @examples
#' A <- matrix(runif(252 * 252), ncol = 252)
#' check_square_matrix(A) # TRUE
#' A <- matrix(runif(10 * 20), ncol = 10)
#' check_square_matrix(A) # FALSE
#'
check_square_matrix <- function(M) {

  # We check that the number of columns and rows is the same.
  return(dim(M)[1] == dim(M)[2])

}

#' Check if a matrix is symmetric
#'
#' @param M   (matrix) The matrix to check
#'
#' @return \code{TRUE} if the matrix is symmetric, \code{FALSE} otherwise.
#'
#' @examples
#' A <- matrix(runif(252 * 252), ncol = 252)
#' check_square_matrix(A) # FALSE
#' A <- matrix(runif(10 * 10), ncol = 10)
#' A <- A + t(A)
#' check_square_matrix(A) # TRUE
#'
#' @details It must be checked that the matrix is squared before running this.
#'
check_symmetric_matrix <- function(M) {

  # Is the matrix equal to it transposed?
  return(all(M == t(M)))

}


#' Remove Autoloops in an Adjacency Matrix
#'
#' @param A   (matrix) The adjacency matrix to remove autoloops
#'
#' @return The matrix with 0s in the diagonal.
#'
#' @examples
#' A <- matrix(runif(252 * 252), ncol = 252)
#' remove_autoloops(A)
#'
#' @details It must be checked that the matrix is squared before running this.
#'
remove_autoloops <- function(A) {

  # Remove the diagonal
  return(A - diag(diag(A)))

}

#' Binarize matrix (make all elements 0 or 1)
#'
#' @param A          (matrix) The adjacency matrix to binarize.
#' @param threshold  (numeric) The value to binarize. Elements with value above it,
#' will be 1 in the output, otherwise they will be 0.
#'
#' @return The binarized matrix.
#'
#' @examples
#' A <- matrix(runif(252 * 252), ncol = 252)
#' binarize_matrix(A, 0.25)
#'
binarize_matrix <- function(A, threshold = 0) {

  return(A > threshold)

}


#' Checks for a Correct Adjacency Matrix
#'
#' @param A (matrix) Matrix to check.
#'
#' @return \code{TRUE} or \code{FALSE} depending on whether the matrix is square
#' and symmetric.
#'
correct_adjacency_matrix <- function(A) {

  correct <- check_square_matrix(A) && check_symmetric_matrix(A)

  return(correct)

}


#' Build a Graph from an Adjacency Matrix
#'
#' @param M          (matrix) The adjacency matrix
#' @param binarize   (logical) Boolean defining if the graph should use the binarized adjacency matrix
#'
#' @return A graph object from the \code{igraph} package.
#'
#' @details First, we check that the adjacency matrix is correct. If it is, it
#' is binarized and autoloops removed and check if there is only one connected component.
#'
build_graph_from_matrix <- function(M, binarize = FALSE) {

  require(igraph)

  if (correct_adjacency_matrix(M)) {

    if (any(is.na(M))) {

      warning("The input adjacency matrix has NAs.")

    }

    if (binarize) {

      if (any(M < 0)) {

        warning("There are negative values in the adjacency matrix.")

      }

      B <- binarize_matrix(M)

    } else {

      B <- M

    }

    B <- remove_autoloops(B)

    G <- graph.adjacency(adjmatrix = B, diag = FALSE, weighted = TRUE, mode = "undirected")

    if (!is_connected(G)) {

      warning("The resulting graph has more than one connected component.")

    }

    return(G)

  } else {

    stop("Wrong adjacency matrix.")

  }


}
