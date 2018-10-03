#' RGL 3D Plot of the Connectome over a Brain Image
#'
#' @param brain_parcellation     (3D array, \code{nifti} or \code{antsImage} object) Parcellation, where each label
#' must correspond to the associated node in the connectome given by \code{graph_representation}.
#' @param graph_representation   (2D square matrix, or \code{igraph} object) Graph representing the connectome, with as many nodes as labels in the parcellation, and in the same ordering.
#'
#' @return The \code{rgl} scene.
#'
plot_connectome_on_brain <- function(brain_parcellation, graph_representation) {

  # Dependencies

  require(ANTsR)
  require(igraph)
  require(oro.nifti)
  require(rgl)

  # Transform the graph into its adjacency matrix
  if (is.igraph(graph_representation)) {

    A <- as_adj(graph_representation)

  } else {

    A <- as.matrix(graph_representation)

  }

  # Compute the centroids of each region in the parcellation, to be used in the 3D plot.
  centroids <- label_centroids(img = brain_parcellation, convex = FALSE, physical = TRUE)$vertices

  # Render brain surface in RGL, with transparency
  # Maybe we need to transform the data type beforehand
  if (is.nifti(brain_parcellation)) {

    P <- as.antsImage(brain_parcellation, spacing = pixdim(brain_parcellation)[2:4])

  } else {

    P <- as.antsImage(brain_parcellation)

  }

  brain <- renderImageLabels(labelsimg = thresholdImage(P, 0, 0), smoothsval = 1, alphasurf = 0.2, col = "lightgray")

  # Overlay the connectome on the brain silhouette
  plotBasicNetwork(centroids, brain, weights = A)

  # Store the whole scene in a variable to be returned
  s <- scene3d()

  return(s)

}
