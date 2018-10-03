# Input will always be a graph (igraph object).
# Number of ROIs are only for the AAL 116 atlas. Different atlases will give a different number of ROIs and, so, of nodes...

# Número de fibras por cada región anatómica: 116 valores.
# This will require the weighted graph.
number_of_fibers <- function(G) {

  require(igraph)

  fibers_per_node <- strength(G, weights = E(G)$weights)

  return(fibers_per_node)

}

# Número de fibras que conectan cada par de regiones: 6670 valores.
# This will require the weighted graph
fibers_between_ROIs <- function(G) {

  require(igraph)

  weights <- E(G)$weight
  M <- as_adj(G, sparse = FALSE)
  M[upper.tri(M)] <- 0
  M <- M - diag(diag(M))
  idx <- which(M > 0)
  M[idx] <- weights

  return(M[lower.tri(M)])

}

# Grado de cada nodo: 116 valores.
node_degree <- function(G) {

  require(igraph)

  return(degree(G))

}

# Grado de la red: 1 valor.
mean_degree <- function(G) {

  require(igraph)

  return(mean(degree(G)))

}

# Densidad de la red: 1 valor.
graph_density <- function(G) {

  require(igraph)

  return(edge_density(G))

}

#  Distancias mínimas entre pares de nodos: 6670 valores.
minimum_distances <- function(G) {

  require(igraph)

  D <- distances(G, weights = NA)

  return(D[upper.tri(D)])

}

minimum_weighted_distances <- function(G) {

  require(igraph)

  D <- distances(G, weights = NULL)

  return(D[upper.tri(D)])

}

# Distancia media de cada nodo al resto: 116 valores.
distance_from_node <- function(G) {

  if (is_connected(G)) {

    require(igraph)

    D <- distances(G, weights = NA)
    D <- apply(D, 1, mean)

  } else {

    D <- rep(Inf, vcount(G))

  }

  return(D)

}

weighted_distance_from_node <- function(G) {

  if (is_connected(G)) {

    require(igraph)

    D <- distances(G, weights = NULL)
    D <- apply(D, 1, mean)

  } else {

    D <- rep(Inf, vcount(G))

  }

  return(D)

}

# CPL: 1 valor.
characteristic_path_length <- function(G) {

  require(igraph)

  distances <- minimum_distances(G)

  return(mean(distances))

}

weighted_characteristic_path_length <- function(G) {

  require(igraph)

  distances <- minimum_weighted_distances(G)

  return(mean(distances))

}

# Eficiencia de cada nodo: 116 valores.
node_efficiency <- function(G) {

  require(igraph)

  D <- distances(G, weights = NA)

  invD <- 1 / D

  diag(invD) <- 0

  efficiency <- apply(invD, 1, sum) / (vcount(G) - 1)

  return(efficiency)

}

node_weighted_efficiency <- function(G) {

  require(igraph)

  D <- distances(G, weights = NULL)

  invD <- 1 / D

  diag(invD) <- 0

  efficiency <- apply(invD, 1, sum) / (vcount(G) - 1)

  return(efficiency)

}

# Eficiencia global: 1 valor.
average_efficiency <- function(G) {

  require(igraph)

  D <- distances(G, weights = NA)

  efficiency <- mean(1 / D[upper.tri(D)])

  return(efficiency)

}

average_weighted_efficiency <- function(G) {

  require(igraph)

  D <- distances(G, weights = NULL)

  efficiency <- mean(1 / D[upper.tri(D)])

  return(efficiency)

}

# Coeficiente de agrupamiento de cada región: 116 valores.
clustering <- function(G) {

  require(igraph)

  return(transitivity(G, "localundirected"))

}

# Coeficiente de agrupamiento global: 1 valor.
global_clustering <- function(G) {

  require(igraph)

  return(transitivity(G, "globalundirected"))

}

# CPL normalizado: 1 valor.
normalized_characteristic_path_length <- function(G, n_random_graphs = 100) {

  require(igraph)

  mean_random_cpl <- 0

  for (i in 1:n_random_graphs) {

    random_graph <- erdos.renyi.game(vcount(G), ecount(G), type = "gnm")

    mean_random_cpl <- mean_random_cpl + characteristic_path_length(random_graph)

  }

  mean_random_cpl <- mean_random_cpl / n_random_graphs

  return(characteristic_path_length(G) / mean_random_cpl)

}

# Coeficiente de agrupamiento global normalizado: 1 valor.
normalized_global_clustering <- function(G, n_random_graphs = 100) {

  require(igraph)

  mean_global_clustering <- 0

  for (i in 1:n_random_graphs) {

    random_graph <- erdos.renyi.game(vcount(G), ecount(G), type = "gnm")

    mean_global_clustering <- mean_global_clustering + global_clustering(random_graph)

  }

  mean_global_clustering <- mean_global_clustering / n_random_graphs

  return(global_clustering(G) / mean_global_clustering)

}

#### Funcional
# diámetro de la red: 1 valor
graph_diameter <- function(G) {

  require(igraph)

  if (is_connected(G)) {

    return(diameter(G, directed = FALSE, weights = NA))

  } else {

    return(Inf)

  }

}

weighted_graph_diameter <- function(G) {

  require(igraph)

  if (is_connected(G)) {

    return(diameter(G, directed = FALSE, weights = NULL))

  } else {

    return(Inf)

  }

}

# medida de centralidad por cercanía: 116 valores
node_closeness <- function(G) {

  closeness(G, mode = "all", weights = NA, normalized = TRUE)

}

weighted_node_closeness <- function(G) {

  closeness(G, mode = "all", weights = NULL, normalized = FALSE)

}

# cercanía media del grafo: 1 valor
average_closeness <- function(G) {

  nodal_closeness <- node_closeness(G)

  return(mean(nodal_closeness))

}

weighted_average_closeness <- function(G) {

  nodal_closeness <- weighted_node_closeness(G)

  return(mean(nodal_closeness))

}

# centralidad por intermediación: 116 valores
node_betweenness <- function(G) {

  betweenness(G, directed = FALSE, weights = NA, normalized = TRUE)

}

weighted_node_betweenness <- function(G) {

  betweenness(G, directed = FALSE, weights = NULL, normalized = FALSE)

}

# índice de intermediación media
average_betweenness <- function(G) {

  nodal_betweenness <- node_betweenness(G)

  return(mean(nodal_betweenness))

}

weighted_average_betweenness <- function(G) {

  nodal_betweenness <- weighted_node_betweenness(G)

  return(mean(nodal_betweenness))

}

# medida de intermediación de cada arista
betweennes_for_edges <- function(G) {

  res <- edge_betweenness(G, directed = FALSE, weights = NA)

  # Falta rellenar
  M <- as_adj(G, sparse = FALSE)
  M[upper.tri(M)] <- 0

  idx <- which(M > 0)

  M[idx] <- res

  return(M[lower.tri(M)])

}

weighted_betweennes_for_edges <- function(G) {

  res <- edge_betweenness(G, directed = FALSE, weights = NULL)

  M <- as_adj(G, sparse = FALSE)
  M[upper.tri(M)] <- 0
  idx <- which(M > 0)

  M[idx] <- res

  return(M[lower.tri(M)])

}

# SmallWorldness
smallworldness <- function(G, n_random_graphs = 100) {

  require(igraph)

  norm_clustering <- normalized_global_clustering(G)

  distances_original <- distances(G)
  D <- distances_original[upper.tri(distances_original)]

  average_path_length <- 0

  for (i in 1:n_random_graphs) {

    random_graph <- erdos.renyi.game(vcount(G), ecount(G), type = "gnm")
    D_tmp <- distances(random_graph)

    average_path_length <- average_path_length + mean(D_tmp[upper.tri(D_tmp)])

  }

  average_path_length <- average_path_length / n_random_graphs

  norm_shortest_path <- mean(D) / average_path_length

  return(norm_clustering / norm_shortest_path)

}


connectome_measures <- function(G) {

  require(igraph)

  # Nodes and edges
  N <- vcount(G)

  M <- 1 + 0 * as_adj(G, sparse = FALSE)
  M <- M - diag(diag(M))
  M[upper.tri(M)] <- 0
  idx <- which(M > 0)
  graph <- arrayInd(idx, .dim = dim(M))

  edge_list <- paste0(graph[, 2], " and ", graph[, 1])
  #
  # edge_list <- paste0(1:(N - 1), " and ", 2:N)

  params <- list()
  params_df <- data.frame(parameter = NULL, value = NULL)

  params$number_of_fibers <- number_of_fibers(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Number of fibers of ROI ", 1:N), value = params$number_of_fibers))

  params$fibers_between_ROIs <- fibers_between_ROIs(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Fibers between ROIs ", edge_list), value = params$fibers_between_ROIs))

  params$nodeDegree <- node_degree(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Degree of ROI ", 1:N), value = params$nodeDegree))

  params$meanDegree <- mean_degree(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Mean Degree of network"), value = params$meanDegree))

  params$graph_density <- graph_density(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Density of network"), value = params$graph_density))


  ##

  params$minimum_distances <- minimum_distances(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Minimum distance between ROIs ", edge_list),
                                value = params$minimum_distances))

  params$minimum_weighted_distances <- minimum_weighted_distances(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Minimum weighted distance between ROIs ", edge_list),
                                value = params$minimum_weighted_distances))

  params$distance_from_node <- distance_from_node(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Mean distance from ROI ", 1:N), value = params$distance_from_node))

  params$weighted_distance_from_node <- weighted_distance_from_node(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Mean weighted distance from ROI ", 1:N),
                                value = params$weighted_distance_from_node))

  params$characteristic_path_length <- characteristic_path_length(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Characteristic Path Length of network"),
                                value = params$characteristic_path_length))

  params$weighted_characteristic_path_length <- weighted_characteristic_path_length(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted characteristic path length of network"),
                                value = params$weighted_characteristic_path_length))

  params$node_efficiency <- node_efficiency(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Efficiency of ROI ", 1:N), value = params$node_efficiency))

  params$node_weighted_efficiency <- node_weighted_efficiency(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted Efficiency of ROI ", 1:N),
                                value = params$node_weighted_efficiency))

  params$average_efficiency <- average_efficiency(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Average efficiency of network"), value = params$average_efficiency))

  params$average_weighted_efficiency <- average_weighted_efficiency(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Average weighted efficiency of network"),
                                value = params$average_weighted_efficiency))

  params$clustering <- clustering(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Clustering Coefficient of ROI ", 1:N), value = params$clustering))

  params$global_clustering <- global_clustering(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Global Clustering coefficient of network"),
                                value = params$global_clustering))

  params$normalized_characteristic_path_length <- normalized_characteristic_path_length(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Normalized Characteristic Path Length of network"),
                                value = params$normalized_characteristic_path_length))

  params$normalized_global_clustering <- normalized_global_clustering(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Normalized Global Clustering coefficient of network"),
                                value = params$normalized_global_clustering))

  params$graph_diameter <- graph_diameter(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Diameter of network"), value = params$graph_diameter))

  params$weighted_graph_diameter <- weighted_graph_diameter(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted diameter of network"), value = params$weighted_graph_diameter))

  params$node_closeness <- node_closeness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Closeness of ROI ", 1:N), value = params$node_closeness))

  params$weighted_node_closeness <- weighted_node_closeness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted closeness of ROI ", 1:N),
                                value = params$weighted_node_closeness))

  params$node_betweenness <- node_betweenness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Betweenness of ROI ", 1:N), value = params$node_betweenness))

  params$weighted_node_betweenness <- weighted_node_betweenness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted betweenness of ROI ", 1:N),
                                value = params$weighted_node_betweenness))

  params$average_closeness <- average_closeness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Average closeness of network"), value = params$average_closeness))

  params$weighted_average_closeness <- weighted_average_closeness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted average closeness of network"),
                                value = params$weighted_average_closeness))

  params$average_betweenness <- average_betweenness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Average betweenness of network"), value = params$average_betweenness))

  params$weighted_average_betweenness <- weighted_average_betweenness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted average betweenness of network"),
                                value = params$weighted_average_betweenness))

  params$betweenness_for_edges <- betweennes_for_edges(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Betweenness for edge between ROIs ", edge_list),
                                value = params$betweenness_for_edges))

  params$weighted_betweenness_for_edges <- weighted_betweennes_for_edges(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Weighted betweenness for edge between ROIs ", edge_list),
                                value = params$weighted_betweenness_for_edges))

  params$smallworldness <- smallworldness(G)
  params_df <- rbind(params_df,
                     data.frame(parameter = paste0("Smallworldness of network"), value = params$smallworldness))

  return(list(parameters = params, parameters_df = params_df))

}




