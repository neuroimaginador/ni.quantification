safe_require <- function(pkg = NULL, username = NULL, package_name = NULL) {

  silent_require <- function(pkg_name) suppressMessages(base::require(package = pkg_name,
                                                                      character.only = TRUE))

  require(devtools)

  if (is.null(pkg) && (is.null(username) || is.null(package_name))) {

    stop("Insufficient data to load package.")

  }

  if (is.null(pkg)) {

    pkg <- paste0(username, "/", package_name)

  }

  # Determine if the package is from github (there is a "/" in the name)
  # In this case, extract user name and package name
  if (grepl(x = pkg, pattern = "/", fixed = TRUE)) {

    from_github <- TRUE
    attributes <- regmatches(x = pkg, regexpr(text = pkg, pattern = "/", fixed = TRUE), invert = TRUE)[[1]]
    username <- attributes[1]
    package_name <- attributes[2]

  } else {

    from_github <- FALSE
    package_name <- pkg

  }

  if (!silent_require(package_name)) {

    cat("Package ", package_name, " not installed yet. Attempting installation.\n")

    if (from_github) {

      suppressMessages(devtools::install_github(repo = pkg, upgrade_dependencies = FALSE))

    } else {

      install.packages(pkg)

    }

    silent_require(package_name)

  }

  return(invisible())

}

plot_adjacency_matrix <- function(mat, diag_off = TRUE, reorder = TRUE, binarize = FALSE) {

  M <- nrow(mat)

  if (diag_off) mat <- mat - diag(diag(mat))

  if (reorder) {

    mapping <- read.csv(system.file("atlases", "ROImapping_v2.csv", package = "npbd.atlases"),
                        header = TRUE, stringsAsFactors = FALSE)

    LH <- which(mapping[["a3"]] == "left_hemisphere")
    RH <- which(mapping[["a3"]] == "right_hemisphere")
    CEREB <- which(mapping[["a3"]] == "cerebellum")

    LH <- LH[LH <= M]
    RH <- RH[RH <= M]
    CEREB <- CEREB[CEREB <= M]

    order <- c(LH, RH, CEREB)

  } else {

    order <- 1:M

  }

  safe_require("jonclayden/shades")

  # colors <- c("white", "red") #unclass(gradient(c("white", "red"), 2, space = "Lab"))

  if (min(mat) < 0)
    colors <- unclass(gradient(c("darkblue", "darkgreen", "yellow", "orange", "red"), 100, space = "Lab"))
  else
    colors <- unclass(gradient(c("darkblue", "yellow", "red"), 100, space = "Lab"))


  old_par <- par(mar = c(0, 0, 0, 0), bg = "gray")

  if (binarize)
    image(1:M, 1:M, sign(mat[order, rev(order)]), col = colors, asp = 1, xaxs = "i")
  else
    image(1:M, 1:M, mat[order, rev(order)], col = colors, asp = 1, xaxs = "i")

  if (reorder) {

    abline(h = M - length(LH) + 0.5, v = length(LH) + 0.5, col = "black")
    abline(h = M - (length(LH) + length(RH)) + 0.5, v = length(LH) + length(RH) + 0.5, col = "black")

  }

  par(old_par)

  invisible()

}

#####################################################################
##                                                                 ##
## - Plots the connectome reading it from the CSV adjacency        ##
##   matrix.                                                       ##
##                                                                 ##
## @Params: conn.csv, the CSV file with the adjacency matrix       ##
## @Params: diag0, boolean that, if it's true, sets the diagonal   ##
##          to 0                                                   ##
##                                                                 ##
#####################################################################

plot.connectome <- function (conn.csv, diag0 = TRUE) {

  conn <- as.matrix(read.csv(conn.csv, header = FALSE, sep = ","))

  if (diag0)
    diag(conn) <- 0

  conn_diag <- diag(conn)
  conn[which(lower.tri(conn))] <- 0
  conn <- conn + t(conn)
  diag(conn) <- conn_diag

  plot_adjacency_matrix(conn)

}