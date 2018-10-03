plot_adjacency_matrix <- function(mat, diag_off = TRUE, reorder = TRUE) {

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

  #safe_require("jonclayden/shades")

  colors <- c("white", "red") #unclass(gradient(c("white", "red"), 2, space = "Lab"))

  old_par <- par(mar = c(0, 0, 0, 0))
  image(1:M, 1:M, mat[order, rev(order)] > 0, col = colors, asp = 1, xaxs = "i")
  abline(h = M - length(LH) + 0.5, v = length(LH) + 0.5, col = "gray")
  abline(h = M - (length(LH) + length(RH)) + 0.5, v = length(LH) + length(RH) + 0.5, col = "gray")
  par(old_par)

  invisible()

}
