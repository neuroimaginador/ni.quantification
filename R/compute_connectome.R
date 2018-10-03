compute_connectome <- function(t1, parcellation, dwi_file) {

  # file is the nifti file with DWI image
  # We need a .bval and a .bvec files with the same name and in the same folder
  # to be able to import the DWI.

  dwi_prefix <- basename(dwi_file)
  dwi_prefix <- gsub(dwi_prefix, pattern = ".nii|.gz", replacement = "")

  dwi_folder <- dirname(dwi_file)

  bval_file <- file.path(dwi_folder, paste0(dwi_prefix), ".bval")
  bvec_file <- file.path(dwi_folder, paste0(dwi_prefix), ".bvec")

  stopifnot(file.exists(bval_file) & file.exists(bvec_file))

  require(dti)
  require(data.table)

  b <- read.table(bval_file) %>% as.vector()
  g <- read.table(bvec_file) %>% as.matrix()

  if (ncol(g) == 3) g <- t(g)

  dtiObj <- readDWIdata(gradient = g, dirlist = dwi_file, format = "NIFTI", bvalue = b)

  tensors <- dtiTensor(dtiObj, method = "linear")

  require(neurobase)
  orig <- readNIfTI(dwi_file, reorient = FALSE)
  s0 <- as.nifti(from = tensors@th0, value = orig)

  require(RNiftyReg)
  nrl <- niftyreg.linear(source = t1, target = s0, scope = "affine")

  mask <- applyAffine(source = as.nifti(parcellation, value = t1),
                      target = s0,
                      affine = nrl$affine[[1]], finalInterpolation = 0)$image

  fibers <- tracking(tensors, subsample = 1, mask = mask > 0 & mask < 300, minfa = 0.3, maxangle = 60)

  # WM_ROIs <- fiber.regions(fibers, labels = mask * (mask < 117), s0 = s0)
  # visited_regions <- fiber_visited_regions(fibers, labels = mask, s0 = s0)
  # connectome <- build_extended_connectome(visited_regions)
  # write_extended_connectome_csv(connectome, filename = file.path(folder.out, "connectome.csv"))

  res <- label.fibers(fibers, labels = mask,
                      voxel.dims = pixdim(mask)[2:4],
                      add.ends = TRUE)

  labels <- res$rois

  # Compute adjacency matrix in tractography
  conn <- matrix(0, nrow = 116, ncol = 116)
  i1 <- (labels[, 1] - 1) * 116 + labels[, 2]
  i2 <- (labels[, 2] - 1) * 116 + labels[, 1]

  s <- tabulate(i1)
  conn[seq(max(i1))] <- s
  s <- tabulate(i2)
  conn[seq(max(i2))] <- s

  conn <- pmax(conn, t(conn))

  # Computing connectomic measures:
  # For classical connectome
  connectome_graph <- build_graph_from_matrix(conn)
  connectome_parameters <- connectome_measures(connectome_graph)

  return(connectome_parameters$parameters)

}

# centroids <- cbind(aal116$x.mni, aal116$y.mni, aal116$z.mni)
# P <- as.antsImage(parcellation)
# brain <- renderImageLabels(labelsimg = thresholdImage(P, 0, 0), smoothsval = 1, alphasurf = 0.2, col = "lightgray")
# get_atlas("MNI")
# brain <-
# plotBasicNetwork(centroids, brain, weights = conn)


import_dwi_data <- function(dwi_file) {

  # file is the nifti file with DWI image
  # We need a .bval and a .bvec files with the same name and in the same folder
  # to be able to import the DWI.

  dwi_prefix <- basename(dwi_file)
  dwi_prefix <- gsub(dwi_prefix, pattern = ".nii|.gz", replacement = "")

  dwi_folder <- dirname(dwi_file)

  bval_file <- file.path(dwi_folder, paste0(dwi_prefix), ".bval")
  bvec_file <- file.path(dwi_folder, paste0(dwi_prefix), ".bvec")

  stopifnot(file.exists(bval_file) & file.exists(bvec_file))

  require(dti)
  require(data.table)

  b <- read.table(bval_file) %>% as.vector()
  g <- read.table(bvec_file) %>% as.matrix()

  if (ncol(g) == 3) g <- t(g)

  dtiObj <- readDWIdata(gradient = g, dirlist = dwi_file, format = "NIFTI", bvalue = b)

  return(dtiObj)

}

compute_tensor <- function(dtiObj) {

  require(dti)
  tensors <- dtiTensor(dtiObj, method = "linear")

  return(tensors)

}

compute_dwi_mask <- function(t1_file, parcellation, dwi_file, tensors) {

  require(neurobase)
  orig <- readNIfTI(dwi_file, reorient = FALSE)
  t1 <- readNIfTI(t1_file, reorient = FALSE)
  s0 <- as.nifti(from = tensors@th0, value = orig)

  require(RNiftyReg)
  nrl <- niftyreg.linear(source = t1, target = s0, scope = "affine")

  mask <- applyAffine(source = as.nifti(parcellation, value = t1),
                      target = s0,
                      affine = nrl$affine[[1]], finalInterpolation = 0)$image

  return(mask)

}

compute_fibers <- function(tensors, mask) {

  require(dti)

  fibers <- tracking(tensors, subsample = 1, mask = mask > 0 & mask < 300, minfa = 0.3, maxangle = 60)

  return(fibers)

}

get_anatomic_adjacency_matrix <- function(fibers, mask) {

  require(neurobase)
  res <- label.fibers(fibers, labels = mask,
                      voxel.dims = pixdim(mask)[2:4],
                      add.ends = TRUE)

  labels <- res$rois

  # Compute adjacency matrix in tractography
  conn <- matrix(0, nrow = 116, ncol = 116)
  i1 <- (labels[, 1] - 1) * 116 + labels[, 2]
  i2 <- (labels[, 2] - 1) * 116 + labels[, 1]

  s <- tabulate(i1)
  conn[seq(max(i1))] <- s
  s <- tabulate(i2)
  conn[seq(max(i2))] <- s

  conn <- pmax(conn, t(conn))

  return(conn)

}

