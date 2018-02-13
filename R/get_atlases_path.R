get_atlases_path <- function() {

  if (nchar(system.file("nifti", package = "npbd.atlases")) > 0)
    return(system.file("nifti", package = "npbd.atlases"))

  if (nchar(system.file("nifti", package = "npbd.atlases.multiatlas")) > 0)
    return(system.file("nifti", package = "npbd.atlases"))

  return(normalizePath("~/Desktop/Repos_102/TestPipelines/npbd.atlases.multiatlas/inst/nifti"))

}