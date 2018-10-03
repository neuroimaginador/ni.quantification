resting_state_connectivity <- function(fmri_file) {

  # 0 - convert data to nii in LAS orientation ( we suggest LAS if you are skipping this step )
  # 1 - slice time correction
  # 2 - motion correction, then regress out motion parameter
  # 3 - skull stripping
  # 4 - normalize data
  # 5 - regress out WM/CSF
  # 6 - lowpass filter
  # 7 - do parcellation and produce correlation matrix from label file
  #   * or split it up:
  #     7a - do parcellation from label file
  #     7b - produce correlation matrix [--func option is ignored if step 7b
  #                                  is run by itself unless --dvarsthreshold is specified, and
  #                                  --corrts overrides default location for input parcellation
  #                                  results (outputpath/corrlabel_ts.txt)]
  # 8 - functional connectivity density mapping

  require(neurobase)
  require(fslr)

  # 1
  img <- fslr::fslslicetimer(file = fmri_file)

  filename <- tempfile(pattern = "slice_timing", fileext = ".nii.gz")
  neurobase::writenii(img, filename = filename)

  # 2
  mc_img <- fslr::mcflirt(file = filename)

  # 3
  mean_func <- utils4ni::build_average_image(mc_img)

  mc_nii <- as.nifti(mean_func, value = mc_img)
  mc_file <- tempfile(pattern = "mc_nii", fileext = ".nii.gz")
  neurobase::writenii(mc_nii, filename = mc_file)

  betted_func <- fslr::fslbet(infile = mc_file)
  func_masked <- utils4ni::mask4D(mc_img@.Data, betted_func@.Data)
  func_masked <- as.nifti(func_masked, value = mc_img)

  # 4
  omat <- tempfile(pattern = "func2mni_", fileext = ".mat")
  func_masked_file <- tempfile(pattern = "func_masked_", fileext = ".nii.gz")
  neurobase::writenii(func_masked, filename = func_masked_file)

  foo <- fslr::flirt(reffile = func_masked_file, infile = get_atlas("MNI"), omat = omat)

  WM_func <- fslr::flirt_apply(infile = get_atlas("atlasWM"), reffile = func_masked_file, initmat = omat, opts = "-interp nearestneighbour")

  aal116_func <- fslr::flirt_apply(infile = get_atlas("AAL116"), reffile = func_masked_file, initmat = omat, opts = "-interp nearestneighbour")

  wm_ts <- mean_ts_by_ROI(V = func_masked, mask = WM_func > 0)

  all_ts <- mean_ts_by_ROI(V = func_masked, mask = aal116_func)

  # 5
  regressed_ts <- all_ts
  for (roi in seq(116)) {

    m <- lm(all_ts[, roi] ~ wm_ts + 1)
    regressed_ts[, roi] <- m$residuals %>% unname()

  }

  # 6 and 7
  conn <- cor(regressed_ts)
  conn <- 0.5 * log((1 + conn) / (1 - conn))
  conn[!is.finite(conn)] <- 0

  return(conn)

}

# devtools::load_all("../ni.datasets/")
# devtools::load_all("../utils4ni/")
# devtools::load_all("../wf4ni/")
# devtools::load_all()


fmri_file <- "~/Desktop/0027306/session_1/rest_1/rest.nii.gz"

# 1
slice_timing_correction <- function(fmri) {

}

# 2
motion_correction <- function(fmri) {


}

# 3
functional_bet <- function(fmri) {


}

# 4...
obtain_time_series <- function(parcellation, fmri) {

  # include regress-out of WM

}

# 7b
get_functional_connectome <- function(timeseries) {



}

functional_connectome_flow <- function() {

  flow <- complete_volumetry_flow()

  flow$add(inputs = "rsfMRI")

  flow$add(what = slice_timing_correction,
           inputs = "rsfMRI",
           output = "fMRI_stc")

  flow$add(what = motion_correction,
           inputs = "fMRI_stc",
           output = "motion_corrected")

  flow$add(what = functional_bet,
           inputs = "motion_corrected",
           output = "fMRI_betted")

  flow$add(what = obtain_time_series,
           inputs = c("parcellation", "fMRI_betted"),
           output = "timeseries")

  flow$add(what = get_functional_connectome,
           inputs = "timeseries",
           output = "functional_connectome")


  return(flow)

}