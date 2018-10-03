connectomic_fingerprint_flow <- function() {

  flow <- complete_volumetry_flow()

  flow$add(inputs = "DWI")

  flow$add(what = import_dwi_data,
           inputs = "DWI",
           output = "DWI_data")

  flow$add(what = compute_tensor,
           inputs = "DWI_data",
           output = "tensors")

  flow$add(what = compute_dwi_mask,
           inputs = c("T1", "parcellation", "DWI_data", "tensors"),
           output = "DWI_mask")

  flow$add(what = compute_fibers,
           inputs = c("tensors", "DWI_mask"),
           output = "fibers")

  flow$add(what = get_anatomic_adjacency_matrix,
           inputs = c("fibers", "DWI_mask"),
           output = "structural_connectome")

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

  flow$add(what = concatenate_connectomes,
           inputs = c("structural_connectome", "functional_connectome"),
           output = "multigraph_connectome")

  flow$add(what = conn_fingerprint,
           inputs = "multigraph_connectome",
           output = "connectome_fingerprint")

  return(flow)

}

multigraph <- function(c1, c2, c3 = NULL) {

  if (!is.null(c3)) {

    res <- abind::abind(c1, c2, c3, along = 3)

  } else {

    res <- abind::abind(c1, c2, along = 3)

  }

  return(res)

}

conn_fingerprint <- function(multi_conn) {

  res <- c()

  for (i in seq(dim(multi_conn)[3])) {

    G <- build_graph_from_matrix(multi_conn[, , i])

    params <- connectome_measures(G)$parameters %>% unlist() %>% unname()

    res <- res %>% append(params)

  }

  return(res)

}