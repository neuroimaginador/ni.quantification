structural_connectome_flow <- function() {

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

  return(flow)

}