#%######################################################%##
#                                                          #
####             Example for Volumetry and              ####
####          Lesion Segmentation using Flows           ####
#                                                          #
##%######################################################%##

devtools::load_all("../ni.datasets/")
devtools::load_all("../utils4ni/")
devtools::load_all("../wf4ni/")
devtools::load_all(".")

#### Flow ####
flow <- complete_volumetry_flow()
flow$plot()

#### Data ####
require(RNifti)
t1_file <- tempfile(fileext = ".nii.gz")
writeNifti(image = import_dicom_folder("~/Downloads/alberto.lozano/Experiment_ID_257/scans/6_COR FSPGR 3D BUENO/")[[1]],
           file = t1_file)


#### Execute Flow ####
system.time(
  result <- flow$execute(inputs = list(T1 = t1_file),
                         desired_outputs = c("T1_bias_field",
                                             "betted_image",
                                             "segmentation",
                                             "basic_volumetry",
                                             "subcortical_parcellation",
                                             "subcortical_volumetry",
                                             "enhanced_segmentation",
                                             "parcellation",
                                             "final_volumetry"))
)

#### Inspect Results ####
num_classes <- 2
col.y <- scales::alpha(colour = scales::hue_pal()(num_classes), alpha = 0.45)

base_image <- result$T1_bias_field

ortho_plot(base_image, text = "ORIGINAL IN MNI")

ortho_plot(x = base_image,
           y = result$betted_image > 0,
           col.y = col.y,
           text = "BET")

num_classes <- 4
col.y <- scales::alpha(colour = scales::viridis_pal()(num_classes), alpha = 0.25)

ortho_plot(x = base_image,
           y = result$enhanced_segmentation,
           col.y = col.y,
           text = "FINAL SEGMENTATION")

num_classes <- 252
col.y <- scales::alpha(colour = scales::hue_pal()(num_classes), alpha = 0.45)

ortho_plot(x = base_image,
           y = result$parcellation,
           col.y = col.y,
           text = "PARCELLATION")
