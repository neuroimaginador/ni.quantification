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

set_atlases_dir("/Volumes/Domingo/ni.atlases/")
set_dataset_dir("/Volumes/Domingo/dlni_data/")

#### Flow ####
flow <- lesion_segmentation_flow()
flow$plot()

#### Data ####
require(RNifti)
t1_file <- "/Volumes/Domingo/UPPMS/images/berta.sebastian/Experiment_ID_23/scans/4_Sag T1 GLOBAL/T1.nii.gz"
flair_file <- "/Volumes/Domingo/UPPMS/images/berta.sebastian/Experiment_ID_23/scans/5_Ax T2 FLAIR/FLAIR.nii.gz"


#### Execute Flow ####
system.time(
  result <- flow$execute(inputs = list(T1 = t1_file, FLAIR = flair_file),
                         desired_outputs = c("T1_bias_field",
                                             "betted_image",
                                             "segmentation",
                                             "basic_volumetry",
                                             "subcortical_parcellation",
                                             "subcortical_volumetry",
                                             "enhanced_segmentation",
                                             "parcellation",
                                             "final_volumetry",
                                             "FLAIR_1mm",
                                             "FLAIR_mni",
                                             "lesion_map"))
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

num_classes <- 2
col.y <- scales::alpha(colour = scales::hue_pal()(num_classes), alpha = 0.85)

ortho_plot(x = base_image,
           y = result$lesion_map > 1,
           col.y = col.y,
           text = "Lesions")

ortho_plot(result$FLAIR_1mm, text = "FLAIR_1mm")
ortho_plot(result$FLAIR_mni, text = "FLAIR_mni")
ortho_plot(result$T1_bias_field)
