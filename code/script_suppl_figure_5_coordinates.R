######################################################################################################
####### Script for OTC part of Supplementary Figure 5 ################################################
####### Evaluation of estimated coordinates            ################################################
######################################################################################################

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
library(OTC)
library(ggplot2)


################# reevaluate yeast data ##############################################################
print("Evaluation of coordinates from yeast data")
data_path <- "../data/real_data/Figure6_Yeast"
output_path <- "../results"
data_sets <- c("Tom40_Cbp3", "Tom40_Mrpl4", "Tom40_Tom20", "Tom40_Tom40")

# calculate transport plans from coordinates
for (i in data_sets){
  print(paste("Processing data set", i))
  data_path_i <- file.path(data_path, i, "coordinates")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans_coordinates(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, "_coord.RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="yeast_coord")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, lower_t = 0, upper_t = 500, output_path =output_path, output_name = "yeast_suppl_fig5")


################## reevaluate HDFa data ##############################################################
print("Evaluation of coordinates from HDFa data")
# get data path
data_path <- "../data/real_data/SupplFig4_HDFa"
data_sets <- c("high_background", "low_background")
output_path <- "../results"

# calculate transport plans from coordinates
for (i in data_sets){
  data_path_i <- file.path(data_path, i, "coordinates")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_A_", files)]
  picsB <- files[grepl("_B_", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans_coordinates(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, "_coord.RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="HDFa_coord")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, lower_t = 0, upper_t = 500, output_path = output_path, output_name = "HDFa_suppl_figure5")


################## reevaluate 2D/3D data ##############################################################
print("Evaluation of coordinates from 2D/3D PSF data")
# get data path
data_path <- "../data/real_data/Figure7_2D_3D"
data_sets <- c("2D", "3D")
output_path <- "../results"

# calculate transport plans from coordinates
for (i in data_sets){
  data_path_i <- file.path(data_path, i, "coordinates")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans_coordinates(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path=output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, "_coord.RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="2D_3D")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, lower_t = 0, upper_t = 500, output_path = output_path, output_name = "2D_3D_suppl_figure5")
