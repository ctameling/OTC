######################################################################################################
####### Script for OTC part of Figure 7 ##############################################################
####### Evalution of STED iamges with 2D and 3D PSF###################################################
######################################################################################################

library(OTC)

# get data path
data_path <- "../data/Real Data/Figure 7 2D 3D"

############ hand picked data #########################################################################

for (i in c("2D", "3D")){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("2D", "3D")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="2D_3D")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "2D_3D_figure7")


############ randomly picked data #########################################################################
seeds <- c(16, 43)
for (i in c("2D", "3D")){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  set.seed(seeds[i])
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, random_sections=TRUE, n_random_sections = 34, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("2D", "3D")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="2D_3D_random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "2D_3D_random_figure7")


