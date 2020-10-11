######################################################################################################
####### Script for OTC part of Supplemenatry Figure 1 ################################################
####### Evalution of simulated random spots ##########################################################
######################################################################################################

library(OTC)

# get data path
data_path <- "../data/Simulated Data/Suppl. Fig. 1 random"

for (i in seq(10, 90, 10)){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_red", files)]
  picsB <- files[grepl("_green", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = "../results", output_name = paste("random", i, sep="_"))
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- seq(10, 90, 10)
data_list <- paste("Tplans_random_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "random_suppl_fig1")


