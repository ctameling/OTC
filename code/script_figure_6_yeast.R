######################################################################################################
####### Script for OTC part of Figure 6 ##############################################################
####### Evalution of yeast data ######################################################################
######################################################################################################

library(OTC)

# get data path
data_path <- "../data/Real Data/Figure 6 Yeast"

############ hand picked data #########################################################################

for (i in c("Tom40_Cbp3", "Tom40_Mrpl4", "Tom40_Tom20", "Tom40_Tom40")){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("Tom40_Cbp3", "Tom40_Mrpl4", "Tom40_Tom20", "Tom40_Tom40")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="yeast")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "yeast_figure6")


############ randomly picked data #########################################################################
seeds <- c(21, 15, 5, 19)
samples <- c(15, 25, 25, 50)
for (i in c("Tom40_Cbp3", "Tom40_Mrpl4", "Tom40_Tom20", "Tom40_Tom40")){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  set.seed(seeds[i])
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, random_sections=TRUE, n_random_sections = samples[i], output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("Tom40_Cbp3", "Tom40_Mrpl4", "Tom40_Tom20", "Tom40_Tom40")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="yeast_random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "yeast_random_figure6")


