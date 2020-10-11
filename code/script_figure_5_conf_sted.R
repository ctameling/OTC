######################################################################################################
####### Script for OTC part of Figure 5 ##############################################################
####### Evalution of confocal and STED images ########################################################
######################################################################################################


library(OTC)

# get data path
data_path <- "../data/Real Data/Figure 5 Conf STED"

############ hand picked data #########################################################################

for (i in c("Conf", "STED")){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_Mic60", files)]
  picsB <- files[grepl("_Tom20", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("Conf", "STED")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="Conf_STED")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "Conf_STED_figure5")


############ randomly picked data #########################################################################

for (i in c("Conf", "STED")){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_Mic60", files)]
  picsB <- files[grepl("_Tom20", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, random_sections=TRUE, n_random_sections = 7, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("Conf", "STED")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="Conf_STED_random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "Conf_STED_random_figure5")


