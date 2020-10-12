rm(list = ls())
RNGversion("3.5.3")  # change to old sampling default
######################################################################################################
####### Script for Supplementary Figure 4#############################################################
####### Evaluation of STED images with low and high background########################################
######################################################################################################

library(OTC)

# get data path
data_path <- "../data/real_data/SupplFig4_HDFa/"

############ hand picked data #########################################################################

for (i in c("high background", "low background")){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_A_", files)]
  picsB <- files[grepl("_B_", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("high background", "low background")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="HDFa")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "HDFa_suppl_figure4")


############ randomly picked data #########################################################################
samples_number <- 10
seed_LB <- 25
seed_HB <- 19

for (i in c("high background", "low background")){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  if(i == "high background") {
    set.seed(seed_HB) 
  } else {
    set.seed(seed_LB) 
  }
  
  # compute tplans
  set.seed(seeds[i])
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, random_sections=TRUE, n_random_sections = samples_number, output_path = "../results", output_name = i)
}

# evaluate OTC
data_path_tplans <- "../results"
data_list <- c("high background", "low background")
data_list <- paste("Tplans_", data_list, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = data_path_tplans, data_list=data_list, pxsize=pxsize, dim=dim, output_path="../results", output_name="HDFa_random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = "../results", output_name = "HDFa_random_suppl_figure4")
