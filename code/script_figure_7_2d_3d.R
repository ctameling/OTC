rm(list = ls())
RNGversion("3.5.3")  # change to old sampling default
######################################################################################################
####### Script for OTC part of Figure 7             ##################################################
####### Evalution of STED iamges with 2D and 3D PSF ##################################################
######################################################################################################


install.packages("OTC_0.1.0.tar.gz", repos = NULL, type = "source")
tryCatch(
  {
    current_path = rstudioapi::getActiveDocumentContext()$path
    setwd(dirname(current_path ))
  }, 
  error=function(cond){
    if (identical(cond, "RStudio not running")){
      this.dir <- dirname(parent.frame(2)$ofile)
      setwd(this.dir)
    }
  })
library(OTC)
library(ggplot2)
RNGversion("3.5.3")  # change to old sampling default

# get data path
data_path <- "../data/real_data/Figure7_2D_3D"
data_sets <- c("2D", "3D")
output_path <- "../results"

############ hand picked data #########################################################################

for (i in data_sets){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path=output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="2D_3D")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = output_path, output_name = "2D_3D_figure7")


############ randomly picked data #########################################################################
seed_2d <- 16
seed_3d <- 43
for (i in data_sets){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  
  # compute tplans
  if (i == "2D"){
    set.seed(seed_2d)
  }else{
    set.seed(seed_3d)
  }
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, random_sections=TRUE, n_random_sections = 34, output_path = output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="2D_3D_random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = output_path, output_name = "2D_3D_random_figure7")



