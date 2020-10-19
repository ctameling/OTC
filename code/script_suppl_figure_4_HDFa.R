rm(list = ls())
RNGversion("3.5.3")  # change to old sampling default
######################################################################################################
####### Script for Supplementary Figure 4#############################################################
####### Evaluation of STED images with low and high background########################################
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

# get data path
data_path <- "../data/real_data/SupplFig4_HDFa"
data_sets <- c("high background", "low background")
output_path <- "../results"

############ hand picked data #########################################################################

for (i in data_sets){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_A_", files)]
  picsB <- files[grepl("_B_", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="HDFa")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = output_path, output_name = "HDFa_suppl_figure4")


############ randomly picked data #########################################################################
samples_number <- 10
seed_LB <- 25
seed_HB <- 19
relmass_factor <- 1

for (i in data_sets){
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
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB,
                                  random_sections=TRUE, n_random_sections = samples_number, relmass_factor = relmass_factor,
                                  output_path = output_path, output_name = i)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="HDFa_random")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = output_path, output_name = "HDFa_random_suppl_figure4")

