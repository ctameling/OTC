#### This script demonstrates how to use the OTC package on own data #######


library(OTC)
library(ggplot2)

# specify here the path to your images 
data_path <- "../data/demo_data"

# get a list of names of the images of protein A
pics_proteinA <- c("proteinA_1.tif", "proteinA_2.tif", "proteinA_3.tif", "proteinA_4.tif")
# get a list of names of the images of protein B
pics_proteinB <- c("proteinB_1.tif", "proteinB_2.tif", "proteinB_3.tif", "proteinB_4.tif")


# calculate transport plans (make sure you have either CPLEX available 
# or the dimension of your images is not larger than 64x64)
tplans <- OTC::calculate_tplans(data_path = data_path, picsA = pics_proteinA, picsB = pics_proteinB, output_path = "../results", output_name = "demo_data")


# evaluate transport plans to get OTC curves
# evaluate OTC
dim <- c(128) # set the dimension of your images
pxsize <- 15 # the size of one pixel in nm
# the data_list will be list(Tplans_<output_name of calculate_tplans>.Rdata), make sure that you have a list format
otc_curves <- OTC::evaluate_tplans(data_path = "../results", data_list=list("Tplans_demo_data.RData"), pxsize=pxsize, dim=dim, output_path="../results", output_name="demo")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path ="../results", output_name = "demo")

