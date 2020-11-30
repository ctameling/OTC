rm(list = ls())
######################################################################################################
####### Script for simulated ring structure data #####################################################
#######  ####################################################
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
library(tiff)
library(ggplot2)
source("../code/corMethods.R") 

# get data path
# setup data and output path as well as data sets
data_path <- "../data/simulated_data/Ring"
output_path <- "../results"
data_set <- seq(40, 240, 40)

# evaluate tplans of all preset levels of colocalization
for (i in data_sets){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_red", files)]
  picsB <- files[grepl("_green", files)]
  #squassh_data <- files[grepl("Squassh_", files)]
  #debias_data <- files[grepl("_GILI", files)]
  #object_coloc <- files[grepl("Object_coloc_", files)]
  #data_obj_coloc <- read.csv(paste0(data_path_i,"/",object_coloc), header = TRUE, sep = ";")
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = output_path, output_name = paste("ring_structure", i, sep="_"))
  
  #------------------------ Pixel based colocalization data ----------------------------------------#
  # n <- length(picsA)
  # 
  # # Initialize coefficient arrays 
  # ICCS_M1 <- rep(NA,n)
  # ICCS_M2 <- rep(NA,n)
  # LI_DeBias <- rep(NA,n)
  # Mander_s_M1 <- rep(NA,n)
  # Mander_s_M2 <- rep(NA,n)
  # Pearson_s_Corr <- rep(NA,n)
  # Pearson_s_with_tresh <- rep(NA,n)
  # Thresholded_Overlap_Coeff_ch1 <- rep(NA,n)
  # Thresholded_Overlap_Coeff_ch2 <- rep(NA,n)
  # 
  # # Include precomputed DeBias LI Coefficient
  # debias_coloc <- read.table(paste0(data_path_i,"/",debias_data), header = TRUE, sep = ",")
  # colnames(debias_coloc) <- c("GI","LI")
  # LI_DeBias <- debias_coloc$LI
  # 
  # # Include precomputed Squassh Colocalization Coefficient
  # squassh_coloc <- read.table(paste0(data_path_i,"/",squassh_data), header = TRUE, sep =",")
  # Thresholded_Overlap_Coeff_ch1 <- squassh_coloc$ColocObjectsNumber[squassh_coloc$Channel != 1]
  # Thresholded_Overlap_Coeff_ch2 <- squassh_coloc$ColocObjectsNumber[squassh_coloc$Channel == 1]
  # 
  # for (j in 1:n) {
  #   # Read simulated Tiff Pictures
  #   picA <- tiff::readTIFF(paste0(data_path_i,"/",picsA[j]))
  #   picB <- tiff::readTIFF(paste0(data_path_i,"/",picsB[j]))
  #   
  #   # Compute Colocalization Coefficients
  #   Pearson_s_Corr[j] <- Pcor(picA, picB)
  #   Pearson_s_with_tresh[j] <- PcorP(picA, picB)$PcorP
  #   Mander_s_M1[j] <- M1(picA, picB)
  #   Mander_s_M2[j] <- M2(picA, picB)
  #   ICCS_M1[j] <- M1_ICCS(picA, picB)
  #   ICCS_M2[j] <- M2_ICCS(picA, picB)
  # }
  # 
  # coloc_data <- data.frame("Picture" = i,
  #                          "ICCS M1" = ICCS_M1,
  #                          "ICCS M2" = ICCS_M2, 
  #                          "LI (DeBias)" = LI_DeBias,
  #                          "Mander's M1" = Mander_s_M1, 
  #                          "Mander's M2" = Mander_s_M2,
  #                          "Pearson's Corr" = Pearson_s_Corr, 
  #                          "Pearson's with tresh." = Pearson_s_with_tresh, 
  #                          "Thresholded Overlap Coeff ch1" = Thresholded_Overlap_Coeff_ch1, 
  #                          "Thresholded Overlap Coeff ch2" = Thresholded_Overlap_Coeff_ch2)
  # assign(paste0("coloc_data_pixel_",i), coloc_data)
  # 
  # coefficients <- c("ICCS M1", "ICCS M2", "LI (DeBias)",
  #                   "Mander's M1", "Mander's M2",
  #                   "Pearson's Corr.", "Pearson's with tresh.",
  #                   "Thresholded Overlap Coeff ch1", "Thresholded Overlap Coeff ch2")
  # 
  # meanvalue <- rep(NA,length(coefficients))
  # stand.error <- function(x) sd(x)/sqrt(length(x))
  # standarderror <- rep(NA,length(coefficients))
  # 
  # # Compute mean Computed Colocalization level
  # meanvalue[1] <- mean(ICCS_M1) 
  # standarderror[1] <- stand.error(ICCS_M1) 
  # meanvalue[2] <- mean(ICCS_M2)  
  # standarderror[2] <- stand.error(ICCS_M2) 
  # meanvalue[3] <- mean(LI_DeBias[complete.cases(LI_DeBias)]) 
  # standarderror[3] <- stand.error(LI_DeBias[complete.cases(LI_DeBias)]) 
  # meanvalue[4] <- mean(Mander_s_M1) 
  # standarderror[4] <- stand.error(Mander_s_M1) 
  # meanvalue[5] <- mean(Mander_s_M2)     
  # standarderror[5] <- stand.error(Mander_s_M2) 
  # meanvalue[6] <- mean(Pearson_s_Corr) 
  # standarderror[6] <- stand.error(Pearson_s_Corr) 
  # meanvalue[7] <- mean(Pearson_s_with_tresh) 
  # standarderror[7] <- stand.error(Pearson_s_with_tresh) 
  # meanvalue[8] <- mean(Thresholded_Overlap_Coeff_ch1) 
  # standarderror[8] <- stand.error(Thresholded_Overlap_Coeff_ch1) 
  # meanvalue[9] <- mean(Thresholded_Overlap_Coeff_ch2) 
  # standarderror[9] <- stand.error(Thresholded_Overlap_Coeff_ch2) 
  # 
  # truecoloc <- i/100
  # frame_coefficients <- data.frame(coefficients, meanvalue, standarderror, 
  #                                  truecoloc, stringsAsFactors=FALSE)
  # assign(paste0("frame_pixel_",i), frame_coefficients)
  # 
  # #------------------------ Object based colocalization data ----------------------------------------#
  # data_obj_coloc <- read.csv(paste0(data_path_i,"/",object_coloc), header = TRUE, sep = ";")
  # 
  # Mask_center_1_inside_Mask_2 <- data_obj_coloc$X..of.coloc.....1.inside.2.
  # Ripley_s_K_of_2_coloc_with_1 <- data_obj_coloc$X..of.1.coloc..with.2..fit.
  # SODA_of_2_coloc._with_1 <- data_obj_coloc$X..of.1.coloc..with.2..SODA.
  # 
  # coloc_data <- data.frame("Picture" = i,
  #                          "Mask center 1 inside Mask 2" = Mask_center_1_inside_Mask_2,
  #                          "Ripley's K of 2 coloc. with 1" = Ripley_s_K_of_2_coloc_with_1, 
  #                          "SODA of 2 coloc. with 1" = SODA_of_2_coloc._with_1)
  # assign(paste0("coloc_data_object_",i), coloc_data)
  # 
  # coefficients <- c("Mask center 1 inside Mask 2", 
  #                   "Ripley's K function of 2 coloc. with 1", 
  #                   "SODA of 2 coloc. with 1")
  # 
  # meanvalue <- rep(NA,length(coefficients))
  # stand.error <- function(x) sd(x)/sqrt(length(x))
  # standarderror <- rep(NA,length(coefficients))
  # 
  # # Compute mean Computed Colocalization level
  # meanvalue[1] <- mean(Mask_center_1_inside_Mask_2) 
  # standarderror[1] <- stand.error(Mask_center_1_inside_Mask_2) 
  # meanvalue[2] <- mean(Ripley_s_K_of_2_coloc_with_1)  
  # standarderror[2] <- stand.error(Ripley_s_K_of_2_coloc_with_1) 
  # meanvalue[3] <- mean(SODA_of_2_coloc._with_1)   
  # standarderror[3] <- stand.error(SODA_of_2_coloc._with_1) 
  # 
  # truecoloc <- i/100
  # frame_coefficients <- data.frame(coefficients, meanvalue, standarderror, 
  #                                  truecoloc, stringsAsFactors=FALSE)
  # assign(paste0("frame_object_",i), frame_coefficients)
}

# evaluate OTC
data_list <- paste("Tplans_ring_structure_", data_sets, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="sparse_structure")

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = output_path, output_name = "sparse_structures_figure4")


# #------------------------ Pixel based colocalization data ----------------------------------------#
# # Combine data of all Colocalization levels
# coloc_pixel_complete <- rbind(coloc_data_pixel_10, coloc_data_pixel_20, coloc_data_pixel_30, coloc_data_pixel_40, coloc_data_pixel_50,
#                         coloc_data_pixel_60, coloc_data_pixel_70, coloc_data_pixel_80, coloc_data_pixel_90)
# mean_pixel_complete <- rbind(frame_pixel_10, frame_pixel_20, frame_pixel_30, frame_pixel_40, frame_pixel_50,
#                        frame_pixel_60, frame_pixel_70, frame_pixel_80, frame_pixel_90)
# 
# # Initialize Lineplot
# lineplot <- ggplot(mean_pixel_complete, aes(x = truecoloc, y = meanvalue, col = coefficients)) +  
#   geom_line() + 
#   labs(x = "True Colocalisation", y = "Computed Colocalisation") +
#   coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
#   guides(color = guide_legend(title = "")) +
#   theme(text = element_text(size = 20), legend.key.height = unit(2, 'lines'), legend.key.width = unit(2, 'lines'))+
#   geom_abline(intercept = 0, linetype = "dotted")
# 
# 
# pdf("../results/sparse_structures_figure4_pixelbased_comparison.pdf", width = 10, height = 7)
# lineplot
# dev.off()
# 
# # Save Colocalization data 
# write.csv(coloc_pixel_complete, "../results/sparse_structures_figure4_pixelbased_comparison_coloc_data.csv") 
# write.csv(mean_pixel_complete, "../results/sparse_structures_figure4_pixelbased_comparison_mean_data.csv")
# 
# #------------------------ Object based colocalization data ----------------------------------------#
# # Combine data of all Colocalization levels
# coloc_object_complete <- rbind(coloc_data_object_10, coloc_data_object_20, coloc_data_object_30, coloc_data_object_40, coloc_data_object_50,
#                         coloc_data_object_60, coloc_data_object_70, coloc_data_object_80, coloc_data_object_90)
# mean_object_complete <- rbind(frame_object_10, frame_object_20, frame_object_30, frame_object_40, frame_object_50,
#                        frame_object_60, frame_object_70, frame_object_80, frame_object_90)
# 
# # Initialize Lineplot
# lineplot <- ggplot(mean_object_complete, aes(x = truecoloc, y = meanvalue, col = coefficients)) +  
#   geom_line() + 
#   labs(x = "True Colocalisation", y = "Computed Colocalisation") +
#   coord_fixed(ratio = 1, xlim = c(0,1), ylim = c(0,1)) +
#   guides(color = guide_legend(title = "")) +
#   theme(text = element_text(size = 20), legend.key.height = unit(2, 'lines'), legend.key.width = unit(2, 'lines'))+
#   geom_abline(intercept = 0, linetype = "dotted")
# 
# 
# pdf("../results/sparse_structures_figure4_objectbased_comparison.pdf", width = 10, height = 7)
# lineplot
# dev.off()
# 
# # Save Colocalization data 
# write.csv(coloc_object_complete, "../results/sparse_structures_figure4_objectbased_comparison_coloc_data.csv") 
# write.csv(mean_object_complete, "../results/sparse_structures_figure4_objectbased_comparison_mean_data.csv")
