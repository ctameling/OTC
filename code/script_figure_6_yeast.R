rm(list = ls())
RNGversion("3.5.3")  # change to old sampling default
######################################################################################################
####### Script for Figure 6 ##########################################################################
####### Evaluation of yeast data #####################################################################
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
library(tiff)
library(tidyr)
source("../code/corMethods.R") 

# get data path
data_path <- "../data/real_data/Figure6_Yeast"
data_sets <- c("Tom40_Cbp3", "Tom40_Mrpl4", "Tom40_Tom20", "Tom40_Tom40")
output_path <- "../results"

############ hand picked data #########################################################################

for (i in data_sets){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  squassh_data <- files[grepl("Squassh_", files)]
  debias_data <- files[grepl("_GILI", files)]
  coloc_tessler_data <- read.csv(file.path(data_path_i, "coloc_tessler.txt"), header= TRUE, sep="\t")
  coloc_tessler_data <- coloc_tessler_data[, -c(1)]
  object_coloc <- files[grepl("Object_coloc_", files)]
  
  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, output_path = output_path, output_name = i)

  #------------------------ Pixel based colocalization data ----------------------------------------#
  n <- length(picsA)

  # Initialize coefficient arrays
  ICCS_M1 <- rep(NA,n)
  ICCS_M2 <- rep(NA,n)
  LI_DeBias <- rep(NA,n)
  Mander_s_M1 <- rep(NA,n)
  Mander_s_M2 <- rep(NA,n)
  Pearson_s_Corr <- rep(NA,n)
  Pearson_s_with_tresh <- rep(NA,n)
  Thresholded_Overlap_Coeff_ch1 <- rep(NA,n)
  Thresholded_Overlap_Coeff_ch2 <- rep(NA,n)

  # Include precomputed DeBias LI Coefficient
  debias_coloc <- read.table(paste0(data_path_i,"/",debias_data), header = TRUE, sep = ",")
  colnames(debias_coloc) <- c("GI","LI")
  LI_DeBias <- debias_coloc$LI

  # Include precomputed Squassh Colocalization Coefficient
  squassh_coloc <- read.table(paste0(data_path_i,"/",squassh_data), header = TRUE, sep =",")
  Thresholded_Overlap_Coeff_ch1 <- squassh_coloc$ColocObjectsNumber[squassh_coloc$Channel != 1]
  Thresholded_Overlap_Coeff_ch2 <- squassh_coloc$ColocObjectsNumber[squassh_coloc$Channel == 1]

  for (j in 1:n) {
    # Read simulated Tiff Pictures
    picA <- readTIFF(paste0(data_path_i,"/",picsA[j]))
    picB <- readTIFF(paste0(data_path_i,"/",picsB[j]))

    # Compute Colocalization Coefficients
    Pearson_s_Corr[j] <- Pcor(picA, picB)
    Pearson_s_with_tresh[j] <- PcorP(picA, picB)$PcorP
    Mander_s_M1[j] <- M1(picA, picB)
    Mander_s_M2[j] <- M2(picA, picB)
    ICCS_M1[j] <- M1_ICCS(picA, picB)
    ICCS_M2[j] <- M2_ICCS(picA, picB)
  }

  coloc_data <- data.frame("Picture" = i,
                           "ICCS M1" = ICCS_M1,
                           "ICCS M2" = ICCS_M2,
                           "LI (DeBias)" = LI_DeBias,
                           "Mander's M1" = Mander_s_M1,
                           "Mander's M2" = Mander_s_M2,
                           "Pearson's Corr" = Pearson_s_Corr,
                           "Pearson's with tresh." = Pearson_s_with_tresh,
                           "Thresholded Overlap Coeff ch1" = Thresholded_Overlap_Coeff_ch1,
                           "Thresholded Overlap Coeff ch2" = Thresholded_Overlap_Coeff_ch2)
  assign(paste0("coloc_data_",i), coloc_data)

  coefficients <- c("ICCS M1", "ICCS M2", "LI (DeBias)",
                    "Mander's M1", "Mander's M2",
                    "Pearson's Corr.", "Pearson's with tresh.",
                    "Thresholded Overlap Coeff ch1", "Thresholded Overlap Coeff ch2")

  meanvalue <- rep(NA,length(coefficients))
  standarderror <- rep(NA,length(coefficients))

  # Compute mean Computed Colocalization level
  meanvalue[1] <- mean(ICCS_M1)
  standarderror[1] <- stand.error(ICCS_M1)
  meanvalue[2] <- mean(ICCS_M2)
  standarderror[2] <- stand.error(ICCS_M2)
  meanvalue[3] <- mean(LI_DeBias[complete.cases(LI_DeBias)])
  standarderror[3] <- stand.error(LI_DeBias[complete.cases(LI_DeBias)])
  meanvalue[4] <- mean(Mander_s_M1)
  standarderror[4] <- stand.error(Mander_s_M1)
  meanvalue[5] <- mean(Mander_s_M2)
  standarderror[5] <- stand.error(Mander_s_M2)
  meanvalue[6] <- mean(Pearson_s_Corr)
  standarderror[6] <- stand.error(Pearson_s_Corr)
  meanvalue[7] <- mean(Pearson_s_with_tresh)
  standarderror[7] <- stand.error(Pearson_s_with_tresh)
  meanvalue[8] <- mean(Thresholded_Overlap_Coeff_ch1)
  standarderror[8] <- stand.error(Thresholded_Overlap_Coeff_ch1)
  meanvalue[9] <- mean(Thresholded_Overlap_Coeff_ch2)
  standarderror[9] <- stand.error(Thresholded_Overlap_Coeff_ch2)

  name_molecule <- i

  frame_coefficients <- data.frame(coefficients, meanvalue, standarderror,
                                   name_molecule, i, stringsAsFactors=FALSE)
  assign(paste0("frame_",i), frame_coefficients)

  #------------------------ Object based colocalization data ----------------------------------------#
  data_obj_coloc <- read.csv(paste0(data_path_i,"/",object_coloc), header = TRUE, sep = ";", dec=",")
  Mask_center_1_inside_Mask_2 <- data_obj_coloc$X..of.1..center.of.mass..inside.masks.2
  Ripley_s_K_of_2_coloc_with_1 <- data_obj_coloc$X...fit.of.K.function..of.2.coloc..with.1
  SODA_of_2_coloc._with_1 <- data_obj_coloc$X...SODA..of.2.coloc..with.1

  coloc_data <- data.frame("Picture" = i,
                           "Mask center 1 inside Mask 2" = Mask_center_1_inside_Mask_2,
                           "Ripley's K of 2 coloc. with 1" = Ripley_s_K_of_2_coloc_with_1,
                           "SODA of 2 coloc. with 1" = SODA_of_2_coloc._with_1,
                           "Coloc-Tesseler Spearmans A" = coloc_tessler_data$Spearmann.A, 
                           "Coloc-Tesseler Spearmans B" = coloc_tessler_data$Spearmann.B,
                           "Coloc-Tesseler Manders A" = coloc_tessler_data$Manders.A,
                           "Coloc-Tesseler Manders B" = coloc_tessler_data$Manders.B)
  assign(paste0("coloc_data_object_",i), coloc_data)

  coefficients <- c("Mask center 1 inside Mask 2",
                    "Ripley's K function of 2 coloc. with 1",
                    "SODA of 2 coloc. with 1",
                    "Coloc-Tesseler Spearmans A", 
                    "Coloc-Tesseler Spearmans B",
                    "Coloc-Tesseler Manders A",
                    "Coloc-Tesseler Manders B")

  meanvalue <- rep(NA,length(coefficients))
  stand.error <- function(x) sd(x)/sqrt(length(x))
  standarderror <- rep(NA,length(coefficients))

  # Compute mean Computed Colocalization level
  meanvalue[1] <- mean(Mask_center_1_inside_Mask_2)
  standarderror[1] <- stand.error(Mask_center_1_inside_Mask_2)
  meanvalue[2] <- mean(Ripley_s_K_of_2_coloc_with_1)
  standarderror[2] <- stand.error(Ripley_s_K_of_2_coloc_with_1)
  meanvalue[3] <- mean(SODA_of_2_coloc._with_1)
  standarderror[3] <- stand.error(SODA_of_2_coloc._with_1)
  # add results from coloc tesseler
  meanvalue[4] <- coloc_tessler_data$Spearmann.A
  standarderror[4] <- stand.error(coloc_tessler_data$Spearmann.A)
  meanvalue[5] <- coloc_tessler_data$Spearmann.B
  standarderror[5] <- stand.error(coloc_tessler_data$Spearmann.B)
  meanvalue[6] <- coloc_tessler_data$Manders.A
  standarderror[6] <- stand.error(coloc_tessler_data$Manders.A)
  meanvalue[7] <- coloc_tessler_data$Manders.B
  standarderror[7] <- stand.error(coloc_tessler_data$Manders.B)

  name_molecule <- i
  frame_coefficients <- data.frame(coefficients, meanvalue, standarderror,
                                   name_molecule, i, stringsAsFactors=FALSE)
  assign(paste0("frame_object_",i), frame_coefficients)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, ".RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="yeast")

# write source data to csv
write.csv(otc_curves, file=file.path(output_path, "yeast_figure6_OTC_source_data.csv"))

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path =output_path, output_name = "yeast_figure6")


#------------------------ Pixel based colocalization data ----------------------------------------#

coloc_complete <- rbind(coloc_data_Tom40_Mrpl4, coloc_data_Tom40_Cbp3,
                        coloc_data_Tom40_Tom20, coloc_data_Tom40_Tom40)
coloc_complete$Picture <- factor(coloc_complete$Picture, levels = unique(coloc_complete$Picture), ordered=TRUE)

# reshape coloc_complete data
coloc_complete <- coloc_complete %>% gather(method, value, colnames(coloc_complete)[-1])

boxplot <- ggplot(coloc_complete, aes(x=Picture, y=value, fill=method)) + 
  geom_boxplot() + ylim(0,2.5) +
  labs(title = "", x="", y="Amount of colocalization")

# save plot data
boxplot_data <- unlist(ggplot_build(boxplot)$data)
write.csv(boxplot_data, "../results/yeast_figure6_pixelbased_comparison_boxplot_data.csv")

pdf("../results/yeast_figure6_pixelbased_comparison.pdf", width = 10, height = 7)
boxplot
dev.off()

# Save Colocalization data
write.csv(coloc_complete, "../results/yeast_figure6_pixelbased_comparison_source_data.csv")

#------------------------ Object based colocalization data ----------------------------------------#
# Combine data of all Colocalization levels
coloc_object_complete <- rbind(coloc_data_object_Tom40_Mrpl4, coloc_data_object_Tom40_Cbp3,
                        coloc_data_object_Tom40_Tom20, coloc_data_object_Tom40_Tom40)
coloc_object_complete$Picture <- factor(coloc_object_complete$Picture, levels = unique(coloc_object_complete$Picture), ordered=TRUE)

# reshape coloc_complete data
coloc_object_complete <- coloc_object_complete %>% gather(method, value, colnames(coloc_object_complete)[-1])

boxplot <- ggplot(coloc_object_complete, aes(x=Picture, y=value, fill=method)) + 
  geom_boxplot() + ylim(-0.5,1) +
  labs(title = "", x="", y="Amount of colocalization")
boxplot

# save plot data
boxplot_data <- unlist(ggplot_build(boxplot)$data)
write.csv(boxplot_data, "../results/yeast_figure6_objectbased_comparison_boxplot_data.csv")

pdf("../results/yeast_figure6_objectbased_comparison.pdf", width = 10, height = 7)
boxplot
dev.off()

# Save Colocalization data
write.csv(coloc_object_complete, "../results/yeast_figure6_objectbased_comparison_source_data.csv")

############ randomly picked data #########################################################################
seed_Tom40_Cbp3 <- 21
seed_Tom40_Mrpl4 <- 15
seed_Tom40_Tom20 <- 5
seed_Tom40_Tom40 <- 19

samples_number_Tom40_Cbp3 <- 15
samples_number_Tom40_Mrpl4 <- 25
samples_number_Tom40_Tom20 <- 20
samples_number_Tom40_Tom40 <- 50

for (i in data_sets){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  squassh_data <- files[grepl("Squassh_", files)]
  debias_data <- files[grepl("_GILI", files)]

  if(i == "Tom40_Cbp3") {
    set.seed(seed_Tom40_Cbp3)
    samples_number <- samples_number_Tom40_Cbp3
  }
  if(i == "Tom40_Mrpl4") {
    set.seed(seed_Tom40_Mrpl4)
    samples_number <- samples_number_Tom40_Mrpl4
  }
  if(i == "Tom40_Tom20") {
    set.seed(seed_Tom40_Tom20)
    samples_number <- samples_number_Tom40_Tom20
  }
  if(i == "Tom40_Tom40") {
    set.seed(seed_Tom40_Tom40)
    samples_number <- samples_number_Tom40_Tom40
  }

  # compute tplans
  tplans <- OTC::calculate_tplans(data_path = data_path_i, picsA = picsA, picsB = picsB, random_sections=TRUE, n_random_sections = samples_number, output_path = output_path, output_name = paste(i, "random", sep="_"))

  #------------------------ Pixel based colocalization data ----------------------------------------#
  n <- length(picsA)*samples_number

  # Initialize coefficient arrays
  ICCS_M1 <- rep(NA,n)
  ICCS_M2 <- rep(NA,n)
  LI_DeBias <- rep(NA,n)
  Mander_s_M1 <- rep(NA,n)
  Mander_s_M2 <- rep(NA,n)
  Pearson_s_Corr <- rep(NA,n)
  Pearson_s_with_tresh <- rep(NA,n)
  Thresholded_Overlap_Coeff_ch1 <- rep(NA,n)
  Thresholded_Overlap_Coeff_ch2 <- rep(NA,n)

  # Include precomputed DeBias LI Coefficient
  debias_coloc <- read.table(paste0(data_path_i,"/",debias_data), header = TRUE, sep = ",")
  colnames(debias_coloc) <- c("GI","LI")
  LI_DeBias <- debias_coloc$LI

  # Include precomputed Squassh Colocalization Coefficient
  squassh_coloc <- read.table(paste0(data_path_i,"/",squassh_data), header = TRUE, sep =",")
  Thresholded_Overlap_Coeff_ch1 <- squassh_coloc$ColocObjectsNumber[squassh_coloc$Channel != 1]
  Thresholded_Overlap_Coeff_ch2 <- squassh_coloc$ColocObjectsNumber[squassh_coloc$Channel == 1]

  for (j in 1:length(picsA)) {
    # Read simulated Tiff Pictures
    picA <- OTC::RGBtoGray(readTIFF(paste0(data_path_i,"/",picsA[j])))
    picB <- OTC::RGBtoGray(readTIFF(paste0(data_path_i,"/",picsB[j])))

    segments <- OTC::randomsegments2(picA, picB, samplesize = samples_number,
                                     segmentlength = 128,
                                     segmentwidth = 128,
                                     relmass = 5*min(length(which(picA > 0)),
                                                     length(which(picB > 0)))/(dim(picA)[1] * dim(picB)[2]))

    for(k in 1:samples_number) {
      snip1 <- as.matrix(segments[[1]][[k]])
      snip2 <- as.matrix(segments[[2]][[k]])

      # Compute Colocalization Coefficients
      Pearson_s_Corr[(j-1)*samples_number + k] <- Pcor(snip1, snip2)
      Pearson_s_with_tresh[(j-1)*samples_number + k] <- PcorP(snip1, snip2)$PcorP
      Mander_s_M1[(j-1)*samples_number + k] <- M1(snip1, snip2)
      Mander_s_M2[(j-1)*samples_number + k] <- M2(snip1, snip2)
      ICCS_M1[(j-1)*samples_number + k] <- M1_ICCS(snip1, snip2)
      ICCS_M2[(j-1)*samples_number + k] <- M2_ICCS(snip1, snip2)
    }
  }

  coloc_data <- data.frame("Picture" = i,
                           "ICCS M1" = ICCS_M1,
                           "ICCS M2" = ICCS_M2,
                           "LI (DeBias)" = LI_DeBias,
                           "Mander's M1" = Mander_s_M1,
                           "Mander's M2" = Mander_s_M2,
                           "Pearson's Corr" = Pearson_s_Corr,
                           "Pearson's with tresh." = Pearson_s_with_tresh,
                           "Thresholded Overlap Coeff ch1" = Thresholded_Overlap_Coeff_ch1,
                           "Thresholded Overlap Coeff ch2" = Thresholded_Overlap_Coeff_ch2)
  assign(paste0("coloc_data_",i), coloc_data)

  coefficients <- c("ICCS M1", "ICCS M2", "LI (DeBias)",
                    "Mander's M1", "Mander's M2",
                    "Pearson's Corr.", "Pearson's with tresh.",
                    "Thresholded Overlap Coeff ch1", "Thresholded Overlap Coeff ch2")

  meanvalue <- rep(NA,length(coefficients))
  stand.error <- function(x) sd(x)/sqrt(length(x))
  standarderror <- rep(NA,length(coefficients))

  # Compute mean Computed Colocalization level
  meanvalue[1] <- mean(ICCS_M1)
  standarderror[1] <- stand.error(ICCS_M1)
  meanvalue[2] <- mean(ICCS_M2)
  standarderror[2] <- stand.error(ICCS_M2)
  meanvalue[3] <- mean(LI_DeBias[complete.cases(LI_DeBias)])
  standarderror[3] <- stand.error(LI_DeBias[complete.cases(LI_DeBias)])
  meanvalue[4] <- mean(Mander_s_M1)
  standarderror[4] <- stand.error(Mander_s_M1)
  meanvalue[5] <- mean(Mander_s_M2)
  standarderror[5] <- stand.error(Mander_s_M2)
  meanvalue[6] <- mean(Pearson_s_Corr)
  standarderror[6] <- stand.error(Pearson_s_Corr)
  meanvalue[7] <- mean(Pearson_s_with_tresh)
  standarderror[7] <- stand.error(Pearson_s_with_tresh)
  meanvalue[8] <- mean(Thresholded_Overlap_Coeff_ch1)
  standarderror[8] <- stand.error(Thresholded_Overlap_Coeff_ch1)
  meanvalue[9] <- mean(Thresholded_Overlap_Coeff_ch2)
  standarderror[9] <- stand.error(Thresholded_Overlap_Coeff_ch2)

  frame_coefficients <- data.frame(coefficients, meanvalue, standarderror,
                                   name_molecule, i, stringsAsFactors=FALSE)
  assign(paste0("frame_",i), frame_coefficients)
}

# evaluate OTC
data_list <- paste("Tplans_", data_sets, "_random.RData", sep="")
dim <- c(128)
pxsize <- 15
otc_curves <- OTC::evaluate_tplans(data_path = output_path, data_list=data_list, pxsize=pxsize, dim=dim, output_path=output_path, output_name="yeast_random")

# write source data to csv
write.csv(otc_curves, file=file.path(output_path, "yeast_random_figure6_OTC_source_data.csv"))

# plot otc curves
OTC::plot_otc_curves(otc_curves = otc_curves, output_path = output_path, output_name = "yeast_random_figure6")


#------------------------ Pixel based colocalization data ----------------------------------------#
coloc_complete <- rbind(coloc_data_Tom40_Mrpl4, coloc_data_Tom40_Cbp3,
                        coloc_data_Tom40_Tom20, coloc_data_Tom40_Tom40)
coloc_complete$Picture <- factor(coloc_complete$Picture, levels = unique(coloc_complete$Picture), ordered=TRUE)

# reshape coloc_complete data
coloc_complete <- coloc_complete %>% gather(method, value, colnames(coloc_complete)[-1])

boxplot <- ggplot(coloc_complete, aes(x=Picture, y=value, fill=method)) + 
  geom_boxplot() + ylim(0,2.5) +
  labs(title = "", x="", y="Amount of colocalization")

# save plot data
boxplot_data <- unlist(ggplot_build(boxplot)$data)
write.csv(boxplot_data, "../results/yeast_random_figure6_pixelbased_comparison_boxplot_data.csv")

pdf("../results/yeast_random_figure6_pixelbased_comparison.pdf", width = 10, height = 7)
boxplot
dev.off()

# Save Colocalization data
write.csv(coloc_complete, "../results/yeast_random_figure6_pixelbased_comparison_source_data.csv")

