rm(list = ls())
RNGversion("3.5.3")  # change to old sampling default
######################################################################################################
####### Script for Supplementary Figure 9 ##########################################################################
####### Evaluation of STED images with 2D and 3D PSF##################################################
######################################################################################################

library(tiff)
library(ggplot2)
source("../code/corMethods.R") 

# get data path
data_path <- "../data/real_data/Figure7_2D_3D/"

name_molecule <- "ATPB/Mic60"

############ hand picked data #########################################################################

for (i in c("2D", "3D")){
  data_path_i <- file.path(data_path, i, "128x128 sections")
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  squassh_data <- files[grepl("Squassh_", files)]
  debias_data <- files[grepl("_GILI", files)]
  
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
  
  frame_coefficients <- data.frame(coefficients, meanvalue, standarderror, 
                                   name_molecule, i, stringsAsFactors=FALSE)
  assign(paste0("frame_",i), frame_coefficients)
}

mean_complete <- rbind(frame_2D, frame_3D)
coloc_complete <- rbind(coloc_data_2D, coloc_data_3D)

mean_complete$name_molecule <- factor(mean_complete$name_molecule, levels = unique(mean_complete$name_molecule), ordered=TRUE)

# Initialize Barplot
barplot <- ggplot(mean_complete, x=i, y=meanvalue, aes(i, meanvalue, fill = coefficients)) +  
  geom_bar(stat = "identity", , position=position_dodge()) + 
  labs(title = "", x = "", y = "") +
  geom_errorbar(aes(ymin=meanvalue-standarderror, ymax=meanvalue+standarderror), 
                position = position_dodge()) +
  coord_cartesian(ylim = c(0, 2.5))

pdf("../results/2D_3D_suppl_figure9_pixelbased_comparison.pdf", width = 10, height = 7)
barplot
dev.off()

# Save Colocalization data 
write.csv(coloc_complete, "../results/2D_3D_suppl_figure9_pixelbased_comparison_coloc_data.csv") 
write.csv(mean_complete, "../results/2D_3D_suppl_figure9_pixelbased_comparison_mean_data.csv")

############ randomly picked data #########################################################################
samples_number <- 34
seed_2D <- 16
seed_3D <- 43

#seeds <- c(16, 43)
for (i in c("2D", "3D")){
  data_path_i <- file.path(data_path, i)
  files <- list.files(data_path_i)
  picsA <- files[grepl("_594", files)]
  picsB <- files[grepl("_640", files)]
  squassh_data <- files[grepl("Squassh_", files)]
  debias_data <- files[grepl("_GILI", files)]
  
  n <- length(picsA)*samples_number
  
  if(i == "2D") {
    set.seed(seed_2D) 
  } else {
    set.seed(seed_3D) 
  }
  
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
                                     relmass = min(length(which(picA > 0)), 
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

mean_complete <- rbind(frame_2D, frame_3D)
coloc_complete <- rbind(coloc_data_2D, coloc_data_3D)

mean_complete$name_molecule <- factor(mean_complete$name_molecule, levels = unique(mean_complete$name_molecule), ordered=TRUE)

# Initialize Barplot
barplot <- ggplot(mean_complete, x=i, y=meanvalue, aes(i, meanvalue, fill = coefficients)) +  
  geom_bar(stat = "identity", , position=position_dodge()) + 
  labs(title = "", x = "", y = "") +
  geom_errorbar(aes(ymin=meanvalue-standarderror, ymax=meanvalue+standarderror), 
                position = position_dodge()) +
  coord_cartesian(ylim = c(0, 2.5))

pdf("../results/2D_3D_random_suppl_figure9_pixelbased_comparison.pdf", width = 10, height = 7)
barplot
dev.off()

# Save Colocalization data 
write.csv(coloc_complete, "../results/2D_3D_random_suppl_figure9_pixelbased_comparison_coloc_data.csv") 
write.csv(mean_complete, "../results/2D_3D_random_suppl_figure9_pixelbased_comparison_mean_data.csv")
