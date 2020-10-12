#' Standard error.
#' 
#' @param x The value vector.
#' @return The standard error.
#' @export
stand.error <- function(x) {
  sd(x)/sqrt(length(x))
}

#' Pearson correlation coefficient.
#'
#' @param R,G The red and green channels (matrices of same dimension).
#' @return The pearson correlation coefficient.
#' @export
Pcor <- function(R, G){
  cor(as.vector(R), as.vector(G))
}

#' Pearson correlation with threshold.
#'
#' @param R,G The red and green channels (matrices of same dimension).
#' @param step The percentile stepwidth. Default is 0.2%.
#' @return The pearson correlation coefficient.
#' @export
PcorP <- function(R, G, step = 0.02){

  R <- as.vector(R)
  G <- as.vector(G)

  resCor <- 1
  percT <- 1

  while(resCor > 0){

    if(percT <= step) break
    # Lower percentile threshold by stepwidth.
    percT <- percT - step

    # Set pixels above threshold to zero.
    Rbelow <- R * (R <= quantile(R, probs = percT))
    Gbelow <- G * (G <= quantile(G, probs = percT))

    # Compute residual correlation.
    if(sd(Rbelow) == 0 | sd(Gbelow) == 0) break
    resCor <- cor(Rbelow, Gbelow)
  }

  # Compute Pearson correlation of pixels above threshold.
  Rabove <- R * (R > quantile(R, probs = percT))
  Gabove <- G * (G > quantile(G, probs = percT))
  PcorP <- cor(Rabove, Gabove)
  return(list(percT = percT, PcorP = PcorP))
}

#' M1 coefficient.
#' @param R,G The red and green channels (matrices of same dimension).
#' @param TG Threshold. Default is NULL; in this case the threshold is
#'           automatically determined.
#' @return The M1 coefficient.
#' @export
M1 <- function(R, G, TG = NULL){

  if(is.null(TG)){
    percT <- PcorP(R, G)$percT
    TG <- quantile(G, probs = percT)
  }

  Rcoloc <- R * (G > TG)
  sum(Rcoloc) / sum(R)
}

#' M2 coefficient.
#' @param R, G The red and green channels (matrices of same dimension).
#' @param TR Threshold. Default is NULL; in this case the threshold is
#'           automatically determined.
#' @return The M2 coefficient.
#' @export
M2 <- function(R, G, TR = NULL){

  if(is.null(TR)){
    percT <- PcorP(R, G)$percT
    TR <- quantile(R, probs = percT)
  }

  Gcoloc <- G * (R > TR)
  sum(Gcoloc) / sum(G)
}

#' M1_ICCS (spatial image cross-correlation spectroscopy) coefficients
#' @param R,G The red and green channels (matrices of same dimension).
#' @return The M1_ICCS coefficient.
#' @export
M1_ICCS <- function(R, G, treshold = FALSE, step_treshold = 0.2){

  R <- as.vector(R)
  G <- as.vector(G)

  if(treshold == FALSE) {

    # <I_R> and <I_G>
    I_R <- mean(R)
    I_G <- mean(G)

    # Compute deltaI_R(x,y) for every pixel - for every entry of R/G matrices
    deltaI_R <- R - I_R
    deltaI_G <- G - I_G

    # spatial cross-correlation amplitude r_rg(0,0), r_rr(0,0)
    r_RG <- mean(deltaI_R * deltaI_G) / (I_R * I_G)
    r_RR <- mean(deltaI_R**2) / (I_R**2)

  } else {

    resCor <- 1
    percT <- 1

    while(resCor > 0){

      if(percT <= step) break
      # Lower percentile threshold by stepwidth.
      percT <- percT - step

      # Set pixels above threshold to zero.
      Rbelow <- R * (R <= quantile(R, probs = percT))
      Gbelow <- G * (G <= quantile(G, probs = percT))

      # Compute residual correlation.
      if(sd(Rbelow) == 0 | sd(Gbelow) == 0) break
      resCor <- cor(Rbelow, Gbelow)
    }

    # Compute coefficient of pixels above threshold.
    Rabove <- R * (R > quantile(R, probs = percT))
    Gabove <- G * (G > quantile(G, probs = percT))

    # <I_R> and <I_G>
    I_Rabove <- mean(Rabove)
    I_Gabove <- mean(Gabove)

    # Compute deltaI_R(x,y) for every pixel - for every entry of R/G matrices
    deltaI_Rabove <- Rabove - I_Rabove
    deltaI_Gabove <- Gabove - I_Gabove

    # spatial cross-correlation amplitude r_rg(0,0), r_gg(0,0)
    r_RG <- mean(deltaI_Rabove * deltaI_Gabove) / (I_Rabove * I_Gabove)
    r_RR <- mean(deltaI_Rabove**2) / (I_Rabove**2)
  }

  return(r_RG / r_RR)
}

#' M2_ICCS (spatial image cross-correlation spectroscopy) coefficients
#' @param R,G The red and green channels (matrices of same dimension).
#' @return The M2_ICCS coefficient.
#' @export
M2_ICCS <- function(R, G, treshold = FALSE, step_treshold = 0.2){

  R <- as.vector(R)
  G <- as.vector(G)

  if(treshold == FALSE) {

    # <I_R> and <I_G>
    I_R <- mean(R)
    I_G <- mean(G)

    # Compute deltaI_R(x,y) for every pixel - for every entry of R/G matrices
    deltaI_R <- R - I_R
    deltaI_G <- G - I_G

    # spatial cross-correlation amplitude r_rg(0,0), r_rr(0,0)
    r_RG <- mean(deltaI_R * deltaI_G) / (I_R * I_G)
    r_GG <- mean(deltaI_G**2) / (I_G**2)

  } else {

    resCor <- 1
    percT <- 1

    while(resCor > 0){

      if(percT <= step) break
      # Lower percentile threshold by stepwidth.
      percT <- percT - step

      # Set pixels above threshold to zero.
      Rbelow <- R * (R <= quantile(R, probs = percT))
      Gbelow <- G * (G <= quantile(G, probs = percT))

      # Compute residual correlation.
      if(sd(Rbelow) == 0 | sd(Gbelow) == 0) break
      resCor <- cor(Rbelow, Gbelow)
    }

    # Compute coefficient of pixels above threshold.
    Rabove <- R * (R > quantile(R, probs = percT))
    Gabove <- G * (G > quantile(G, probs = percT))

    # <I_R> and <I_G>
    I_Rabove <- mean(Rabove)
    I_Gabove <- mean(Gabove)

    # Compute deltaI_R(x,y) for every pixel - for every entry of R/G matrices
    deltaI_Rabove <- Rabove - I_Rabove
    deltaI_Gabove <- Gabove - I_Gabove

    # spatial cross-correlation amplitude r_rg(0,0), r_gg(0,0)
    r_RG <- mean(deltaI_Rabove * deltaI_Gabove) / (I_Rabove * I_Gabove)
    r_GG <- mean(deltaI_Gabove**2) / (I_Gabove**2)
  }

  return(r_RG / r_GG)
}
