rm(list=ls())
library(GenEval)
library(cgwtools)
library(gdata)
library(combinat)
library(BGLR)
library(data.table)



map_df <- read.csv("path/map_df_Ratio0.csv")
p1 = map_df$AlleleFreq_basePB1
p2 = map_df$AlleleFreq_basePB2


read_and_combine <- function(prefix, generations, path, ratio) {
  file_list <- lapply(generations, function(i) {
    as.matrix(fread(paste0(path, prefix, "_gen", i, "_Ratio", ratio, ".txt"))[, -1])
  })
  do.call(rbind, file_list)
}

# Define common parameters
path <- "path/"
generations <- 31:39
ratio <- 0.1

# Read and combine data for each dataset
geno <- read_and_combine("geno", generations, path, ratio)
BOA1 <- read_and_combine("BOA1", generations, path, ratio)
phase1 <- read_and_combine("phase1", generations, path, ratio)
BOA2 <- read_and_combine("BOA2", generations, path, ratio)
phase2 <- read_and_combine("phase2", generations, path, ratio)


scale_geno <- function(breed_index, p, BOA, Phase) {
  # Initialize an empty list to store the scaled matrices
  scaled_phase <- vector("list", length(BOA))
  
  # Iterate over the two breed matrices (assuming BOA and Phase have same dimensions)
  for (i in seq_along(BOA)) {
    mask <- (BOA[[i]] == breed_index)  # Create a logical matrix mask
    scaled_phase[[i]] <- matrix(0, nrow = nrow(Phase[[i]]), ncol = ncol(Phase[[i]]))  # Initialize zero matrix
    scaled_phase[[i]][mask] <- Phase[[i]][mask] - p[col(mask)]  # Vectorized subtraction using col indices
    
    # Print progress every 10 SNPs
    if (ncol(Phase[[i]]) >= 10) {
      message("Processing breed ", i, " - SNPs processed in batches of 10")
    }
  }
  
  # Combine all matrices (element-wise sum)
  scaled_geno <- Reduce(`+`, scaled_phase)
  
  return(scaled_geno)
}


# Lists containing the matrices
BOA <- list(BOA1, BOA2)
Phase <- list(phase1, phase2)

# Compute scaled genotype
scaled_geno_B1 <- scale_geno(1, p1, BOA, Phase)
scaled_geno_B2 <- scale_geno(2, p2, BOA, Phase)

alpha_B1 <- as.data.frame(read.csv("path/alpha_B1_Ratio01.csv", 
                                   row.names=1))

scaled_geno_B1_noQTL = scaled_geno_B1[,-alpha_B1$QTLpos]
scaled_geno_B2_noQTL = scaled_geno_B2[,-alpha_B1$QTLpos]

pheno <- read.csv("path/pheno_gen45_Ratio01.csv", row.names=1)
ped <- read.csv("path/ped_gen45_Ratio01.csv", row.names=1)

pheno_train <- pheno[ped[ped$generation >= 31 & ped$generation <= 39, "ID"],]



BayesA_train_Trait_1 <- BGLR(y= pheno_train$V1,saveAt="car1_BA_all",
                             ETA=list(list(X=scaled_geno_B1_noQTL,model="BayesA"),
                                      list(X=scaled_geno_B2_noQTL,model="BayesA")),
                             nIter = 25000, burnIn = 10000)

save(BayesA_train_Trait_1, file = 'path/trait1bayesA_R01.RData')

