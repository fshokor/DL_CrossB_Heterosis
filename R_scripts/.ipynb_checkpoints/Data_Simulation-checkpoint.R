

rm(list=ls());gc()
memory.limit(size=25000)
Sys.setenv(LANG="eng")
setwd("/travail/fshokor/CrossBred")
library(GenEval)
library(cgwtools)
library(gdata)
library(combinat)
library(BGLR)
library(ggplot2)
library(reshape2)
library(factoextra)
library(data.table)
library(tidyr)
library(dplyr)
library(matrixStats)
library(MASS)

#############################################################################################################
# function to simulate data for a base population                                                           #
#############################################################################################################

SimBasePop <- function(n,m,nQTL,mu=0,Vy=1,h2=0.5,nbreed=1,
                       corQTLbreed=1,rho_g=0,rho_e=0,zero_prob=0.3,lambda=10,
                       namedir=NULL,namefile=NULL,phenofile="",verbose=TRUE){
  if(verbose){
    message("\n|-----------------|-------------------------------------------------|")
    message("|                 | Thank you for using SimCattle to simulate your  |")
    message("|          (____) | cattle population. This package is still under  |")
    message("|           (oo)  | development and some bugs may still be present. |")
    message("|   /-------(^^)  | We kindly ask you to report issues to the email |")
    message("|  / |        |   | address: beatriz.castro-dias-cuyabano@inrae.fr  |")
    message("| *  |  ______|   |                                                 |")
    message("|    ||W     ||   | Developers:                                     |")
    message("|    ^^      ^^   | Beatriz Cuyabano and Pascal Croiseau            |")
    message("|-----------------|-------------------------------------------------|")
    message("|     STEP 01     |           SIMULATING BASE POPULATION            |")
    message("|-----------------|-------------------------------------------------|\n")
  }
  if(is.null(namedir)){
    namedir <- paste0(m/1000,"k_",nQTL,"QTL_",n,"nPB")
  }
  if(is.null(namefile)){
    namefile <- "00BasePop"
  }
  if(!file.exists(namedir)){
    dir.create(namedir)
  }
  if(!file.exists(paste0(namedir,"/pheno",namefile,phenofile,".Rdata"))){
    #----------#
    # pedigree #
    #----------#
    if(!file.exists(paste0(namedir,"/ped",namefile,".Rdata"))){
      if(verbose){
        message("Initializing pedigree information:")
      }
      ped <- data.frame()
      for(b in 1:nbreed){
        ped <- rbind(ped,cbind(simPopInfo(n=n,n.gen=0,pre.id=paste0("PB",b)),breed=rep("pure",n)))
      };rm(b)
      ped <- cbind(ped,kronecker(diag(1,nbreed),rep(1,n)))
      names(ped)[7:ncol(ped)] <- paste0("B",1:nbreed)
      # save simulated data
      if(verbose){
        message("  Saving simulated data ... ",appendLF=FALSE)
      }
      save(ped,file=paste0(namedir,"/ped",namefile,".Rdata"))
      if(verbose){
        message("DONE")
      }
    } else{
      if(verbose){
        message("\nPedigree file '",namedir,"/ped",namefile,".Rdata' found:")
        message("  Loading data ... ",appendLF=FALSE)
      }
      load(paste0(namedir,"/ped",namefile,".Rdata"))
      message("DONE")
    }
    #-----------#
    # genotypes #
    #-----------#
    if(!file.exists(paste0(namedir,"/geno",namefile,".Rdata"))){
      if(verbose){
        message("Genotypes in LD for base population:")
      }
      map <- data.frame()
      # proportions of SNPs per chromosome, following the cattle genome
      p <- c(0.063,0.053,0.048,0.048,0.040,0.050,0.043,0.046,0.040,0.043,0.042,0.034,0.034,0.035,0.032,0.033,0.032,0.026,0.027,0.030,0.027,0.024,0.022,0.025,0.020,0.022,0.019,0.020,0.022)
      # two objects to save genome, one phased for each chromosome, and another with allele counts (0,1,2)
      geno <- list(PHASED=vector("list",29),BOA=vector("list",29),M012=numeric(0))
      names(geno$PHASED) <- paste0("chr",c(rep("0",9),rep("",20)),1:29)
      names(geno$BOA) <- paste0("chr",c(rep("0",9),rep("",20)),1:29)
      # simulate data chromosome by chromosome and generate map
      for(chr in 1:length(p)){
        if(verbose){
          message("  Chromosome ",ifelse(chr < 10," ",""),chr," ... ",appendLF=FALSE)
        }
        geno$PHASED[[chr]] <- list(numeric(0),numeric(0))
        geno$BOA[[chr]] <- list(numeric(0),numeric(0))
        aux <- data.frame(SNP=nrow(map)+1:round(m*p[chr]),chrom=rep(chr,round(m*p[chr])),pos=((1:round(m*p[chr]))*300+round(runif(round(m*p[chr]),1,5))))
        for(b in 1:nbreed){
          # simulate genotypes evolving five generations under random mating for stabilization
          if(!any(ls() == "pedtmp")){
            pedtmp <- vector("list",nbreed)
          }
          if(is.null(pedtmp[[b]])){
            pedtmp[[b]] <- simPopInfo(n,n.gen=5)
          }
          tmp <- simGeno(n,runif(round(m*p[chr]),0.05,0.5),LD=TRUE,phased=TRUE)
          tmp <- lapply(popBreed(tmp,pedtmp[[b]]),subset,subset=((1:nrow(pedtmp[[b]])) > 5*n))
          geno$PHASED[[chr]][[1]] <- rbind(geno$PHASED[[chr]][[1]],tmp[[1]])
          geno$PHASED[[chr]][[2]] <- rbind(geno$PHASED[[chr]][[2]],tmp[[2]])
          geno$BOA[[chr]][[1]] <- rbind(geno$BOA[[chr]][[1]],matrix(10*b,n,round(m*p[chr])))
          geno$BOA[[chr]][[2]] <- rbind(geno$BOA[[chr]][[2]],matrix(10*b,n,round(m*p[chr])))
          aux <- cbind(aux,apply(tmp[[1]]+tmp[[2]],2,mean)/2)
          names(aux)[ncol(aux)] <- paste0("AlleleFreq_basePB",b)
          rm(tmp);gc()
        };rm(b)
        map <- rbind(map,aux);rm(aux)
        geno$BOA[[chr]][[1]] <- geno$BOA[[chr]][[1]]+geno$PHASED[[chr]][[1]]
        geno$BOA[[chr]][[2]] <- geno$BOA[[chr]][[2]]+geno$PHASED[[chr]][[2]]
        geno$M012 <- cbind(geno$M012,geno$PHASED[[chr]][[1]]+geno$PHASED[[chr]][[2]])
        if(verbose){
          message("DONE")
        }
      };rm(chr,p,pedtmp)
      if(verbose){
        message("  Map has been created.")
      }
      geno$QTL <- sort(sample(ncol(geno$M012),nQTL))
      if(verbose){
        message("  QTL were sampled for all phenotypes.")
      }
      # better organize genotype info
      for(chr in 2:29){
        geno$PHASED[[1]][[1]] <- cbind(geno$PHASED[[1]][[1]],geno$PHASED[[2]][[1]])
        geno$PHASED[[1]][[2]] <- cbind(geno$PHASED[[1]][[2]],geno$PHASED[[2]][[2]])
        geno$PHASED <- geno$PHASED[-2];gc()
        geno$BOA[[1]][[1]] <- cbind(geno$BOA[[1]][[1]],geno$BOA[[2]][[1]])
        geno$BOA[[1]][[2]] <- cbind(geno$BOA[[1]][[2]],geno$BOA[[2]][[2]])
        geno$BOA <- geno$BOA[-2];gc()
      };rm(chr)
      geno$PHASED <- geno$PHASED[[1]]
      geno$BOA <- geno$BOA[[1]]
      colnames(geno$PHASED[[1]]) <- paste0("SNP",1:nrow(map))
      colnames(geno$PHASED[[2]]) <- paste0("SNP",1:nrow(map))
      colnames(geno$BOA[[1]]) <- paste0("SNP",1:nrow(map))
      colnames(geno$BOA[[2]]) <- paste0("SNP",1:nrow(map))
      colnames(geno$M012) <- paste0("SNP",1:nrow(map))
      rownames(geno$PHASED[[1]]) <- ped$ID
      rownames(geno$PHASED[[2]]) <- ped$ID
      rownames(geno$BOA[[1]]) <- ped$ID
      rownames(geno$BOA[[2]]) <- ped$ID
      rownames(geno$M012) <- ped$ID
      # save simulated data
      if(verbose){
        message("  Saving simulated data ... ",appendLF=FALSE)
      }
      save(map,geno,file=paste0(namedir,"/geno",namefile,".Rdata"))
      if(verbose){
        message("DONE")
      }
    } else{
      if(verbose){
        message("Genotype file '",namedir,"/geno",namefile,".Rdata' found:")
        message("  Loading data ... ",appendLF=FALSE)
      }
      load(paste0(namedir,"/geno",namefile,".Rdata"))
      message("DONE")
    }
    #------------#
    # phenotypes #
    #------------#
    if(verbose){
      message("\nPhenotypes",appendLF=FALSE)
      if(class(rho_g)[1] == "matrix"){
        if(ncol(rho_g) > 1){
          if(any(rho_g != diag(1,ncol(rho_g)))){
            message(" with correlated breeding values",appendLF=FALSE)
          }
          message(" for multiple traits",appendLF=FALSE)
        }
      }
      message(" ... ",appendLF=FALSE)
    }
    pheno <- list(a=vector("list",nbreed),g=numeric(0),e=NULL,y=NULL)
    names(pheno$a) <- paste0("breed",1:nbreed)
    pheno$mu <- mu
    pheno$Vy <- Vy
    pheno$h2 <- h2
    # pheno$varDVarGratio <- varDVarGratio
    # pheno$h2d <- H2*varDVarGratio
    # pheno$h2a <- H2-pheno$h2d
    if(class(rho_g)[1] == "matrix"){
      if(ncol(rho_g) > 1){
        pheno$rho_g <- rho_g
      }
    } else rho_g <- diag(1,length(Vy))
    if(class(rho_e)[1] == "matrix"){
      if(ncol(rho_e) > 1){
        pheno$rho_e <- rho_e
      }
    } else rho_e <- diag(1,length(Vy))
    pheno$corQTLbreed <- corQTLbreed
    rho_QTL <- diag(1,nbreed)
    rho_QTL[upper.tri(rho_QTL)] <- corQTLbreed
    rho_QTL[lower.tri(rho_QTL)] <- corQTLbreed
    rho_g_sim <- kronecker(rho_QTL,rho_g)
    rm(rho_QTL)
    s2pq <- numeric(0)
    atmp <- mvrnorm(length(geno$QTL),rep(0,ncol(rho_g_sim)),rho_g_sim)
    for(b in 1:nbreed){
      pheno$a[[b]] <- atmp[,1:ncol(rho_g)]
      atmp <- atmp[,-(1:ncol(rho_g))]
      s2pq <- sum(2*map[geno$QTL,3+b]*(1-map[geno$QTL,3+b]))
      tmp <- matrix(sqrt(Vy*pheno$h2/s2pq),nrow(pheno$a[[b]]),ncol(pheno$a[[b]]),byrow=TRUE)
      pheno$a[[b]] <- matrix(scale(pheno$a[[b]]),nrow(scale(pheno$a[[b]])),ncol(scale(pheno$a[[b]])))
      pheno$a[[b]] <- pheno$a[[b]]*tmp
      rm(s2pq,tmp)
      pheno$g <- rbind(pheno$g,(geno$M012[ped[,6+b] == 1,geno$QTL]-matrix(2*map[geno$QTL,3+b],sum(ped[,6+b] == 1),length(geno$QTL),byrow=TRUE))%*%pheno$a[[b]])
    };rm(b,atmp,rho_g_sim)
    
    # Heterosis effects for all QTLs and traits
    # Extract all breeds additive effects
    breeds_list <- pheno[["a"]]
    
    # Ensure all breeds have the same dimensions
    breed_dims <- lapply(breeds_list, dim)
    stopifnot(length(unique(breed_dims)) == 1)  # Ensure all breeds have the same shape
    breed_dfs <- lapply(breeds_list, as.data.frame)
    
    # Function to generate a sampled column
    generate_sample_column <- function(...) {
      breed_cols <- list(...)  # Collect all breed columns as a list
      n <- length(breed_cols[[1]])  # Number of rows in the column
      max_val <- max(unlist(breed_cols))  # Get max value across all breeds for this column
      
      # Sample from exponential distribution
      nonzero_values <- rexp(n, rate = lambda)
      nonzero_values <- nonzero_values / max(nonzero_values) * max_val  # Scale
      
      # Introduce zeros
      sample_vector <- ifelse(runif(n) < zero_prob, 0, nonzero_values)
      return(sample_vector)
    }
    
    # Apply function across all columns dynamically using do.call()
    sampled_list <- do.call(Map, c(f = generate_sample_column, breed_dfs))
    sampled_df <- as.data.frame(sampled_list)  # Convert back to dataframe
    
    # Convert sampled_df to a numeric matrix for matrix multiplication
    heterosis_effects <- as.matrix(sapply(sampled_df, as.numeric))
    pheno$h <- heterosis_effects
    
    if(verbose){
      message("DONE")
    }
    # simulate residuals correlated between traits, but independent for each individual
    pheno$e <- matrix(mvrnorm(n=nbreed*n,mu=rep(0,ncol(rho_e)),Sigma=rho_e),nbreed*n,ncol(rho_e))
    # re-scale residuals to satisfy Var(e) = (1-h2)*Var(y)
    pheno$e <- matrix(sqrt(Vy*(1-pheno$h2)),nbreed*n,ncol(pheno$e),byrow=TRUE)*pheno$e/matrix(apply(pheno$e,2,sd),nbreed*n,ncol(pheno$e),byrow=TRUE)
    # final phenotypes
    pheno$y <- pheno$g+pheno$e+matrix(mu,nbreed*n,ncol(pheno$g),byrow=TRUE)
    # save simulated data
    if(verbose){
      message("  Saving simulated data ... ",appendLF=FALSE)
    }
    save(pheno,file=paste0(namedir,"/pheno",namefile,phenofile,".Rdata"))
    tmp <- c(paste0(c("geno","ped"),namefile,".Rdata"),paste0("pheno",namefile,phenofile,".Rdata"))
    if(verbose){
      message("DONE\n")
      message("|-------------------------------------------------------------|")
      message("| Summary of created files                                    |")
      message("|-------------------------------------------------------------|\n")
      message("Folder in which generated files are saved:")
      message("  ",getwd(),"/",namedir,"\n")
      message("Files created:")
      message("  ",tmp[1])
      message("  ",tmp[2])
      message("  ",tmp[3],"\n")
    }
    return(list(path=paste0(getwd(),"/",namedir),files=list(geno=tmp[1],ped=tmp[2],pheno=tmp[3])))
  } else{
    tmp <- c(paste0(c("geno","ped"),namefile,".Rdata"),paste0("pheno",namefile,phenofile,".Rdata"))
    message("Folder already exists with simulated data:")
    message("  ",getwd(),"/",namedir,"\n")
    message("Files summary:")
    message("  ",tmp[1]," - ",ifelse(any(dir(namedir) == tmp[1]),"found,","MISSING!"))
    message("  ",tmp[2]," - ",ifelse(any(dir(namedir) == tmp[2]),"found,","MISSING!"))
    message("  ",tmp[3]," - ",ifelse(any(dir(namedir) == tmp[3]),"found.","MISSING!"),"\n")
    if(sum(tmp[1:2] %in% dir(namedir)) < 2){
      message("|-------------------------------------------------------------|")
      message("| !!! WARNING !!!                                             |")
      message("| NECESSARY FILE(S) MAY HAVE BEEN MOVED, RENAMED, OR DELETED  |")
      message("|-------------------------------------------------------------|")
      message("| YOU WILL NOT BE ABLE TO PROCEED WITHOUT THE MISSING FILE(S) |")
      message("|-------------------------------------------------------------|\n")
    }
    return(list(path=paste0(getwd(),"/",namedir),files=list(geno=tmp[1],ped=tmp[2],pheno=tmp[3])))
  }
}

#############################################################################################################
LoadPop <- function(BasePop,ngen,varDVarGratio,verbose=TRUE){
  if(ngen == 0){
    load(paste0(BasePop$path,"/",BasePop$files$geno),verbose=verbose)
    load(paste0(BasePop$path,"/",BasePop$files$ped),verbose=verbose)
    load(paste0(BasePop$path,"/",BasePop$files$pheno),verbose=verbose)
    if(length(geno) > 2){
      geno$PHASED <- NULL
      geno$M012 <- NULL
      gc()
    }
    return(list(map=map,geno=geno,ped=ped,pheno=pheno,BasePop=BasePop))
  } else{
    load(paste0(BasePop$path,"/geno",ifelse(ngen < 10,"0",""),ngen,"_Ratio",varDVarGratio,".Rdata"))
    load(paste0(BasePop$path,"/ped",ifelse(ngen < 10,"0",""),ngen,"_Ratio",varDVarGratio,".Rdata"))
    load(paste0(BasePop$path,"/pheno",ifelse(ngen < 10,"0",""),ngen,"_Ratio",varDVarGratio,".Rdata"))
    if(length(geno) > 2){
      geno$PHASED <- NULL
      geno$M012 <- NULL
      gc()
    }
    return(list(map=map,geno=geno,ped=ped,breedProp=breedProp,pheno=pheno,BasePop=BasePop))
  }
}


#############################################################################################################
EvolvePop <- function(PopData,ngen=1,startXBGen=1,selTraitPB=NULL,
                      selTraitXB=NULL,selIntensityPBFemales=NULL,
                      selIntensityPBMales=NULL,selIntensityXB=NULL,
                      male_gen = 4, female_gen = 2, propMalesPB4XB = 0.1, 
                      propFemalesPB4XB = 0.1, varDVarGratio=0, verbose=TRUE){
  if(verbose){
    message("\n|-----------------|-------------------------------------------------|")
    message("|                 | Thank you for using SimCattle to simulate your  |")
    message("|          (____) | cattle population. This package is still under  |")
    message("|           (oo)  | development and some bugs may still be present. |")
    message("|   /-------(^^)  | We kindly ask you to report issues to the email |")
    message("|  / |        |   | address: beatriz.castro-dias-cuyabano@inrae.fr  |")
    message("| *  |  ______|   |                                                 |")
    message("|    ||W     ||   | Developers:                                     |")
    message("|    ^^      ^^   | Beatriz Cuyabano and Pascal Croiseau            |")
    message("|-----------------|-------------------------------------------------|")
    message("|     STEP 02     |              EVOLVING GENERATIONS               |")
    message("|-----------------|-------------------------------------------------|\n")
  }
  for(j in 1:length(PopData)){
    assign(names(PopData)[j],PopData[[j]])
    PopData[[j]] <- NA;gc()
  };rm(j,PopData)
  #######################
  # loop on generations #
  #######################
  genfiles <- max(ped$generation)+(1:ngen)
  for(k in genfiles){
    if(verbose){
      message("Generation ",ifelse(k < 10," ",""),k,":")
    }
    #------------------------#
    # pedigree and genotypes #
    #------------------------#
    if(verbose){
      message("  Pedigree information ... ",appendLF=FALSE)
    }
    # Get the available generations in the data
    available_gens <- sort(unique(ped$generation))
    
    # Determine the generations to consider for males and females
    male_gen_start <- if (any(available_gens <= (k - male_gen))) {
      max(available_gens[available_gens <= (k - male_gen)])
    } else {
      min(available_gens)  # Fallback to the earliest generation
    }
    
    female_gen_start <- if (any(available_gens <= (k - female_gen))) {
      max(available_gens[available_gens <= (k - female_gen)])
    } else {
      min(available_gens)  # Fallback to the earliest generation
    }
    
    # Select individuals based on the determined generations
    pedtmp <- ped[
      (ped$generation >= female_gen_start & ped$generation < k & ped$sex == 0) |  # Females: 2 generations
        (ped$generation >= male_gen_start & ped$generation < k & ped$sex == 1),     # Males: 4 generations
    ]
    pedtmp$generation <- 0
    # perform selection of pure-breeds, on both males and females
    # Pre-initialize the pedaux structure for all b values
    if (!exists("pedaux")) {
      # Initialize the PB and XB lists with appropriate length (based on the number of columns)
      pedaux <- list(PB = vector("list", ncol(ped) - 6),
                     XB = vector("list", ncol(ped) - 6))  
    }
    
    # Loop through columns representing different traits/breeds
    for (b in 1:(ncol(ped) - 6)) {
      
      # Sort individuals from the previous generation based on their phenotypes for the selected trait
      selected_trait <- selTraitPB[b]  # Get the trait to be used for sorting
      phenotypes <- pheno$y[ped$ID %in% pedtmp$ID, selected_trait]  # Extract phenotypes for individuals in the previous generation
      sorted_indices <- sort(phenotypes, decreasing = TRUE, index.return = TRUE)$ix  # Sort individuals by phenotype in descending order
      
      # Filter individuals by sex and breed, and group them into males and females
      males <- sorted_indices[pedtmp[sorted_indices, 6 + b] == 1 & pedtmp$sex[sorted_indices] == 1]  # Get indices of selected males
      females <- sorted_indices[pedtmp[sorted_indices, 6 + b] == 1 & pedtmp$sex[sorted_indices] == 0]  # Get indices of selected females
      
      # Select the top individuals based on selection intensity
      top_count_males <- round(selIntensityPBMales[b] * length(males))  # Calculate the number of top males to select
      top_count_females <- round(selIntensityPBFemales[b] * length(females))  # Calculate the number of top females to select
      slected_males_PB <- males[1:top_count_males]  # Select the top individuals males 
      slected_females_PB <- females[1:top_count_females]  # Select the top individuals females
      
      # If the current generation is eligible for crossbreeding
      if (k >= startXBGen) {
        crossbreed_females_count <- round(propFemalesPB4XB * length(slected_females_PB))  # Calculate n% of the females for crossbreeding
        aux4X <- sample(length(slected_females_PB), crossbreed_females_count)  # Sample n% of the top females for crossbreeding
      }
      
      # Update pedaux structure based on the current generation and crossbreeding status
      if (k < startXBGen) {
        # For generations before crossbreeding, assign purebred males and females to pedaux$PB
        pedaux$PB[[b]] <- rbind(pedtmp[slected_males_PB, ], pedtmp[slected_females_PB, ])
      } else {
        # For generations after crossbreeding, assign purebred males and remaining purebred females to pedaux$PB
        pedaux$PB[[b]] <- rbind(pedtmp[slected_males_PB, ], pedtmp[slected_females_PB[-aux4X], ])
        # Assign the selected crossbred females to pedaux$XB
        pedaux$XB[[b]] <- pedtmp[slected_females_PB[aux4X], ]
      }
    };rm(b)
    
    if (k >= startXBGen) {
      for(b in 1:(ncol(ped)-6)){
        tmp <- data.frame()
        for(bx in (1:(ncol(ped)-6))[-b]){
          # Calculate the number of rows to sample (n% of the male rows)
          num_to_sample <- round(propMalesPB4XB * sum(pedaux$PB[[bx]]$sex == 1))
          # Randomly sample n% of the male rows
          sampled_rows <- sample(which(pedaux$PB[[bx]]$sex == 1), num_to_sample)
          # Use the sampled row indices to subset the data frame
          tmp <- rbind(tmp, pedaux$PB[[bx]][sampled_rows,])
        };rm(bx)
        pedaux$XB[[b]] <- rbind(pedaux$XB[[b]],tmp)
        rm(tmp)
      };rm(b)
    }
    if(k > startXBGen){
      # sort crossbred individuals by their phenotypes
      if(selTraitXB != 0){
        tmp <- sort(pheno$y[ped$ID %in% pedtmp$ID,selTraitXB],decreasing=TRUE,index.return=TRUE)$ix
        # select females
        aux <- tmp[pedtmp$breed[tmp] != "pure" & pedtmp$sex[tmp] == 0]
        aux <- aux[1:round(selIntensityXB*length(aux))]
      } else{
        for(j in 1:ncol(pheno$y)){
          tmp <- sort(pheno$y[ped$ID %in% pedtmp$ID,j],decreasing=TRUE,index.return=TRUE)$ix
          if(!any(ls() == "aux")){
            aux <- tmp[pedtmp$breed[tmp] != "pure" & pedtmp$sex[tmp] == 0]
          } else aux <- cbind(aux,tmp[pedtmp$breed[tmp] != "pure" & pedtmp$sex[tmp] == 0])
          rm(tmp)
        };rm(j)
        tmpfun <- function(x){
          return(length(unique(as.numeric(aux[1:x,])[duplicated(as.numeric(aux[1:x,]))])))
        }
        tmp <- max(which(apply(matrix(1:nrow(aux),ncol=1),1,tmpfun) <= round(selIntensityXB*nrow(aux))))
        aux <- unique(as.numeric(aux[1:tmp,])[duplicated(as.numeric(aux[1:tmp,]))])
        rm(tmpfun,tmp)
      }
      pedaux$XBF2 <- pedtmp[aux,]
      if (!exists("current_breed_index")) {
        current_breed_index <- 1
      }
      pedaux$XBF2 <- rbind(
        pedaux$XBF2,
        pedaux$PB[[current_breed_index]][pedaux$PB[[current_breed_index]]$sex == 1, ]
      )
      current_breed_index <- (current_breed_index %% length(pedaux$PB)) + 1
    }
    list2data <- function(x){
      y <- data.frame()
      for(i in 1:length(x)){
        y <- rbind(y,x[[i]])
      };rm(i)
      rownames(y) <- 1:nrow(y)
      return(y)
    }
    # define matings for new pure-breeds
    if (k < startXBGen) {
      size_off_PB <- sum(ped$generation == 0)/length(pedaux$PB)
    } else if (k == startXBGen) {
      size_off_XB <- length(pedaux$XB[[1]]$ID[pedaux$XB[[1]][, 6 + 1] == 1])
      size_off_PB <- sum(ped$generation == 0)/length(pedaux$PB) - size_off_XB
    } else if (k > startXBGen) {
      size_off_XB <- length(pedaux$XB[[1]]$ID[pedaux$XB[[1]][, 6 + 1] == 1])
      size_off_PB <- sum(ped$generation == 0)/length(pedaux$PB) * 0.7
      size_off_XBF2 <- sum(ped$generation == 0) - size_off_PB*length(pedaux$PB) - size_off_XB*length(pedaux$PB)
    }
    for(b in 1:length(pedaux$PB)){
      rownames(pedaux$PB[[b]]) <- 1:nrow(pedaux$PB[[b]])
      tmp <- simPopInfo(c(nrow(pedaux$PB[[b]]), size_off_PB),n.gen=1,pop.info=pedaux$PB[[b]][,1:5],pre.id=paste0("PB",b))
      tmp <- tmp[tmp$generation == 1,]
      tmp$ID <- paste0(substr(tmp$ID,1,4),ifelse(k < 10,"00",ifelse(k < 100,"0","")),k,substr(tmp$ID,8,14))
      tmp <- cbind(tmp,breed=rep("pure",nrow(tmp)))
      aux <- as.data.frame(matrix(0,nrow(tmp),ncol(pedaux$PB[[b]])-ncol(tmp)))
      names(aux) <- paste0("B",1:ncol(aux))
      aux[,b] <- 1
      tmp <- cbind(tmp,aux)
      pedaux$PB[[b]] <- tmp
      rm(tmp,aux)
    };rm(b)
    pedaux$PB <- list2data(pedaux$PB)
    # define matings for new F1 crosses
    if (k >= startXBGen) {  # If the current generation is greater than or equal to the start of crossbreeding
      # Loop through each crossbred population (XB) for each breed 'b'
      for (b in 1:length(pedaux$XB)) {
        # Update the row names of the current XB population to ensure they are unique
        rownames(pedaux$XB[[b]]) <- 1:nrow(pedaux$XB[[b]])
        # Simulate the next generation of crossbred individuals
        tmp <- simPopInfo(
          c(nrow(pedaux$XB[[b]]), size_off_XB),  
          n.gen = 1,                  
          pop.info = pedaux$XB[[b]][, 1:5],  # Use the current crossbred population (XB) as the base for new generation
          pre.id = "XF1"              # Prepend the ID with "XF1" for the new crossbred generation
        )
        # Remove the individuals from the base population (founders) and keep only the new generation
        tmp <- cbind(
          tmp[-(1:nrow(pedaux$XB[[b]])), 1:2],  # Keep the new individuals and their first two columns
          pedaux$XB[[b]]$ID[pedaux$XB[[b]][, 6 + b] == 1],  # Assign the ID of parent 2 based on breed-specific selection
          tmp[-(1:nrow(pedaux$XB[[b]])), 4:5]  # Keep the other columns (excluding parent1 info from base)
        )
        # Rename the third column to "parent2"
        names(tmp)[3] <- "parent2"
        # Adjust the IDs of the new generation individuals by adding the current generation number to the ID
        tmp$ID <- paste0(
          substr(tmp$ID, 1, 4),         # Keep the first part of the ID
          ifelse(k < 10, "00", ifelse(k < 100, "0", "")),  # Add leading zeros to ensure correct ID formatting
          k,                            # Add the generation number to the ID
          substr(tmp$ID, 8, 14)          # Keep the remaining part of the ID
        )
        # Add a new column labeling these individuals as 'cross' breed
        tmp <- cbind(tmp, breed = rep("cross", nrow(tmp)))
        # Create an empty data frame to store breed-specific information (B1, B2, etc.)
        aux <- as.data.frame(matrix(0, nrow(tmp), ncol(pedaux$XB[[b]]) - ncol(tmp)))
        names(aux) <- paste0("B", 1:ncol(aux))  # Name the columns as B1, B2, etc.
        # For each new individual, calculate the mean of breed information (B1, B2, etc.) from its parents
        for (i in 1:nrow(tmp)) {
          # Extract the breed information from the parents and calculate the mean for each trait (B1, B2, etc.)
          aux[i, ] <- apply(
            pedaux$XB[[b]][pedaux$XB[[b]]$ID %in% c(tmp$parent1[i], tmp$parent2[i]), 7:ncol(pedaux$XB[[b]])],
            2, mean
          )
        }
        rm(i)  # Remove the loop variable 'i' after the loop completes
        # Combine the new individuals with the calculated breed information
        tmp <- cbind(tmp, aux)
        # Update the current XB population for breed 'b' with the new generation individuals
        pedaux$XB[[b]] <- tmp
        # Clean up temporary variables to free up memory
        rm(tmp, aux)
      };rm(b)  # Remove the loop variable 'b' after the loop completes
      # Convert the list of crossbred populations (XB) into a single data frame
      pedaux$XB <- list2data(pedaux$XB)
      # Update the IDs of the crossbred individuals in a sequential manner
      pedaux$XB$ID <- paste0(substr(pedaux$XB$ID, 1, 8), substr(1000000 + (1:nrow(pedaux$XB)), 2, 7))
      # Combine the new purebred (PB) and crossbred (XB) individuals into the main pedigree (pedtmp)
      pedtmp <- rbind(pedtmp, pedaux$PB, pedaux$XB)
    }
    if (k < startXBGen) {
      pedtmp <- rbind(pedtmp,pedaux$PB)
    }
    if (k > startXBGen) {  # If the current generation is after the start of crossbreeding
      # Update the row names of the XBF2 population to ensure they are unique
      rownames(pedaux$XBF2) <- 1:nrow(pedaux$XBF2)
      # Simulate the next generation of F2 crossbred individuals
      tmp <- simPopInfo(
        c(nrow(pedaux$XBF2), size_off_XBF2),
        n.gen = 1,                  # Simulate only 1 generation after the current one
        pop.info = pedaux$XBF2[, 1:5],  # Use the current XBF2 population as the base for new generation
        pre.id = "XF2"              # Prepend the ID with "XF2" for the crossbred generation
      )
      # Remove the individuals from the base population (founders) and keep only the new generation
      tmp <- tmp[-(1:nrow(pedaux$XBF2)),]
      # Adjust the ID of the individuals in the new generation by appending the generation number to the ID
      tmp$ID <- paste0(
        substr(tmp$ID, 1, 4),       # Keep the first part of the ID
        ifelse(k < 10, "00", ifelse(k < 100, "0", "")),  # Add leading zeros for generation number (for proper formatting)
        k,                          # Add the generation number to the ID
        substr(tmp$ID, 8, 14)        # Keep the remaining part of the ID
      )
      
      # Add a new column to label these individuals as 'cross' breed
      tmp <- cbind(tmp, breed = rep("cross", nrow(tmp)))
      
      # Create a data frame to store the breed information (columns B1, B2, etc.)
      aux <- as.data.frame(matrix(0, nrow(tmp), ncol(pedaux$XBF2) - ncol(tmp)))
      names(aux) <- paste0("B", 1:ncol(aux))  # Name the columns as B1, B2, etc.
      # For each individual in the new generation, calculate the average breed information from the parents
      for (i in 1:nrow(tmp)) {
        # Extract the breed information (columns B1, B2, etc.) from the parents and calculate the mean for the offspring
        aux[i, ] <- apply(
          pedaux$XBF2[pedaux$XBF2$ID %in% c(tmp$parent1[i], tmp$parent2[i]), 7:ncol(pedaux$XBF2)], 2, mean
        )
      }; rm(i)  # Remove the loop index 'i' after the loop completes
      # Combine the new individuals with the calculated breed information
      tmp <- cbind(tmp, aux)
      # Update the F2 crossbred population (pedaux$XBF2) with the new generation data
      pedaux$XBF2 <- tmp
      # Remove temporary variables to free up memory
      rm(tmp, aux)
      # Add the new F2 crossbred individuals to the main pedigree table (pedtmp)
      pedtmp <- rbind(pedtmp, pedaux$XBF2)
    }
    rm(pedaux,list2data)
    if(verbose){
      message("DONE")
    }
    # new genotypes according to matings defined
    if(verbose){
      message("  Genotypes ... ",appendLF=FALSE)
    }
    
    gc()
    if (k == 1) {
      geno_tmp <- geno$BOA  # Initialize geno_tmp for the first generation
    } else {
      geno_tmp <- combined_BOA  # Use combined_BOA for subsequent generations
    }
    
    # Filter geno_tmp to include only individuals in pedtmp
    geno_tmp <- list(
      geno_tmp[[1]][rownames(geno_tmp[[1]]) %in% pedtmp[pedtmp$generation == 0, ]$ID, ],
      geno_tmp[[2]][rownames(geno_tmp[[2]]) %in% pedtmp[pedtmp$generation == 0, ]$ID, ]
    )
    # Generate genotypes for the current generation
    tmp <- list(numeric(0), numeric(0))  # Initialize tmp
    for (chr in 1:max(map$chrom)) {
      aux <- popBreed(
        list(geno_tmp[[1]][, map$chrom == chr], geno_tmp[[2]][, map$chrom == chr]), 
        pedtmp
      )
      # Extract genotypes for individuals in the current generation
      tmp[[1]] <- cbind(tmp[[1]], aux[[1]][pedtmp$generation == 1, ])
      tmp[[2]] <- cbind(tmp[[2]], aux[[2]][pedtmp$generation == 1, ])
      rm(aux); gc()  # Clean up temporary data
    };rm(chr)
    
    rownames_newgen <- pedtmp$ID[pedtmp$generation == 1]
    # Assign row names to tmp based on pedtmp IDs
    rownames(tmp[[1]]) <- rownames_newgen
    rownames(tmp[[2]]) <- rownames_newgen
    
    # Combine genotypes into combined_BOA
    if (k == 1) {
      combined_BOA <- list(
        rbind(geno$BOA[[1]], tmp[[1]]),  # Combine first matrices
        rbind(geno$BOA[[2]], tmp[[2]])   # Combine second matrices
      )
    } else {
      # Append new genotypes to combined_BOA for subsequent generations
      combined_BOA <- list(
        rbind(combined_BOA[[1]], tmp[[1]]),  # Combine first matrices
        rbind(combined_BOA[[2]], tmp[[2]])   # Combine second matrices
      )
    }
    
    geno$BOA <- tmp
    rm(tmp)
    pedtmp <- pedtmp[pedtmp$generation == 1,]
    pedtmp$generation <- k
    ped <- rbind(ped,pedtmp[pedtmp$generation == k,])
    rownames(ped) <- 1:nrow(ped)
    # rm(pedtmp)
    genoQTL <- list(subset(geno$BOA[[1]],select=((1:nrow(map)) %in% geno$QTL)),subset(geno$BOA[[2]],select=((1:nrow(map)) %in% geno$QTL)))
    pbaseQTL <- map[geno$QTL,-(1:3)]
    if(k == 1){
      breedProp <- list(PED=ped[,-c(2:6)],GENO=ped[ped$generation == 0,-c(2:6)],QTL=ped[ped$generation == 0,-c(2:6)])
    } else breedProp$PED <- rbind(breedProp$PED,ped[ped$generation == k,-c(2:6)])
    tmp <- list(GENO=numeric(0),QTL=numeric(0))
    for(b in 1:(ncol(ped)-6)){
      tmp$GENO <- cbind(tmp$GENO,(apply(geno$BOA[[1]]%/%10 == b,1,sum) + apply(geno$BOA[[2]]%/%10 == b,1,sum))/(2*ncol(geno$BOA[[1]])))
      tmp$QTL <- cbind(tmp$QTL,(apply(genoQTL[[1]]%/%10 == b,1,sum) + apply(genoQTL[[2]]%/%10 == b,1,sum))/(2*length(geno$QTL)))
    };rm(b)
    tmp$GENO <- cbind(ped$ID[ped$generation == k],as.data.frame(tmp$GENO));names(tmp$GENO) <- names(ped[,-(2:6)])
    tmp$QTL <- cbind(ped$ID[ped$generation == k],as.data.frame(tmp$QTL));names(tmp$QTL) <- names(ped[,-(2:6)])
    breedProp$GENO <- rbind(breedProp$GENO,tmp$GENO)
    breedProp$QTL <- rbind(breedProp$QTL,tmp$QTL)
    rm(tmp)
    rownames(breedProp$GENO) <- 1:nrow(breedProp$GENO)
    rownames(breedProp$QTL) <- 1:nrow(breedProp$QTL)
    
    if(verbose){
      message("DONE")
      message("  Saving simulated data ... ",appendLF=FALSE)
    }
    # save new pedigree data
    save(ped,breedProp,file=paste0(BasePop$path,"/ped",ifelse(k < 10,"0",""),k,"_Ratio",varDVarGratio,".Rdata"))
    # save new genotype data
    save(map,geno,file=paste0(BasePop$path,"/geno",ifelse(k < 10,"0",""),k,"_Ratio",varDVarGratio,".Rdata"))
    if(verbose){
      message("DONE")
      message("  Phenotypes ... ",appendLF=FALSE)
    }
    #------------#
    # phenotypes #
    #------------#
    BOApartialBV <- function(x) {
      # Compute Mtmp and Mboa only once
      Mtmp <- cbind(
        genoQTL[[1]][, x] %% 10 - as.numeric(pbaseQTL[x, genoQTL[[1]][, x] %/% 10]),
        genoQTL[[2]][, x] %% 10 - as.numeric(pbaseQTL[x, genoQTL[[2]][, x] %/% 10])
      )
      Mboa <- cbind(genoQTL[[1]][, x] %/% 10, genoQTL[[2]][, x] %/% 10)
      Mgeno <- genoQTL[[1]][,x]-10*genoQTL[[1]][,x]%/%10+genoQTL[[2]][,x]-10*genoQTL[[2]][,x]%/%10
      # Initialize gQTL and placeholders for effects
      num_individuals <- nrow(Mtmp)
      num_traits <- ncol(pheno$g)
      
      additive_effects <- matrix(0, num_individuals, num_traits)
      heterosis_effects <- matrix(0, num_individuals, num_traits)
      gQTL <- matrix(0, num_individuals, num_traits)
      
      # Compute additive effects
      for (j in 1:num_traits) {
        atmp <- matrix(0, num_individuals, ncol(Mtmp))
        for (b in seq_along(pheno$a)) {
          atmp[Mboa == b] <- pheno$a[[b]][x, j]
        }
        additive_effects[, j] <- rowSums(Mtmp * atmp)
      }
      # Compute heterosis effects
      for (j in 1:num_traits) {
        # Identify individuals with heterozygous genotype (Mgeno == 1)
        heterozygous <- which(Mgeno == 1)
        # Among heterozygous individuals, check if Mboa columns are different
        heterosis_condition <- heterozygous[Mboa[heterozygous, 1] != Mboa[heterozygous, 2]]
        # Add heterosis effect to the corresponding individuals
        heterosis_effects[heterosis_condition, j] <- pheno[["h"]][x, j]
      }
      
      # Combine additive and heterosis effects
      rownames(additive_effects) <- rownames_newgen
      rownames(heterosis_effects) <- rownames_newgen
      
      return(list(additive = additive_effects, heterosis = heterosis_effects))
    }
    # Initialize gNEW with the first result of BOApartialBV
    gNEW <- BOApartialBV(1)
    
    # Loop over the remaining rows of pbaseQTL and accumulate results
    for (j in 2:nrow(pbaseQTL)) {
      result <- BOApartialBV(j)
      # Accumulate each component separately
      gNEW$additive <- gNEW$additive + result$additive
      gNEW$heterosis <- gNEW$heterosis + result$heterosis
    };rm(j)
    pheno$varDVarGratio = varDVarGratio
    # Identify rows that start with "XF"
    xf_rows <- grepl("^XF", rownames(gNEW$additive))
    # If any such rows exist, apply the scaling transformation
    if (any(xf_rows)) {
      # Compute varD for selected individuals
      varD <- (pheno$varDVarGratio * colVars(gNEW$additive[xf_rows, , drop = FALSE])) / 
        (1 - pheno$varDVarGratio)
      # Compute scaling factors
      scaling_factors <- sqrt(varD / colVars(gNEW$heterosis[xf_rows, , drop = FALSE]))
      # Scale only the selected rows, keeping the rest unchanged
      gNEW$heterosis[xf_rows, ] <- sweep(gNEW$heterosis[xf_rows, , drop = FALSE], 2, scaling_factors, "*")
    }
    gNEW$total <- gNEW$additive + gNEW$heterosis
    eNEW <- matrix(mvrnorm(n=sum(ped$generation == k),mu=rep(0,ncol(pheno$rho_e)),Sigma=pheno$rho_e),
                   sum(ped$generation == k),ncol(pheno$rho_e))
    eNEW <- matrix(sqrt(pheno$Vy*(1-pheno$h2)),sum(ped$generation == k),
                   ncol(eNEW),byrow=TRUE)*eNEW/matrix(apply(eNEW,2,sd),sum(ped$generation == k),ncol(eNEW),byrow=TRUE)
    yNEW <- gNEW$total+eNEW+matrix(pheno$mu,sum(ped$generation == k),ncol(gNEW$total),byrow=TRUE)
    
    rownames(yNEW) <- rownames_newgen
    pheno$g_add <- rbind(pheno$g_add,gNEW$additive)
    pheno$g_h <- rbind(pheno$g_h,gNEW$heterosis)
    pheno$g <- rbind(pheno$g,gNEW$total)
    pheno$e <- rbind(pheno$e,eNEW)
    pheno$y <- rbind(pheno$y,yNEW)
    rm(genoQTL,pbaseQTL,gNEW,eNEW,yNEW,BOApartialBV,rownames_newgen)
    if(verbose){
      message("DONE")
      message("  Saving simulated data ... ",appendLF=FALSE)
    }
    # save new phenotype data
    save(pheno,file=paste0(BasePop$path,"/pheno",ifelse(k < 10,"0",""),k,"_Ratio",varDVarGratio,".Rdata"))
    if(verbose){
      message("DONE\n")
    }
  };rm(k)
  if(verbose){
    message("DONE\n")
    message("|-------------------------------------------------------------|")
    message("| Summary of created files                                    |")
    message("|-------------------------------------------------------------|\n")
    message("Folder in which generated files are saved:")
    message("  ",BasePop$path,"\n")
    message("Files created:")
    for(k in genfiles){
      message("  geno",ifelse(k < 10,"0",""),k,"_Ratio",varDVarGratio,".Rdata")
    };rm(k)
    for(k in genfiles){
      message("  ped",ifelse(k < 10,"0",""),k,"_Ratio",varDVarGratio,".Rdata")
    };rm(k)
    for(k in genfiles){
      message("  pheno",ifelse(k < 10,"0",""),k,"_Ratio",varDVarGratio,".Rdata")
    };rm(k)
  }
  tmp <- list(dir(BasePop$path)[substr(dir(BasePop$path),1,4) == "geno" & dir(BasePop$path) != BasePop$files$geno & nchar(dir(BasePop$path)) == 12])
  tmp[[2]] <- dir(BasePop$path)[substr(dir(BasePop$path),1,3) == "ped" & dir(BasePop$path) != BasePop$files$ped & nchar(dir(BasePop$path)) == 11]
  tmp[[3]] <- dir(BasePop$path)[substr(dir(BasePop$path),1,5) == "pheno" & dir(BasePop$path) != BasePop$files$pheno & nchar(dir(BasePop$path)) == 13]
  gc()
  return(list(path=BasePop$path,files=list(geno=tmp[[1]],ped=tmp[[2]],pheno=tmp[[3]])))
}

split_snps <- function(x) {
  # Ensure x is character, then split
  first_digit <- as.numeric(substr(as.character(x), 1, 1))
  second_digit <- as.numeric(substr(as.character(x), 2, 2))
  return(list(BOA = first_digit, phase = second_digit))
}
# Function to split an entire dataframe into BOA and phase and return numeric results
split_BOA_phase <- function(BOA_Phase_df) {
  # Initialize empty data frames for BOA and phase
  BOA <- data.frame(matrix(ncol = ncol(BOA_Phase_df), nrow = nrow(BOA_Phase_df)))
  phase <- data.frame(matrix(ncol = ncol(BOA_Phase_df), nrow = nrow(BOA_Phase_df)))
  
  # Name the columns the same as original dataframe
  colnames(BOA) <- colnames(BOA_Phase_df)
  colnames(phase) <- colnames(BOA_Phase_df)
  
  # Iterate over SNP columns and split values
  for (i in 1:ncol(BOA_Phase_df)) {
    snp_split <- lapply(BOA_Phase_df[, i], split_snps)
    
    # Extract the first and second elements from each split SNP
    BOA[[i]] <- sapply(snp_split, function(x) x$BOA)
    phase[[i]] <- sapply(snp_split, function(x) x$phase)
  }
  
  # Return both BOA and phase as a list
  return(list(BOA = BOA, Phase = phase))
}

split_save <- function(varDVarGratio, out_path, stratgen, endgen){
  for (i in stratgen:endgen) {
    print(paste('Generation', i))
    if (i < 10) { 
      load(paste0(out_path, 'geno0', i, '_Ratio',varDVarGratio,'.Rdata'))
    } else { 
      load(paste0(out_path, 'geno', i, '_Ratio',varDVarGratio,'.Rdata'))
    }
    print('split phase1')
    res_BOA_phase_1 <- split_BOA_phase(as.data.frame(geno$BOA[[1]]))
    print('split phase2')
    res_BOA_phase_2 <- split_BOA_phase(as.data.frame(geno$BOA[[2]]))
    print('geno generation')
    geno_gen <- res_BOA_phase_1[['Phase']] + res_BOA_phase_2[['Phase']]
    
    print('geno save')
    write.table(as.data.frame(geno_gen), 
                file = paste0(out_path, "geno_gen", i,'_Ratio',varDVarGratio,".txt"), 
                row.names = TRUE, col.names = TRUE)
    print('phase1 save')
    write.table(as.data.frame(res_BOA_phase_1[['Phase']]), 
                file = paste0(out_path, "phase1_gen", i,'_Ratio',varDVarGratio,".txt"), 
                row.names = TRUE, col.names = TRUE)
    print('phase2 save')
    write.table(as.data.frame(res_BOA_phase_2[['Phase']]), 
                file = paste0(out_path, "phase2_gen", i,'_Ratio',varDVarGratio,".txt"), 
                row.names = TRUE, col.names = TRUE)
    print('BOA1 save')
    write.table(as.data.frame(res_BOA_phase_1[['BOA']]), 
                file = paste0(out_path, "BOA1_gen", i,'_Ratio',varDVarGratio,".txt"), 
                row.names = TRUE, col.names = TRUE)
    print('BOA2 save')
    write.table(as.data.frame(res_BOA_phase_2[['BOA']]), 
                file = paste0(out_path, "BOA2_gen", i,'_Ratio',varDVarGratio,".txt"), 
                row.names = TRUE, col.names = TRUE)
  }
  print('Saving Done')
}

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

read_and_combine <- function(prefix, generations, path, ratio) {
  file_list <- lapply(generations, function(i) {
    as.matrix(fread(paste0(path, prefix, "_gen", i, "_Ratio", ratio, ".txt"))[, -1])
  })
  do.call(rbind, file_list)
}


save_pheno_data <- function(BasePop, lastgen, varDVarGratio, output_folder) {
  # Load data
  PopDataEvol <- LoadPop(BasePop, lastgen, varDVarGratio, verbose=TRUE)
  
  alpha_B1 <- as.data.frame(PopDataEvol[["pheno"]][["a"]][["breed1"]])
  alpha_B2 <- as.data.frame(PopDataEvol[["pheno"]][["a"]][["breed2"]])
  alpha_B1$QTLpos <- as.numeric(PopDataEvol[["geno"]][["QTL"]])
  alpha_B2$QTLpos <- as.numeric(PopDataEvol[["geno"]][["QTL"]])
  
  # List of common data frames
  data_frames <- list(
    alpha_B1 = alpha_B1,
    alpha_B2 = alpha_B2,
    ped_gen = PopDataEvol[["ped"]],
    gv_gen = PopDataEvol[["pheno"]][["g"]],
    pheno_gen = as.data.frame(PopDataEvol[["pheno"]][["y"]]),
    breedProp_Ped = as.data.frame(PopDataEvol[["breedProp"]][["PED"]]),
    breedProp_Geno = as.data.frame(PopDataEvol[["breedProp"]][["GENO"]]),
    breedProp_QTL = as.data.frame(PopDataEvol[["breedProp"]][["QTL"]]),
    map_df = as.data.frame(PopDataEvol[["map"]])
    
  )
  
  # Add heterosis effects and gv_hetero_gen only if varDVarGratio is greater than 0
  if (varDVarGratio > 0) {
    data_frames$Heterosis_effects <- as.data.frame(PopDataEvol[["pheno"]][["h"]])
    data_frames$gv_hetero_gen <- as.data.frame(PopDataEvol[["pheno"]][["g_h"]])
    data_frames$gv_add_gen <- as.data.frame(PopDataEvol[["pheno"]][["g_add"]])
  }
  
  # Save each data frame as a CSV file
  for (name in names(data_frames)) {
    file_path <- paste0(output_folder, name, "_Ratio", varDVarGratio, ".csv")
    write.csv(data_frames[[name]], file_path, row.names = TRUE)
    
    # Optional: print status to monitor progress
    print(paste("Saved", name, "to", file_path))
  }
}

run_bayes_analysis <- function(varDVarGratio, out_path, startgenTrain=21, endgenTrain=29) {
  map_df <- read.csv(paste0(out_path, "map_df_Ratio", varDVarGratio, ".csv"))
  alpha_B1 <- as.data.frame(read.csv(paste0(out_path, "alpha_B1_Ratio", varDVarGratio, ".csv"), row.names=1))
  pheno <- read.csv(paste0(out_path, "pheno_gen_Ratio", varDVarGratio, ".csv"), row.names=1)
  ped <- read.csv(paste0(out_path, "ped_gen_Ratio", varDVarGratio, ".csv"), row.names=1)
  
  p1 <- map_df$AlleleFreq_basePB1
  p2 <- map_df$AlleleFreq_basePB2
  
  geno <- read_and_combine("geno", startgenTrain:endgenTrain, out_path, varDVarGratio)
  BOA1 <- read_and_combine("BOA1", startgenTrain:endgenTrain, out_path, varDVarGratio)
  phase1 <- read_and_combine("phase1", startgenTrain:endgenTrain, out_path, varDVarGratio)
  BOA2 <- read_and_combine("BOA2", startgenTrain:endgenTrain, out_path, varDVarGratio)
  phase2 <- read_and_combine("phase2", startgenTrain:endgenTrain, out_path, varDVarGratio)
  
  BOA <- list(BOA1, BOA2)
  Phase <- list(phase1, phase2)
  
  scaled_geno_B1 <- scale_geno(1, p1, BOA, Phase)
  scaled_geno_B2 <- scale_geno(2, p2, BOA, Phase)
  
  scaled_geno_B1_noQTL <- scaled_geno_B1[, -alpha_B1$QTLpos]
  scaled_geno_B2_noQTL <- scaled_geno_B2[, -alpha_B1$QTLpos]
  
  pheno_train <- pheno[ped[ped$generation >= startgenTrain & ped$generation <= endgenTrain, "ID"], ]
  
  for (i in 1:4) {
    BayesA_train <- BGLR(y=pheno_train[[paste0("V", i)]], saveAt=paste0("car", i, "_BA_all"),
                         ETA=list(list(X=scaled_geno_B1_noQTL, model="BayesA"),
                                  list(X=scaled_geno_B2_noQTL, model="BayesA")),
                         nIter=30000, burnIn=15000)
    SNPeffBA_1 <- BayesA_train$ETA[[1]]$b
    SNPeffBA_2 <- BayesA_train$ETA[[2]]$b
    SNPeff_df <- data.frame(SNPeffBA_1, SNPeffBA_2)
    write.csv(SNPeff_df, paste0(out_path, "SNPeff_trait", i, "_Ratio", varDVarGratio, ".csv"), 
              row.names=FALSE, quote=FALSE)
  }
}

load_gblup_data <- function(out_path, varDVarGratio, startgen=21, endgen=35) {
  print(paste("Loading GBLUP data for Ratio:", varDVarGratio))
  geno_list <- list()
  for (i in startgen:endgen) {
    print(paste("Loading genotype data for generation:", i))
    geno_list[[as.character(i)]] <- as.matrix(fread(paste0(out_path, "geno_gen", i, "_Ratio", varDVarGratio, ".txt"))[, -1])
  }
  geno <- do.call(rbind, geno_list)
  
  print("Loading alpha and phenotype data")
  alpha_B1 <- as.data.frame(read.csv(paste0(out_path, "alpha_B1_Ratio", varDVarGratio, ".csv"), row.names=1))
  geno1 <- geno[, -alpha_B1$QTLpos]
  G1 <- mkGRM(geno1)
  
  pheno <- read.csv(paste0(out_path, "pheno_gen_Ratio", varDVarGratio, ".csv"), row.names=1)
  ped <- read.csv(paste0(out_path, "ped_gen_Ratio", varDVarGratio, ".csv"), row.names=1)
  
  print("GBLUP data loading completed")
  list(G1=G1, Geno=geno1, pheno=pheno, ped=ped)
}

run_purebred_analysis <- function(G1, Geno, pheno, ped, out_path, varDVarGratio, 
                                  breed_indicator, startgen=21, endgen=35, testgenstart=30) {
  print(paste("Running purebred analysis for breed:", breed_indicator, "and Ratio:", varDVarGratio))
  
  pheno_train_test <- pheno[ped[ped$generation >= startgen & ped$generation <= endgen, "ID"], ]
  ids_gen_test <- ped[ped$generation >= testgenstart & ped$generation <= endgen, "ID"]
  ids_cross <- ped[ped$breed == "cross", "ID"]
  ids_breed <- ped[ped[[breed_indicator]] == 0, "ID"]
  
  ids_to_replace <- unique(c(ids_gen_test, ids_cross, ids_breed))
  pheno_train_test[rownames(pheno_train_test) %in% ids_to_replace, ] <- NA
  
  print("Running Multitrait model")
  MT_PB <- Multitrait(y=as.matrix(pheno_train_test), ETA=list(list(K=G1, model="RKHS")), 
                      nIter=30000, burnIn=15000, thin=10, verbose=TRUE)
  
  print("Saving results")
  mu_values <- MT_PB[["mu"]]
  write.csv(mu_values, paste0(out_path, "mu_MTGBLUP_", breed_indicator, "_Ratio", varDVarGratio, ".csv"), 
            row.names=FALSE, quote=FALSE)
  
  PBV <- as.data.frame(sapply(1:ncol(MT_PB$ETA[[1]]$u), function(i) MT_PB$ETA[[1]]$u[, i]))
  write.csv(PBV, paste0(out_path, "PGV_MTGBLUP_", breed_indicator, "_Ratio", varDVarGratio, ".csv"), 
            row.names=FALSE, quote=FALSE)
  
  
  print("Extract SNP effects Started")
  # Filter the pedigree to exclude rows with unwanted IDs
  ped_filtered <- ped[!ped$ID %in% ids_to_replace & ped$generation >= 21, ]
  rownames(ped_filtered) <- NULL
  
  # Use the row indices of the filtered pedigree to filter geno
  geno_train <- Geno[as.numeric(rownames(ped_filtered)),]
  # center genotype matrix
  genoCentered <- scale(geno_train,center=TRUE,scale=FALSE)
  
  invM = ginv(genoCentered)
  
  PGV_train <- PBV[as.numeric(rownames(ped_filtered)),]
  
  SNPeff <- invM %*% as.matrix(PGV_train)
  write.csv(SNPeff, paste0(out_path, "GBLUP_SNPeff_", breed_indicator, "_Ratio", varDVarGratio, ".csv"),
            row.names=FALSE, quote=FALSE)
  
  print("Purebred analysis completed")
}

run_gblup_all_breeds <- function(G1, pheno, ped, out_path, varDVarGratio, startgen=21, endgen=35, testgenstart=30) {
  print(paste("Running GBLUP for all breeds with Ratio:", varDVarGratio))
  pheno_train_test <- pheno[ped[ped$generation >= startgen & ped$generation <= endgen, "ID"],]
  ids_gen_test <- ped[ped$generation >= testgenstart & ped$generation <= endgen, "ID"]
  pheno_train_test[rownames(pheno_train_test) %in% ids_gen_test, ] <- NA
  
  print("Running Multitrait model for all breeds")
  MT <- Multitrait(y=as.matrix(pheno_train_test), ETA=list(list(K=G1, model="RKHS")), 
                   nIter=30000, burnIn=15000, thin=10, verbose=TRUE)
  
  PBV <- as.data.frame(sapply(1:ncol(MT$ETA[[1]]$u), function(i) MT$ETA[[1]]$u[, i]))
  write.csv(PBV, paste0(out_path, "PGV_MTGBLUP_PB12CB_Ratio", varDVarGratio, ".csv"), row.names=FALSE, quote=FALSE)
  print("GBLUP for all breeds completed")
}

run_gblup_cross_bred <- function(G1, pheno, ped, out_path, varDVarGratio, startgen=21, endgen=35, testgenstart=30) {
  print(paste("Running GBLUP for crossbred with Ratio:", varDVarGratio))
  pheno_train_test <- pheno[ped[ped$generation >= startgen & ped$generation <= endgen, "ID"],]
  ids_gen_test <- ped[ped$generation >= testgenstart & ped$generation <= endgen, "ID"]
  ids_b1 <- ped[ped$B1 == 1, "ID"]
  ids_b2 <- ped[ped$B2 == 1, "ID"]
  ids_to_replace <- unique(c(ids_gen_test, ids_b1, ids_b2))
  pheno_train_test[rownames(pheno_train_test) %in% ids_to_replace, ] <- NA
  
  print("Running Multitrait model for crossbred")
  MT_CB <- Multitrait(y=as.matrix(pheno_train_test), ETA=list(list(K=G1, model="RKHS")), 
                      nIter=30000, burnIn=15000, thin=10, verbose=TRUE)
  
  PBV <- as.data.frame(sapply(1:ncol(MT_CB$ETA[[1]]$u), function(i) MT_CB$ETA[[1]]$u[, i]))
  write.csv(PBV, paste0(out_path, "PGV_MTGBLUP_CB_Ratio", varDVarGratio, ".csv"), row.names=FALSE, quote=FALSE)
  print("GBLUP for crossbred completed")
}


# ## Sim Base pop
rho_g <- diag(0.5,4)
rho_g[lower.tri(rho_g,diag=FALSE)] <- c(0.6,-0.2,0.15,-0.1,0.4,0.3)
rho_g <- rho_g+t(rho_g)

namedir <- '10k_500QTL_1000nPB_rep1'
BasePop <- SimBasePop(n=1000,m=10000,nQTL=500,mu=c(0,0,0,0),Vy=c(1,1,1,1),
                      h2=c(0.3,0.05,0.1,0.5),nbreed=2,corQTLbreed=0.3,
                      rho_g=rho_g,rho_e=diag(1,4), namedir=namedir)

PopData <- LoadPop(BasePop,0,verbose=TRUE)
summary(PopData)


out_path = paste0(getwd(),'/', namedir,'/')
### Evolve Pop
varDVarG_ratios <- seq(0.1, 0.3, by=0.1)

for (varDVarGratio in varDVarG_ratios) {
  EvolvePop(PopData, ngen=35, startXBGen=25, selTraitPB=c(1,3),
            selTraitXB=0, selIntensityPBFemales=c(1,1),
            selIntensityPBMales=c(0.5, 0.65), selIntensityXB=0.5,
            varDVarGratio=varDVarGratio)

  save_pheno_data(BasePop, 35, varDVarGratio, out_path)
  split_save(varDVarGratio, out_path, 20, 35)
}

### Run Bayes A
for (varDVarGratio in varDVarG_ratios) {
  run_bayes_analysis(varDVarGratio, out_path)
}

### Run GBLUP 

breed_indicators <- c("B1", "B2")

for (varDVarGratio in varDVarG_ratios) {
  print(paste("Processing Ratio:", varDVarGratio))
  data <- load_gblup_data(out_path, varDVarGratio)
  
  for (breed_indicator in breed_indicators) {
    print(paste("Running purebred analysis for breed:", breed_indicator, "with Ratio:", varDVarGratio))
    run_purebred_analysis(data$G1, data$Geno, data$pheno, data$ped, out_path, varDVarGratio, breed_indicator)
  }
  
  print(paste("Running GBLUP for all breeds with Ratio:", varDVarGratio))
  run_gblup_all_breeds(data$G1, data$pheno, data$ped, out_path, varDVarGratio)
  
  print(paste("Running GBLUP for crossbred with Ratio:", varDVarGratio))
  run_gblup_cross_bred(data$G1, data$pheno, data$ped, out_path, varDVarGratio)
  
  print(paste("Completed processing for Ratio:", varDVarGratio))
}






