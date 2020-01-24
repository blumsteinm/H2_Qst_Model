###############################################################################
###############################################################################
## 
## 11/30/17
##
## Written by: Meghan Blumstein (blumsteinm@gmail.com)
##
## Calculate Heritability and Qst of NSC data
##
## Input: NSC Data by genoytpe and population (NSC_Dataset_Step3...)
##
## Output: H2 matrix and Qst matrix; NSC Data Step4 - Includes the Population and Genotype IDS
##
## ***Note for all models, raw data is gamma distributed (long right tail) and run as such. For this analysis, it is as important to
##  replicate the patterns of the tails of the distribution as finding the mean, since we are looking for genotypes/populations that diverge from the group.
##  If the minimum value for any set of values was 0, 0.01 was added instead of running a hurdle model, as there isn't much of a biological difference between 
##  0 and 0.01 mg/g of sugar/starch.
###############################################################################
###############################################################################
rm(list = ls())
###############################################################################
##                              SETUP
###############################################################################

## Packages
require(plyr)
library(lme4)
require(rstanarm)
require(rstan)

## Directories
dropbox_folder <- "/Users/meghs/Dropbox/PhD_Dissertation/"
# dropbox_folder <- "C:/Users/BlumsteinMeghan/Dropbox/PhD_Dissertation/"
wd <- paste0(dropbox_folder, "Data_Files/NSC_Data/")
od <- paste0(dropbox_folder, "H2_Qst/")
setwd(wd)

## Read in data
nsc_data <- read.csv("NSC_Dataset_Step3_01.14.19.csv")

## Functions
source("/Users/Meghs/Dropbox/Code_From_HF/Favorite_Functions.R")
# source("C:/Users/BlumsteinMeghan/Dropbox/Code_From_HF/Favorite_Functions.R")
stan_code_qst <- paste0(dropbox_folder, "Code/H2_Qst_Plasticity/Gamma_Hierarchical_no_predictors_qst.stan")
stan_code_h2 <- paste0(dropbox_folder, "Code/H2_Qst_Plasticity/Gamma_Hierarchical_no_predictors_h2.stan")
stan_code_h2_normal <- paste0(dropbox_folder, "Code/Other/STAN_Gamma_Hierarchical/Normal_Hierarchical_01.09.2018.stan")

## Constants
num_cores <- 6
num_iter <- 6000

## Generate a third replicate?
rep3 <- F

###############################################################################
##                     SETUP ITERATORS  
###############################################################################

## Iteration Values
gardens <- c("Clatskanie")
tissues <- c("Roots", "Stems")
trait_columns <- c("TNC")

#############################################################################
##                    Calculate Heritability and Qst
###############################################################################
for(j in 1:length(gardens)){
     
     for(k in 1:length(tissues)){
          
          ## Only run branches in corvallis, not other tissues
          if(gardens[j] == "Corvallis" & tissues[k] != "Branches"){next}
          # if(gardens[j] == "Clatskanie" & tissues[k] == "Branches"){next} ## If you want to skip Clatskanie Branches
          
          ## Subset to tissue and garden of outer loop
          sub_data <- nsc_data[nsc_data$Garden == gardens[j] & nsc_data$Tissue == tissues[k],]
          
          ## Remove Populations with 2 or less genotypes
          remove_pops <- ddply(sub_data, .(Population), summarize, n = length(unique(Genotype)))
          rm_pops_list <- as.character(remove_pops$Population[remove_pops$n < 3])
          sub_data <- sub_data[!sub_data$Population %in% rm_pops_list,]
          
          ## Reconvert columns from factor to character - do this so factor levels are reset
          sub_data <- unfactor_df(sub_data)
          sub_data <- factor_df(sub_data)
          
          ## Add a numeric code for Populations and Genotypes
          identifiers <- ddply(sub_data, .(Population, Genotype), summarize, 
                               River_ID  = as.numeric(Population)[1],
                               Genotype_ID = as.numeric(Genotype)[1])
          sub_data <- merge(sub_data, identifiers)
          
          ## Save out the data used for the analysis
          nme <- paste0(gardens[j], "_", tissues[k])
          write.csv(sub_data,  paste0(od, nme, "_values_", format(Sys.Date(), "%d.%m.%y"), ".csv"), row.names = F)
          
          for(i in 1:length(trait_columns)){
               
               ############################
               ## Data initialization
               ############################
               
               ## Indicate Iteration
               nme <- paste0(gardens[j], "_", tissues[k], "_", trait_columns[i])
               print(nme)
               
               ## Set y value of models
               y <- sub_data[,trait_columns[i]]

               ## Add 0.01 if minimum of y is 0
               if(min(y) <= 0){y <- y + 0.01}

               ##################################
               ## Format Data for Bayesian Runs
               ##################################

               ## Set up looping structure for stan model - population ID must be in genotype (1:ngen) order
               sub_level_b <- ddply(identifiers, .(River_ID, Genotype_ID), summarize, n = 1)
               pops2gens <- data.frame(River_ID = sub_data$River_ID, Genotype_ID = sub_data$Genotype_ID)
               pops2gens <- pops2gens[order(pops2gens$Genotype_ID),]

               ## Format for Stan - H2
               bayesData_H2 <- list(Ni = nrow(sub_data),
                                    Ng = length(unique(sub_data$Genotype)),
                                    genotype_ids = sub_data$Genotype_ID,
                                    Y_ig = y)
               
               ## Format for Stan - Qst
               bayesData_Qst <- list(Ni = nrow(sub_data),
                                    Ng = length(unique(sub_data$Genotype)),
                                    Np = length(unique(sub_data$Population)),
                                    genotype_ids = sub_data$Genotype_ID,
                                    population_ids = pops2gens$River_ID,
                                    PopsToGens = sub_level_b$River_ID,
                                    Y_igp = y)

               ############################
               ## Heritability & Qst
               ############################
               
               ## Output filenames
               datetime <-  format(Sys.time(), "%m.%d.%y")
               of1 <- paste0(od, "Model_Outputs/", nme, "_method_2_bayes_h2_", datetime, ".rds")
               of2 <- paste0(od, "Model_Outputs/", nme, "_method_2_bayes_qst_", datetime, ".rds")
               of3 <- paste0(od, "Model_Outputs/", nme, "_method_2_bayes_h2_normal_", datetime, ".rds")

               ## Method 1  - Evans et al. 2014 H2 y_{ij} = beta0 + alpha_{j} +  e_{ij}, fit bayesian
               m2a <- stan(file = stan_code_h2, data = bayesData_H2, chains = num_cores, cores = num_cores, iter = num_iter, control = list(adapt_delta = 0.99, max_treedepth = 15))
               saveRDS(m2a, file = of1)

               ## Method 2  - Evans et al. 2014 Qst y_{ij} = beta0 + alpha_{j} + e_{ijk}, fit bayesian
               m2b <- stan(file = stan_code_qst, data = bayesData_Qst, chains = num_cores, cores = num_cores, iter = num_iter, control = list(adapt_delta = 0.99, max_treedepth = 15))
               saveRDS(m2b, file = of2)
               
               ## Normal Comparison 
               m2c <- stan_lmer(y ~ (1|Genotype), data = sub_data)
               saveRDS(m2c, file = of3)


          }
     }
}

#############################################################################
##                   COMBINE DATSETS AND SAVE OUT
###############################################################################

final_nsc_data <- c()
for(j in 1:length(gardens)){

     for(k in 1:length(tissues)){

          ## Only run branches in corvallis, not other tissues
          if(gardens[j] == "Corvallis" & tissues[k] != "Branches"){next}
          # if(gardens[j] == "Clatskanie" & tissues[k] == "Branches"){next} ## If you want to skip Clatskanie Branches


               ############################
               ## Data initialization
               ############################

               ## Indicate Iteration
               nme <- paste0(gardens[j], "_", tissues[k])
               print(nme)

               ## Pull in file to add to dataframe & add
               dat <- read.csv(paste0(od, nme, "_values_", format(Sys.Date(), "%d.%m.%y"), ".csv"))
               final_nsc_data <- rbind(final_nsc_data, dat)

               ## Delete File
               file.remove(paste0(od, nme, "_values_", format(Sys.Date(), "%d.%m.%y"), ".csv"))

     }
}


## Save out new file
write.csv(final_nsc_data, paste0(wd, "NSC_Dataset_Step4_", format(Sys.time(), "%m.%d.%y"), ".csv"), row.names = F)









