###############################################################################
###############################################################################
## 
## 11/30/17
##
## Written by: Meghan Blumstein (blumsteinm@gmail.com)
##
## Calculate Heritability and Qst of NSC data
##
## Input: NSC Data by genoytpe and population
##
## Output: H2 matrix and Qst matrix; create files of h2 and qst estimates 
##         from each posterior distribution draw
##
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
dropbox_folder <- "/Users/Meghs/Dropbox/"
# dropbox_folder <- "D:/Dropbox/"
od <- paste0(dropbox_folder, "PhD_Dissertation/H2_Qst/")

## Functions
source(paste0(dropbox_folder, "Code_From_HF/Favorite_Functions.R"))

## Time stamp of models
datetime <- "01.14.19" #format(Sys.time(), "%m.%d.%y")

###############################################################################
##               Setup storage arrays for outputs of analysis
###############################################################################

## Iteration Values
gardens <- c("Clatskanie")
tissues <- c("Roots", "Stems")
trait_columns <- c("TNC")

## Value Matrix
parameters <- c("Population Variation", "Genetic Variation", "MicroEnv Variation", "Qst_Population", "H2")
cnames <- c("Garden", "Tissue", "Trait", "Parameter", "H2_code", "Variance_h2", "Qst_Code", "Variance_qst")
output_df <- matrix(ncol = length(cnames), nrow = 80) ## Rows == number of gardens/tissues (4) * number of parameters (5) * number traits (3)
colnames(output_df) <- cnames


###############################################################################
##                    Calculate Heritability and Qst
###############################################################################

count <- 1
for(j in 1:length(gardens)){
  
  for(k in 1:length(tissues)){
    
    ## Only run branches in corvallis, not other tissues
    if(gardens[j] == "Corvallis" & tissues[k] != "Branches"){next}
    if(gardens[j] == "Clatskanie" & tissues[k] == "Branches"){next}
  
  for(i in 1:length(trait_columns)){
      
      ## Indicate Iteration
      nme     <- paste0(gardens[j], "_", tissues[k], "_", trait_columns[i])
      print(nme)
 
      ## Load model results
      r_moda <- readRDS(file = paste0(od, "Model_Outputs/", nme, "_method_2_bayes_h2_", datetime, ".rds"))
      r_modb <- readRDS(file = paste0(od, "Model_Outputs/", nme, "_method_2_bayes_qst_", datetime, ".rds"))
      
      ##Get Variance/Covariance Matrix
      mod_vcor_a <- as.data.frame( r_moda )
      mod_vcor_b <- as.data.frame( r_modb )
      
      ###############################################
      ## Calculate Mean H2 and Qst
      ###############################################
      
      rowa <- count
      rowb <- count + (length(parameters)) - 1
      
      ## Identifier Info
      output_df[rowa:rowb, 1] <- gardens[j]
      output_df[rowa:rowb, 2] <- tissues[k]
      output_df[rowa:rowb, 3] <- trait_columns[i]
      output_df[rowa:rowb, 4] <- parameters
      
      ## Parse Variation
      microEnv_variation_a <- mod_vcor_a$phi_ig
      microEnv_variation_b <- mod_vcor_b$phi_igp
      
      genetic_variation_a <- mod_vcor_a$phi_g
      genetic_variation_b <- mod_vcor_b$phi_gp
      
      population_variation_a <- NA
      population_variation_b <- mod_vcor_b$phi_p

      ## Output Parameters
      H2_a <- genetic_variation_a / (genetic_variation_a + microEnv_variation_a)
      H2_b <- genetic_variation_b / (genetic_variation_b + microEnv_variation_b)
      
      Qst_population_a <- NA
      Qst_population_b <- population_variation_b / (population_variation_b + (2*genetic_variation_b) ) #Population Qst - as in Whitlock and Gilbert 2012 Mol Ecol Res)

      ## Fill in Variance Parameters
      output_df[rowa + 0, 4 + 1] <- NA
      output_df[rowa + 0, 4 + 2] <- NA
      output_df[rowa + 0, 4 + 3] <- mean(population_variation_b)
      output_df[rowa + 0, 4 + 4] <- var(population_variation_b)
      
      output_df[rowa + 1, 4 + 1] <- mean(genetic_variation_a)
      output_df[rowa + 1, 4 + 2] <- var(genetic_variation_a)
      output_df[rowa + 1, 4 + 3] <- mean(genetic_variation_b)
      output_df[rowa + 1, 4 + 4] <- var(genetic_variation_b)
      
      output_df[rowa + 2, 4 + 1] <- mean(microEnv_variation_a)
      output_df[rowa + 2, 4 + 2] <- var(microEnv_variation_a)
      output_df[rowa + 2, 4 + 3] <- mean(microEnv_variation_b)
      output_df[rowa + 2, 4 + 4] <- var(microEnv_variation_b)
      
      ## Fill in Output Parameters
      output_df[rowa + 3, 4 + 1] <- NA
      output_df[rowa + 3, 4 + 2] <- NA
      output_df[rowa + 3, 4 + 3] <- mean(Qst_population_b)
      output_df[rowa + 3, 4 + 4] <- var(Qst_population_b)
      
      output_df[rowa + 4, 4 + 1] <- mean(H2_a)
      output_df[rowa + 4, 4 + 2] <- var(H2_a)
      output_df[rowa + 4, 4 + 3] <- mean(H2_b)
      output_df[rowa + 4, 4 + 4] <- var(H2_b)
      
      count <- count + length(parameters)
      
      ###############################################
      ## Plot Credible Intervals 
      ###############################################
      
      #####
      ## Random Effects 
      #####
      pdf(paste0(dropbox_folder, "PhD_Dissertation/Figures/H2_Qst_Plasticity/Credible_Intervals_", tissues[k], ".pdf"))
      
      par(mfrow = c(1, 2), mar = c(3.5, 3.5, 0.1, 0.1), mgp = c(2, 1, 0))
      plot(NULL, ylim = c(1, 5), xlim = c(0, 100), xlab = expression(paste(sigma^2, " 95% Credible Intervals")), ylab = "", yaxt = "n")

      ## Alpha g
      y <- 1
      ag <- quantile(genetic_variation_a, c(.05, .50, .95))
      points(ag[2], y, pch = 16)
      lines(c(ag[1], ag[3]), c(y, y), lwd = 2)
      
      ## Alpha ge
      y <- 2
      eg <- quantile(microEnv_variation_a, c(.05, .50, .95))
      points(eg[2], y, pch = 16)
      lines(c(eg[1], eg[3]), c(y, y), lwd = 2)
      
      ## Alpha p
      y <- 3
      ap <- quantile(population_variation_b, c(.05, .50, .95))
      points(ap[2], y, pch = 16)
      lines(c(ap[1], ap[3]), c(y, y), lwd = 2)
      
      ## Alpha pg
      y <- 4
      apg <- quantile(genetic_variation_b, c(.05, .50, .95))
      points(apg[2], y, pch = 16)
      lines(c(apg[1], apg[3]), c(y, y), lwd = 2)
      
      ## Alpha pge
      y <- 5
      epg <- quantile(microEnv_variation_b, c(.05, .50, .95))
      points(epg[2], y, pch = 16)
      lines(c(epg[1], epg[3]), c(y, y), lwd = 2)
      
      labels <- c(expression(paste(sigma^2, alpha["g"])), 
                  expression(paste(sigma^2, epsilon["g"])), 
                  expression(paste(sigma^2, alpha["p"])), 
                  expression(paste(sigma^2, alpha["pg"])), 
                  expression(paste(sigma^2, epsilon["pg"])))
      axis(side = 2, at = c(1:5), labels = labels, las = 2)
      
      #####
      ## H2/Qst
      #####
      
      plot(NULL, ylim = c(1, 4), xlim = c(-0.1, 1), xlab = expression(paste(sigma^2, " 95% Credible Intervals")), ylab = "", yaxt = "n")
      
      ## Heritability 
      h2 <- quantile(H2_a, c(.05, .50, .95))
      points(h2[2], 2, pch = 16)
      lines(c(h2[1], h2[3]), c(2, 2), lwd = 2)
      
      ## Qst
      qst <- quantile(Qst_population_b, c(.05, .50, .95))
      points(qst[2], 3, pch = 16)
      lines(c(qst[1], qst[3]), c(3, 3), lwd = 2)
      
      labels <- c(expression(paste('H'^2)), 
                  expression(paste("Q"["st"])) )
      axis(side = 2, at = c(2:3), labels = labels, las = 2)
      abline(v = 0, col = "gray90", lty = 2)
      
      dev.off()
      
      #####
      ## Save Credible Intervals 
      #####
      
      credible_intervals <- matrix(nrow = 7, ncol = 3)
      
      credible_intervals[1, 1:3] <- ag[1:3]
      credible_intervals[2, 1:3] <- eg[1:3]
      credible_intervals[3, 1:3] <- ap[1:3]
      credible_intervals[4, 1:3] <- apg[1:3]
      credible_intervals[5, 1:3] <- epg[1:3]
      credible_intervals[6, 1:3] <- h2[1:3]
      credible_intervals[7, 1:3] <- qst[1:3]
      
      rnames <- c("ag", "eg", "ap", "apg", "epg", "h2", "qst")
      cnames <- c("q05", "q50", "q95")
      rownames(credible_intervals) <- rnames
      colnames(credible_intervals) <- cnames
      write.csv(credible_intervals, paste0(dropbox_folder, "PhD_Dissertation/H2_Qst/Credible_Intervals_", tissues[k], ".csv"), row.names = T)
    }
  }
}
output_df <- as.data.frame(output_df)  
output_df <- output_df[complete.cases(output_df$Qst_Code),]
output_df[c( "H2_code", "Variance_h2", "Qst_Code", "Variance_qst")] <- apply(output_df[c( "H2_code", "Variance_h2", "Qst_Code", "Variance_qst")], 1:2, function(x) as.numeric(x))
round_df(output_df, 2)

write.csv(output_df, paste0(od, "H2_QST_mean_estimates_", format(Sys.time(), "%m.%d.%y"),".csv"), row.names = F)


# root_runs <- cbind(H2_a, Qst_population_b)
# write.csv(root_runs, paste0(od, "Roots_H2_Qst_predicted_runs.csv"),row.names = F)
