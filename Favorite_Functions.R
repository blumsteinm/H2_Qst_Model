##################################################
## Useful Functions
##################################################

## Running Mean
running_filter <- function(x, y, n = 10, func = "mean"){
  
  ## Order the data based on x so running mean is logical
  df <- data.frame(x = x, y = y)
  df <- df[order(df$x),]
  x <- df$x
  y <- df$y
  
  if(func == "mean"){
    y_out <- c()
    for(i in 1:nrow(df)){
      if(i < n){
        val <- mean(y[1:i], na.rm = T)
      }else{
        val <- mean(y[(i-n):i], na.rm = T)
      }
      y_out <- c(y_out, val)
    }
    output_df <- data.frame(x = x, y = y_out)
    return(output_df)
  }
  if(func == "sd"){
    y_out <- c()
    for(i in 1:nrow(df)){
      if(i < n){
        val <- sd(y[1:i], na.rm = T)
      }else{
        val <- sd(y[(i-n):i], na.rm = T)
      }
      y_out <- c(y_out, val)
    }
    output_df <- data.frame(x = x, y = y_out)
    output_df[1:2,2] <- mean(output_df$y, na.rm = T)
    return(output_df)
  }
  
  
}

## Round DF Function
round_df <- function(x, digits) {
  numeric_columns <- sapply(x, function(x) class(x)) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

## Histograms as Proportions
hist_prop <- function(x, ...){
  h <- hist(x, plot = F)
  vals <- h$counts/sum(h$counts)
  b <- barplot(vals, border = "white", col = "forestgreen", ylim = range(pretty(vals)), ...)
  # axis(side = 1, at = b, h$mids )
  text(b, -0.01, xpd = T,h$mids, srt = 45)
}

## Switch all factors to characters
unfactor_df <- function(x, numeric_or_character = "character", ignore.name = NULL) {
     
    ignore_cols <- which(colnames(x) %in% ignore.name)
    factor_columns <- sapply(x, function(x) class(x)) == 'factor'
    factor_columns[ignore_cols] <- FALSE ## Change columns that are factors to is.factor = F if you want to leave untouched
     
     if(numeric_or_character == "character"){
       x[factor_columns] <-  apply( x[factor_columns], 1:2, function(x) as.character(x) )
     }
     if(numeric_or_character == "numeric"){
       x[factor_columns] <-  apply( x[factor_columns], 1:2, function(x) as.numeric(x) )
     }
     
     return(x)
}

factor_df <- function(x) {
     character_columns <- sapply(x, function(x) class(x)) == 'character'
     x[character_columns] <-  lapply( x[character_columns], as.factor )
     return(x)
}

## Instead of plotting frequency or density, plot the proportion of data in each bin, for data comparison
hist_percentage <- function(data1 = data1, data2 = data2, col1 = col1, border1 = border1, col2 = col2, border2 = border2){
  
  ## Load Data and find densities
  h1 <- hist(data1, plot = F)
  h2 <- hist(data2, plot = F)

  ## Find the distance between breaks to generate the x-length of each bar
  xlen_1 <- mean(diff(h1$breaks))
  xlen_2 <- mean(diff(h2$breaks))
  
  ## If the x-values are on different scales, readjust smaller scale
  if(xlen_1 < xlen_2){
    n_breaks <- length(h1$breaks)
    new_breaks <- seq( h1$breaks[1], h1$breaks[n_breaks], xlen_2)
    n_breaks <- length(new_breaks)
    h1 <- hist(data1, breaks = n_breaks, plot = F)
  }
  if(xlen_2 < xlen_1){
    n_breaks <- length(h2$breaks)
    new_breaks <- seq( h2$breaks[1], h2$breaks[n_breaks], xlen_1)
    n_breaks <- length(new_breaks)
    h2 <- hist(data2, breaks = n_breaks, plot = F)
  }
  
  ## Convert to proportion of data
  h1$density <- h1$counts / sum(h1$counts)
  h2$density <- h2$counts / sum(h2$counts)
  
  ## Recalculate x-bar lengths
  xlen_1 <- mean(diff(h1$breaks))
  xlen_2 <- mean(diff(h2$breaks))
  
  xlen <- max(xlen_1, xlen_2)/2
  
  ## Set up rectangle coordinates
  coord_df_1 <- data.frame(x1 = h1$mids - xlen, 
                           x2 = h1$mids - xlen, 
                           x3 = h1$mids + xlen,
                           x4 = h1$mids + xlen,
                           y1 = 0, y2 = h1$density, y3 = h1$density, y4 = 0)
  coord_df_2 <- data.frame(x1 = h2$mids - xlen, 
                           x2 = h2$mids - xlen, 
                           x3 = h2$mids + xlen,
                           x4 = h2$mids + xlen,
                           y1 = 0, y2 = h2$density, y3 = h2$density, y4 = 0)
  
  ## Draw rectangles
  apply(coord_df_1, 1, function(x) polygon(x[1:4], x[5:8], col = col1, border = border1))
  apply(coord_df_2, 1, function(x) polygon(x[1:4], x[5:8], col = col2, border = border2))
}


## Color Ramp for rasters
clr_ramp <- function(n, plotclr, text, horizontal = F, xl = 280000, xr = 290000, yb =860000 , yt = 920000, type = "categorical"){
  
  if(horizontal == F){
    
    ## Set up Plot Parameters
    par(mar = c(0,0,0,0))
    color_ramp <- colorRampPalette(plotclr)(n)
    
    ## Draw rectangles 
    jump <- (yt - yb) / n
    y_corners <- seq(yb, yt, by = jump)
    if(type == "categorical"){y_midpoints <- y_corners[-(n+1)] + jump/2}
    if(type == "continuous"){y_midpoints <- y_corners}
    for(i in 1:(length(y_corners)-1)){
      rect(xl, y_corners[i], xr, y_corners[i+1], col = color_ramp[i], border = NA)
    }
    
    ## Add Text
    if(length(text) == 2){y_midpoints <- c(y_midpoints[1], y_midpoints[n])}
    text(xr, y_midpoints, text, family = "serif",adj = c(0,0), cex = 1)
  }
  par(mar = c(4.5, 4.5, 0, 0))
}

sunCalc<-function(d,lat,long){
  if(any(is.null(c(d,lat,long)))){
    stop("day, latitude and Longitude must be supplied")
  }
  
  date<-as.Date(d,format="%Y-%m-%d")
  
  decimal.day<-as.numeric(format(date,format="%j"))
  
  ## Function to convert degrees to radians
  rad<-function(x)pi*x/180
  
  ##Radius of the earth (km)
  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)
  
  ##Convert observer's latitude to radians
  L=rad(lat)
  
  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  timezone = -4*(abs(long)%%15)*sign(long)
  
  ## The earth's mean distance from the sun (km)
  r = 149598000
  
  theta = 2*pi/365.25*(decimal.day-80)
  
  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)
  
  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun
  that = t0+5 
  
  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(decimal.day-80)/365.25)+8*sin(2*pi*decimal.day/365.25)
  
  ## now sunrise and sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  
  srH<-floor(sunrise)
  srM<-floor((sunrise-srH)*60)
  
  ssH<-floor(sunset)
  ssM<-floor((sunset-ssH)*60)
  
  SR<-paste(date,paste(srH, srM,sep=":"),sep=" ")
  SS<-paste(date,paste(ssH, ssM,sep=":"),sep=" ")
  
  return(list("sunrise" = SR,"sunset" = SS, "daylength" = sunset - sunrise))
}

##################################################
## Color Pallettes 
##################################################
choose_colors <- function(alpha){
  
                default_par <- par()
                alpha <- alpha
                
                color_list <- list(maxres(alpha), cinmint(alpha), instaCl(alpha), skittles(alpha), landuse(alpha), porcelin(alpha), castle(alpha),
                                   beachGlass(alpha), gainingHeat(alpha), fishtank(alpha), sunsetCamping(alpha), comics(alpha), 
                                   dark_ocean(alpha), ocean(alpha), long(alpha))
                color_nmes <- c("maxres", "cinmint", "instaCl", "skittles", "landuse", "porcelin", "castle", "beachGlass", "gainingHeat", 
                                "fishtank", "sunsetCamping", "comics", "dark_ocean", "ocean", "long")
                par(mfrow = c(4, 3), mar = c(0, 0.2, 2, 0.2))
                for(i in 1:length(color_list)){
                  col_subset <- color_list[[i]]
                  barplot(rep(1, length(col_subset)), col = col_subset, border = "white", axes = F, xlim = c(0, (length(col_subset) + 1.5)), ylim = c(-0.5, 1.5))
                  title(main = color_nmes[i])
                  axis(side = 4, lwd = 5, lwd.ticks = 0, labels = F)
                }
                par <- default_par
}

maxres <- function(alpha){
            alpha <- alpha 
            col_pal <- c(rgb(2, 135, 169, max = 255, alpha = alpha),
            rgb(127, 211, 194, max = 255, alpha = alpha),
            rgb(243, 241, 209, max = 255, alpha = alpha),
            rgb(254, 147, 91, max = 255, alpha = alpha),
            rgb(240, 104, 75, max = 255, alpha = alpha))
            return(col_pal)
}
cinmint <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(20, 166, 151, max = 255, alpha = alpha),
               rgb(242, 193, 46, max = 255, alpha = alpha),
               rgb(242, 157, 53, max = 255, alpha = alpha),
               rgb(242, 118, 73, max = 255, alpha = alpha),
               rgb(242, 82, 82, max = 255, alpha = alpha))
  return(col_pal)
}
instaCl <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(37, 39, 42, max = 255, alpha = alpha),
               rgb(23, 90, 139, max = 255, alpha = alpha),
               rgb(102, 147, 178, max = 255, alpha = alpha),
               rgb(198, 200, 201, max = 255, alpha = alpha),
               rgb(240, 240, 241, max = 255, alpha = alpha))
  return(col_pal)
}
skittles <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(207, 66, 50, max = 255, alpha = alpha),
               rgb(235, 127, 35, max = 255, alpha = alpha),
               rgb(250, 192, 35, max = 255, alpha = alpha),
               rgb(6, 134, 117, max = 255, alpha = alpha),
               rgb(107, 32, 46, max = 255, alpha = alpha))
  return(col_pal)
}
landuse <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(162, 0, 1, max = 255, alpha = alpha),
               rgb(198, 76, 37, max = 255, alpha = alpha),
               rgb(249, 239, 74, max = 255, alpha = alpha),
               rgb(75, 119, 68, max = 255, alpha = alpha),
               rgb(29, 54, 74, max = 255, alpha = alpha))
  return(col_pal)
}
porcelin <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(15, 31, 56, max = 255, alpha = alpha), 
               rgb(253, 60, 60, max = 255, alpha = alpha),
               rgb(255, 183, 76, max = 255, alpha = alpha),
               rgb(19, 141, 144, max = 255, alpha = alpha), 
               rgb(253, 246, 246, max = 255, alpha = alpha))
  return(col_pal)
}
castle <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(123, 42, 59, max = 255, alpha = alpha), 
               rgb(229, 118, 97, max = 255, alpha = alpha),
               rgb(248, 197, 140, max = 255, alpha = alpha),
               rgb(248, 231, 162, max = 255, alpha = alpha), 
               rgb(134, 221, 178, max = 255, alpha = alpha))
  return(col_pal)
}
beachGlass <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(255, 246, 201, max = 255, alpha = alpha), 
               rgb(200, 232, 199, max = 255, alpha = alpha),
               rgb(164, 222, 171, max = 255, alpha = alpha),
               rgb(133, 204, 159, max = 255, alpha = alpha), 
               rgb(73, 158, 141, max = 255, alpha = alpha))
  return(col_pal)
}
gainingHeat <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(255, 230, 230, max = 255, alpha = alpha), 
               rgb(255, 172, 172, max = 255, alpha = alpha),
               rgb(255, 115, 115, max = 255, alpha = alpha),
               rgb(255, 58, 58, max = 255, alpha = alpha), 
               rgb(255, 0, 0, max = 255, alpha = alpha))
  return(col_pal)
}
fishtank <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(48, 119, 156, max = 255, alpha = alpha), 
               rgb(97, 232, 215, max = 255, alpha = alpha),
               rgb(104, 170, 71, max = 255, alpha = alpha),
               rgb(255, 220, 89, max = 255, alpha = alpha), 
               rgb(232, 98, 67, max = 255, alpha = alpha))
  return(col_pal)
}
sunsetCamping <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(46, 27, 45, max = 255, alpha = alpha), 
               rgb(84, 0, 50, max = 255, alpha = alpha),
               rgb(130, 3, 51, max = 255, alpha = alpha),
               rgb(201, 40, 62, max = 255, alpha = alpha), 
               rgb(240, 67, 58, max = 255, alpha = alpha))
  return(col_pal)
}
comics <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(191, 48, 86, max = 255, alpha = alpha), 
               rgb(79, 73, 115, max = 255, alpha = alpha),
               rgb(67, 189, 217, max = 255, alpha = alpha),
               rgb(95, 191, 80, max = 255, alpha = alpha), 
               rgb(217, 208, 89, max = 255, alpha = alpha))
  return(col_pal)
}
dark_ocean <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(16, 19, 38, max = 255, alpha = alpha), 
               rgb(43, 53, 65, max = 255, alpha = alpha),
               rgb(84, 102, 114, max = 255, alpha = alpha),
               rgb(148, 167, 163, max = 255, alpha = alpha), 
               rgb(182, 191, 186, max = 255, alpha = alpha))
  return(col_pal)
}
ocean <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(0, 56, 64, max = 255, alpha = alpha), 
               rgb(0, 90, 91, max = 255, alpha = alpha),
               rgb(0, 115, 105, max = 255, alpha = alpha),
               rgb(0, 140, 114, max = 255, alpha = alpha), 
               rgb(2, 166, 118, max = 255, alpha = alpha))
  return(col_pal)
}
long <- function(alpha){
  alpha <- alpha 
  col_pal <- c(rgb(28, 51, 81, max = 255, alpha = alpha), 
               rgb(28, 76, 81, max = 255, alpha = alpha),
               rgb(28, 100, 81, max = 255, alpha = alpha),
               rgb(28, 124, 81, max = 255, alpha = alpha), 
               rgb(148, 160, 81, max = 255, alpha = alpha),
               
               rgb(245, 201, 117, max = 255, alpha = alpha), 
               rgb(245, 168, 96, max = 255, alpha = alpha),
               rgb(245, 136, 88, max = 255, alpha = alpha),
               rgb(245, 94, 74, max = 255, alpha = alpha), 
               rgb(245, 52, 58, max = 255, alpha = alpha),
               
               rgb(180, 46, 65, max = 255, alpha = alpha), 
               rgb(143, 29, 44, max = 255, alpha = alpha),
               rgb(90, 20, 42, max = 255, alpha = alpha),
               rgb(64, 13, 42, max = 255, alpha = alpha), 
               rgb(20, 10, 37, max = 255, alpha = alpha)
               )
  return(col_pal)
}

##################################################
## SNP Subsetting
##################################################

subset_by_SNP <- function(GWAS_file, Nsnps, output_file){
  
  ## Read and oranize GWAS results
  snps <- read.table(GWAS_file)
  colnames(snps) <- c("Chr", "Pos", "P-Val")
  snps <- snps[order(snps$`P-Val`),]
  
  ## Just the number of snps you are iterating over 
  snp_rows <- snps[1:Nsnps,] 
  snp_rows <- snp_rows[order(snp_rows$Chr, as.numeric(snp_rows$Pos)),] ## order by chromosome, then position for efficiency
  # write.csv(snp_rows, snp_info, row.names = F)
  
  ## Save out the values for selecting rows to textfile to pass to bash
  cat("unique_chr=(",paste(unique(snp_rows$Chr), collapse = " "), ")" , file = output_file, sep = "\n")
  cat("snp_chr=(",paste(snp_rows$Chr, collapse = " "), ")"            , file = output_file, sep = "\n", append = T)
  cat("snp_order=(",paste(snp_rows$Order, collapse = " "), ")"        , file = output_file, sep = "\n", append = T)
  cat("snp_row=(",paste(snp_rows$Pos, collapse = " "), ")"            , file = output_file, sep = "\n", append = T)
  
  # system("/Users/Meghs/Dropbox/PhD_Dissertation/Code/Genotype_Predictions/Select_Rows_From_SNPs.txt")
}

pseudo_manhattan_base <- function (tped, correlations, values = "p", ...) {
  require(fields)
  if (is.character(tped)) {
    message("loading tped file")
    tped <- fread(tped, data.table = F)
  }
  snps_analyzed <- which(rownames(correlations$Pvalues) %in% 
                           tped[, 2])
  map <- tped[snps_analyzed, 1:2]
  if (values == "p") {
    likelyhood_sum <- colSums(-log(correlations$Pvalues))
    data_p <- as.data.frame(cbind(likelyhood_sum, c(1:length(likelyhood_sum)),map[, 1]))
    colnames(data_p) <- c("like_sum", "Nvar", "chr")
    chr.color <- color.scale(data_p$chr, col = rainbow(19, alpha = .7))
    par(mar = c(3,3, 0.5, 0.5), mgp = c(1.2,0.5, 0))
    plot(like_sum ~ Nvar, data = data_p, col = chr.color, pch = 16, family = "serif", xlab = "Index", ylab = "Sum of -Log Likelihood", ...)
    
  }
  if (values == "c") {
    
    likelyhood_sum <- colSums(abs(correlations$Coefficients))
    data_p <- as.data.frame(cbind(likelyhood_sum, c(1:length(likelyhood_sum)), map[, 1]))
    colnames(data_p) <- c("like_sum", "Nvar", "chr")
    chr.color <- color.scale(data_p$chr, col = rainbow(19, alpha = .7))
    par(mar = c(3,3, 0.5, 0.5), mgp = c(1.2,0.5, 0))
    plot(like_sum ~ Nvar, data = data_p, col = chr.color, pch = 16, family = "serif", xlab = "Index", ylab = "Sum of Effect Size", ...)
    
  }
}


##################################################
## Plots
##################################################


## Customized Biplot for PCAs
mb_biplot <- function(pca, arrow_col, pnt_col, pcs_to_plot){
  
  par_settings <- par()
  par(mgp = c(1.5, 0.5, 0), mar = c(2.5, 3.5, 3.5, 2.5))
  
  ## Rescale data
  pcs_to_plot = 1:2                              # Selecting first two PC's
  scale = 1                                      # Default
  scores= pca$x                                  # The scores
  lam = pca$sdev[pcs_to_plot]                    # Sqrt e-vals (lambda) 2 PC's
  n = nrow(scores)                               # no. rows scores
  lam = lam * sqrt(n)                            # See below.
  x = t(t(scores[,pcs_to_plot])       / lam)     # scaled scores
  y = t(t(pca$rotation[,pcs_to_plot]) * lam)     # scaled eigenvecs (loadings)
  n = nrow(x)                                    # Same as dataset 
  p = nrow(y)                                    # Num. of variables
  imp <- summary(pca)$importance
  print(paste("Lam scaling variable is: ", lam))
  
  # Function to get the range:
  unsigned.range = function(x) c(-abs(min(x, na.rm = TRUE)), abs(max(x, na.rm = TRUE)))
  rangx1 = unsigned.range(x[, 1L])               # Range first col x
  rangx2 = unsigned.range(x[, 2L])               # Range second col x
  rangy1 = unsigned.range(y[, 1L])               # Range 1st scaled evec
  rangy2 = unsigned.range(y[, 2L])               # Range 2nd scaled evec
  xlim = ylim = rangx1 = rangx2 = range(rangx1, rangx2)
  
  # And the critical value is the maximum of the ratios of ranges of 
  # scaled e-vectors / scaled scores:
  ratio = max(rangy1/rangx1, rangy2/rangx2)
  
  # Plotting Genotypes
  par(pty = "s")                                 
  plot(x, type = "n", xlim = xlim, ylim = ylim, xlab = paste0("PC ", pcs_to_plot[1], ": ", round(imp[2,pcs_to_plot[1]],2)), ylab = paste0("PC ", pcs_to_plot[2], ": ", round(imp[2,pcs_to_plot[2]],2)))  # No points
  points(x[, 1:2], pch = 16, col = pnt_col)
  
  ## Plotting climate lines
  par(new = TRUE)                            
  
  # Setting x and y limits for the arrows:
  xlim = xlim * ratio  # We multiply the original limits x ratio
  ylim = ylim * ratio  # ... for both the x and y axis
  
  ## make blank plot on top of og plot
  plot(y, axes = FALSE, type = "n", xlim = xlim, ylim = ylim, xlab = "", ylab = "")
  axis(3); axis(4)
  # text(y, labels = ylabs, col = 2)  # This just prints the species
  arrows(0, 0, y[, 1L] * 0.8, y[, 2L] * 0.8, length = 0.1, col = arrow_col, lwd = 2)
  
  for(j in 1:nrow(pca$rotation)){
    
    Label_Name <- rownames(pca$rotation)[j]
    hp <- c(y[j,1],y[j,2]) * 0.8
    ## Label
    if(hp[1] >= 0){
      adj_var <- c(0,0)
    }
    if(hp[1] < 0){
      adj_var <- c(1,0)
    }
    
    if(Label_Name == "DD0"){hp[2] <- hp[2] - 0.02}
    text(hp[1], hp[2], Label_Name, adj = adj_var, srt = 0)
  }
  
  # par(par_settings)
  
}





















