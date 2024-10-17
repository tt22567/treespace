#Xiumei Lu code
#Modified from Mongiardino Koch (2021).

#-------------------# 
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls() )

#---------------#
# LOAD PACKAGES #
#---------------#
library(ggplot2)
library(dplyr)
library(ape)

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package for automatically finding the path to this R script
library( rstudioapi ) 
# Get the path to current open R script
path_to_file <- rstudioapi::getActiveDocumentContext()$path
full_path_to_file <- base::dirname( path_to_file )

#csv directory
wd <- gsub( pattern = "Scripts",
                 replacement = "csv",
                 x = full_path_to_file )

# find and replace characters in your path
setwd( wd )
			 
#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/",
                 x = full_path_to_file )

#csv directory
csv_dir <- gsub( pattern = "Scripts",
                 replacement = "csv/",
                 x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )
				   
#-------#
# TASKS #
#-------#

#-----------------------#
# remove periphery loci #
#-----------------------#

#read properties
properties <- read.csv('properties_summary_2041.csv')

#-------------------------------#
# #remove na and infinite value #
#-------------------------------#

#check whether there are absent value and infinite value
any(is.na(properties))
which(is.na(properties), arr.ind = TRUE)
NA_col <- properties[1475,]

#remove absent value, also remove from tree distance matrix
properties2 <- properties[,-which(is.na(properties), arr.ind = T)[,1]]  #remove columns
properties3 <- properties2[-which(is.na(properties), arr.ind = T)[,1],]  #remove rows

#save rownames
rowname <- properties3[ , 1:2]
  
#save
output <- paste0(csv_dir,"properties_noNA_2040.csv")
#write.csv(properties3, file = output)

#-------------------------------#
#         Perform PCA           #
#-------------------------------#

#data for analysis
properties4 <- properties3[ , 3:17]
  
# Perform PCA using prcomp
pca <- princomp(properties4, scale. = TRUE)
str(pca)

# Summary to interpret variance explained
summary(pca)
#95.5, 2.5

# Add PCA scores to the properties data frame
properties4$PC_1 = pca$scores[, 1]
properties4$PC_2 = pca$scores[, 2]

# check Number of dimensions
n_dims <- ncol(pca$scores)
# Mean vector
mean_vector <- rep(0, n_dims)

# Reduce dimensionality if necessary
pca_reduced <- pca$scores[, 1:min(14, ncol(pca$scores))]

# Recompute covariance matrix
cov_matrix <- cov(pca_reduced)
if (det(cov_matrix) == 0) {
  # Handle singular covariance matrix
  cov_matrix <- cov_matrix + diag(n_dims) * 1e-5
}

# Calculate Mahalanobis distances
distances <- mahalanobis(pca_reduced, rep(0, ncol(pca_reduced)), cov_matrix)
outliers <- order(distances, decreasing = TRUE)[1:(nrow(properties4) / 100)]
properties4$outlier <- ifelse(1:nrow(properties4) %in% outliers, 'red', 'black')

# To plot "outside box", do not send `ggplot` to an object
grDevices::windows()

# Plot PCA scores with outliers highlighted
plot <- ggplot(properties4, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.3, aes(color = outlier), shape = 16, size = 3) + theme_bw() +
  scale_color_manual(values = c('#BFD3E6', '#E93658')) + theme(legend.position = "none")+
  theme(aspect.ratio = 1)+
  xlab("PC1 (95.5%)") +
  ylab("PC2 (2.5%)")

plot

# save plot
output <- paste0(fig_dir, "property_pca_outlier.tiff") 
#ggsave(output, plot, width = 6, height = 6, dpi = 300)

#-------------------------------#
#         remove outlier         #
#-------------------------------#

#detect outliers, remove (top 1%) and repeat
outliers = sort(outliers)
re_pca = princomp(properties4[-outliers,c(1:16)], cor = T, scores = T)

scores1 = re_pca$scores[,1]
scores2 = re_pca$scores[,2]

for(i in 1:length(outliers)) {
  if(outliers[i] == 1) {
    scores1 = c(NA, scores1)
    scores2 = c(NA, scores2)
  } else {
    if(outliers[i] < length(scores1)) {
      scores1 = c(scores1[1:(outliers[i]-1)], NA, scores1[outliers[i]:length(scores1)])
      scores2 = c(scores2[1:(outliers[i]-1)], NA, scores2[outliers[i]:length(scores2)])
    } else {
      scores1 = c(scores1, NA)
      scores2 = c(scores2, NA)
    }
  }
}

#add to property matrix the scores of loci along PCs 1 and 2
properties4$PC1 = scores1
properties4$PC2 = scores2

# Assign these row names to properties
properties4$rownames <- rowname

# Add additional NA columns
NA_col <- cbind(NA_col[,3:17], PC_1=NA, PC_2=NA,	outlier=NA,	PC1=NA, PC2=NA,	
                rownames='BSC_100162_mft.bmg.faa,BSC_100162_mft.bmg.faa.treefile')
length(NA_col)
length(colnames(properties4))

# Append the row to the end of properties4
properties5 <- rbind(properties4, NA_col)

#save
output <- paste0(csv_dir, "properties_final_2041.csv") 
#write.csv(properties5, file = output)

# Save the workspace to a file
output <- paste0(Rdata_dir, "1_propterty_space.RData") 
save.image(file = output)

