#Xiumei Lu and Joseph keating code

#-------------------# 
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls() )

#---------------#
# LOAD PACKAGES #
#---------------#
#install.packages('phangorn')
#install.packages("Quartet")
#install.packages("future.apply")
#install.packages("ggplot2")
#install.packages("factoextra")
#install.packages("phytools")
#install.packages("ggfortify")

library(phangorn)
library(Quartet)
library(future.apply)
library(ggplot2)
library(factoextra)
library(phytools)
library(ggfortify)

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package for automatically finding the path to this R script
library( rstudioapi ) 
# Get the path to current open R script
path_to_file <- rstudioapi::getActiveDocumentContext()$path
full_path_to_file <- base::dirname( path_to_file )

#work directory
wd <- gsub( pattern = "Scripts",
         replacement = "data/MFP_BMGE_TREE_2048",
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

# Load your `Functions.R` script
fun_dir <- full_path_to_file
source( paste( fun_dir, "/Functions.R", sep = "" ))

#-------#
# TASKS #
#-------#

#---------------------------#
# #calculate tree distance ##
#---------------------------#

#list all of the tree filenames
filenames <- list.files( path = wd, pattern = ".treefile" )
length(filenames)

#make an empty vector of mode "list"
#trees <- list() # incorrectly management of computer memory
trees <- vector( mode = "list", length = length( filenames ) )

#read every tree into a single object (called trees)
for(i in 1:length(filenames)){
  trees[[i]] <- phytools::read.newick( file = filenames[i] )
}

# Change to class "multiPhylo" so that subsequent analyses
# can take place
class(trees) <- "multiPhylo"

# calculate tree distances
res <- tree_dist(trees = trees, method = "quartet", normalise = T)  #this is a function

#change the row and coln name of res to whatever you want
# Find whatever matches the regular expression in argument
##replace "pattern" into "replacement"
name_tt <- filenames

# Replace row names
rownames(res) <- name_tt

# Replace column names
colnames(res) <- name_tt

#save
output <- paste0(csv_dir, "res_tree_distance_BMGE_Normolise_2048.csv")
#write.csv(res,file = output)

#-------------------------------#
# #remove NA and infinite value #
#-------------------------------#

#Remove the same species with NA from both columns and rows.
#Be cautious with this step, as it may remove rows and columns unnecessarily
#For instance, if multiple NAs are present for a single species, directly removing all rows or columns containing NAs could lead to data loss.
#Therefore, it might be easier to remove them manually.
#check whether there are absent value and infinite value
any(is.na(res))
any(is.infinite(res))
which(is.na(res), arr.ind = TRUE)

#remove absent value
res1 <- res[,-which(is.na(res), arr.ind = T)[,1]]  #remove columns
res2 <- res[-which(is.na(res), arr.ind = T)[,1],]  #remove rows

#save
output <- paste0(csv_dir, "res_tree_distance_BMGE_Normalise_noNA_2041.csv")
#write.csv(res2,file = output)

#-------------------------------#
# #       PCA analysis        # #
#-------------------------------#

#read the distance matrix
input <- paste0(csv_dir, "res_tree_distance_BMGE_Normalise_noNA_2041.csv")
#res <- read.csv( file = input, row.names = 1)

#PCA analysis for dimension reduction
pca_res <- prcomp(res) 

#get a summary of the PCA results
sum_res <- summary(pca_res)
attributes(sum_res)

# Print rows of 'importance'
rownames(sum_res$importance)
# First row: SD
# Second row: Prop. of variance
# Third row: Cum. prop.
# Print columns
colnames(sum_res$importance)
# Show the information for the first 5 PCs
sum_res$importance[,1:5]
#look at the proportion of variance explained by each of the PC's
sum_res$importance[2,]
#by the first two PC's
sum_res$importance[2,1:2]

#check scores 
pca_res$x[1:5, 1:5]
#save the first two PCs
PCs <- pca_res$x[, 1:2]
#write.csv(PCs,file="pc1_pc2.csv")

#-------------------------------#
# #       PCA visialize       # #
#-------------------------------#

# To plot "outside box", do not send `ggplot` to an object
grDevices::windows()

# To save the output of ggplot in an object, run below
plot <- ggplot(pca_res$x, aes(PC1, PC2)) +
  geom_point(color = "black", size = 1.0, alpha = 1.0) +
  theme_gray() +  # This sets the default grey theme
  xlab("PC1 (73.7%)") +
  ylab("PC2 (3.2%)")

plot

# save plot
output <- paste0(fig_dir, "plot.tiff") 
#ggsave(output, plot, width = 6, height = 4, dpi = 300)

# Save the workspace to a file
output <- paste0(Rdata_dir, "1_TreeDistance.RData")
save.image(file = output)
