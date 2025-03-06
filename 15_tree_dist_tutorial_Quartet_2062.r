#Joseph keating code for distance matrix
#Xiumei Lu code for plotting

#-------------------# 
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls() )

#---------------#
# LOAD PACKAGES #
#---------------#

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
         replacement = "data/MFP_BMGE_TREE_2062",
         x = full_path_to_file )

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

# find and replace characters in your path
setwd( wd )

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
length(trees)

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
output <- paste0(csv_dir, "res_tree_distance_BMGE_Normolise_2063.csv")
write.csv(res,file = output)

#-------------------------------#
# #remove NA and infinite value #
#-------------------------------#

res <- read.csv(paste0(csv_dir, "res_tree_distance_BMGE_Normolise_2063.csv"))

#Remove the same species with NA from both columns and rows.
#Be cautious with this step, as it may remove rows and columns unnecessarily
#For instance, if multiple NAs are present for a single species, directly removing all rows or columns containing NAs could lead to data loss.
#Therefore, it might be easier to remove them manually.
any(is.na(res))
any(is.infinite(res))
which(is.na(res), arr.ind = TRUE)

#remove absent value
res1 <- res[,-which(is.na(res), arr.ind = T)[,1]]  #remove columns
res2 <- res[-which(is.na(res), arr.ind = T)[,1],]  #remove rows

#save
output <- paste0(csv_dir, "res_tree_distance_BMGE_Normolise_2063_NONA_2062.csv")
#write.csv(res2,file = output)

#-------------------------------#
# #       NMDS analysis       # #
#-------------------------------#

#read the distance matrix
input <- paste0(csv_dir, "res_tree_distance_BMGE_Normolise_2063_NONA_2062.csv")
res <- read.csv( file = input, row.names = 1)

#distance matrix
dist.matrix<- res

#calculate nmds
nmds2 <- vegan::metaMDS(dist.matrix, k = 2, trymax = 100)
nmds2$points
#pca_res$x → nmds2$points

#check whether the rank is the same
head(row.names(res))
nmds2$points[1:5,]

#get data
data = data.frame(nmds2$points)
data$group = row.names(dist.matrix)

# To plot "outside box", do not send `ggplot` to an object
grDevices::windows()

#plot(data)
plot <- ggplot(data, aes(x = MDS1, y = MDS2)) +
  geom_point(size = 1.0) +
  theme_grey() +
  #theme_light() +  # light background, with line
  #theme_classic() + #no background line
  #geom_text(                
  #aes(label = rownames(data)),
  #vjust = 1.5,
  #size = 2,
  #color = "black"
  #) +
  labs(                     # add stress
    subtitle = paste("Stress =", round(nmds2$stress, 3))
  )

plot

# Save the workspace to a file
output <- paste0(Rdata_dir, "15_TreeDistance_2062_NMDS.RData")
save.image(file = output)
