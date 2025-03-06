#Joseph keating code for distance matrix
#Xiumei Lu code for plotting

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
any(is.na(res))
any(is.infinite(res))
which(is.na(res), arr.ind = TRUE)

#remove absent value
res1 <- res[,-which(is.na(res), arr.ind = T)[,1]]  #remove columns
res2 <- res[-which(is.na(res), arr.ind = T)[,1],]  #remove rows

#save
output <- paste0(csv_dir, "res_tree_distance_BMGE_Normalise_noNA_2041.csv")
write.csv(res2,file = output)


#-------------------------------#
# #       stress analysis     # #
#-------------------------------#

#read the distance matrix
input <- paste0(csv_dir, "res_tree_distance_BMGE_Normalise_noNA_2041.csv")
res <- read.csv( file = input, row.names = 1)

dist.matrix<- res

#calculate nmds based on distance matrix
#Get stress values for NMDS plots with axes 1:10
#check the stress, a measure used to evaluate how match between discrete matrix and distance matrix:
#less than 0.1, good; between 0.1-0.2 acceptable; more than 0.2, not good.
stress.values <- sapply(1:10, function(k) {
  tryCatch({
    vegan::metaMDS(dist.matrix, k = k, trymax = 100)$stress
  }, error = function(e) {
    NA
  })
})

# tiff
tiff("NMDS_stress.tiff", width = 6, height = 5, units = "in", res = 300)
barplot(stress.values, names.arg = 1:10, ylim = c(0, max(stress.values, na.rm = TRUE)), 
        ylab = "Stress", xlab = "Dimensions", col = "grey", main = "Stress vs. NMDS Dimensions")
dev.off()  

#-------------------------------#
# #       NMDS analysis       # #
#-------------------------------#

#read the distance matrix
input <- paste0(csv_dir, "res_tree_distance_BMGE_Normalise_noNA_2041.csv")
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

# save plot
output <- paste0(fig_dir, "NMDS.tiff") 
ggsave(output, plot, width = 5.5, height = 5, dpi = 300)

# Save the workspace to a file
output <- paste0(Rdata_dir, "1_TreeDistance_NMDS.RData")
save.image(file = output)
