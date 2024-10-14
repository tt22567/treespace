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
         replacement = "data/criteria_treefile_rooted/",
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
output <- paste0(csv_dir, "res_tree_distance_BMGE_Normolise_20.csv")
#write.csv(res,file = output)

#-------------------------------#
#  check NA and infinite value  #
#-------------------------------#

#take care this step, which might remove unnecessary rows and cols
#check whether there are absent value and infinite value
any(is.na(res))
any(is.infinite(res))
which(is.na(res), arr.ind = TRUE)

# Save the workspace to a file
output <- paste0(Rdata_dir, "12_TreeDistance.RData")
save.image(file = output)
