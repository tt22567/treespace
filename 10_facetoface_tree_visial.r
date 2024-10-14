#Xiumei Lu

#-------------------# 
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls() )

#---------------#
# LOAD PACKAGES #
#---------------#

library(shadowtext)
library(ggnewscale)
library(dplyr)
library(ggtree)

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#

# Load package for automatically finding the path to this R script
library( rstudioapi ) 
# Get the path to current open R script
path_to_file <- rstudioapi::getActiveDocumentContext()$path
full_path_to_file <- base::dirname( path_to_file )

# working directory
wd <- gsub( pattern = "Scripts",
            replacement = "data/criteria_treefile_rooted",
            x = full_path_to_file )

#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/supplementary_Fig7_facetofaceTree",
                 x = full_path_to_file )

# find and replace characters in your path
setwd( wd )

#-------#
# TASKS #
#-------#

#------------------------------#
#      face to face tree       #
#------------------------------#

# Read the Newick tree file
tree1 <- read.tree(file.path(wd, "Closest_1000.treefile"))
tree1_name <- "Center1000"

# Read the tree file names in the working directory
filenames <- list.files(path = wd, pattern = "\\.treefile$")

# Loop over all the tree files
for (i in 1:length(filenames)) {
  tree2 <- read.tree(file = file.path(wd, filenames[i]))
  
  # Extract the tree name from the file name
  tree2_name <- sub("\\.treefile$", "", filenames[i])
  
  # Generate the output file name
  output_file <- file.path(fig_dir, paste(tree1_name, "_vs_", tree2_name, ".pdf", sep = ""))
  
  # Set the font size for species names
  par(cex = 0.05)
  
  # Increase the plotting area margins (bottom, left, top, right)
  par(mar = c(0.2, 0.2, 0.2, 0.8))
  
  # Save the image to a PDF file
  pdf(output_file, width = 80 / 2.54, height = 60 / 2.54)  # Convert cm to inches
  
  # Ensure that the tip labels are not empty and match between the trees
  if (length(tree1$tip.label) > 0 && length(tree2$tip.label) > 0) {
    # Create the association matrix
    common_tips <- intersect(tree1$tip.label, tree2$tip.label)
    assoc <- cbind(common_tips, common_tips)
    
    # Plot the trees with associated tips
    cophyloplot(tree1, tree2, 
                assoc = assoc,
                Font = 1,
                tip.color = "blue",
                length.line = 5,
                space = 300,
                gap = 40,
                lwd = 3,
                col = "orange",
                pos = 1)
  } else {
    print(tree2_name)
  }
  
  # Close the graphics device
  dev.off()
}
