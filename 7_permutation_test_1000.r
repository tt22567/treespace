#Xiumei Lu code

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
library(future)
library(viridis)

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
            replacement = "csv",
            x = full_path_to_file )

#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/",
                 x = full_path_to_file )

#csv directory
csv_dir <- gsub( pattern = "Scripts",
                 replacement = "csv/criteria_1000/",
                 x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )

# find and replace characters in your path
setwd( wd )

#-------#
# TASKS #
#-------#

#---------------------------------------------------#
#      Randomly select 1000 genes 100 times         #
#    calculate the frequency-centre gene plot       #
#---------------------------------------------------#

# Load necessary libraries
library(dplyr)

# Load the list of tree files, center genes, number of center genes captured by properties
alltree2041 <- read.csv('alltree2041.csv', header = TRUE)
centre_1000 <- read.csv('closest1000.csv', header = FALSE)
centre_pros <- read.csv('num_shared_gene_top1000.csv', header = TRUE)

head(alltree2041)
head(centre_1000)
head(centre_pros)

set.seed(42) # For reproducibility
# Create an empty dataframe to store gene names and their selection counts
gene_selection_data <- data.frame(gene = character(), selection_count = integer(), stringsAsFactors = FALSE)

# Initialize a vector to store the frequency of selecting centre genes
selection_frequencies <- integer(100)

# Perform 100 samplings
for (i in 1:100) {
  # Randomly select 1000 genes from the 'tree' column
  sampled_genes <- sample(alltree2041$tree, 1000)
  
  # Update the selection count for each selected gene
  for (gene in sampled_genes) {
    if (gene %in% gene_selection_data$gene) {
      gene_selection_data$selection_count[gene_selection_data$gene == gene] <- 
        gene_selection_data$selection_count[gene_selection_data$gene == gene] + 1
    } else {
      gene_selection_data <- rbind(gene_selection_data, data.frame(gene = gene, selection_count = 1, stringsAsFactors = FALSE))
    }
  }
  
  # Calculate how many of these sampled genes are in centre_1000
  num_centre_genes <- sum(sampled_genes %in% centre_1000$V1)
  
  # Store the count in the frequencies vector
  selection_frequencies[i] <- num_centre_genes
}

# View the results
head(gene_selection_data)

# Sort the dataframe by selection_count in descending order and select the top 1000
top_1000_genes <- gene_selection_data %>%
  arrange(desc(selection_count)) %>%
  head(1000)

#select the top 1000 genes with most frequency
randomness_MAX_1000 <- top_1000_genes$gene

#save
output <- paste0(csv_dir, "Randomness_MAX_1000.csv")
#write.csv(randomness_MAX_1000, file = output)

# Calculate the frequency distribution
frequency_distribution <- table(selection_frequencies)

# Print the frequency distribution
print(frequency_distribution)

# generate colors for legends
generate_desaturated_colors <- function(n, saturation = 1.0) {
  hues <- seq(0, 360, length.out = n + 1)[-1]
  colors <- hcl(h = hues, c = saturation * 100, l = 60)
  return(colors)
}

# set the color
n <- nrow(centre_pros)  # Number of colors needed
colors <- generate_desaturated_colors(n)

#Plot to show the generated colors
#barplot(rep(1, n), col = colors, border = NA, space = 0, main = "Desaturated Colors")

# generate colors for Histogram
#col <- viridis(n = 10, alpha = 1)

#set the types for legends
tip.type <- c(1:19,1:2)

# Set the output file path for PDF
output <- paste0(fig_dir, "Permutation_test_1000.pdf")

# Open PDF device with wider width
pdf(output, width = 14, height = 8.5)  # Increased width to 14 inches, height remains 8.5 inches (standard A4 size)

# Set graphical parameters
par(#lwd = 3.0,
    cex.lab = 1.5,   # Axis labels (xlab, ylab)
    cex.main = 1.5,  # Main title (main)
    cex.axis = 1.5)  # Axis tick labels

# Create histogram
r <- hist(selection_frequencies, breaks = 30, 
          col = "#82AFF9", # #82AFF9
          #border = "blue",  #pink
          main = "Frequency of Random Selection of Central Genes",
          xlab = "Number of Center Genes Captured by Properties", ylab = "Frequency",
          xlim = c(min(selection_frequencies, centre_pros$central_gene), 
                   max(selection_frequencies, centre_pros$central_gene) + 50))

# Add points
points(centre_pros$central_gene, rep(0, nrow(centre_pros)), col = colors, pch = tip.type, cex = 1.5)

# Add legend
legend("topright", legend = centre_pros$criteria, col = colors, 
       pch = tip.type, cex = 1.5, title = " ", box.lty = 0, inset = c(0.1, 0.02))

# Customize axis lines
axis(side = 1, lwd = 2)  # X axis
axis(side = 2, lwd = 2)  # Y axis

# Close the PDF device
dev.off()

# Save the workspace to a file
output <- paste0(Rdata_dir,"Permutation_test.RData")
#save.image(file = output)

