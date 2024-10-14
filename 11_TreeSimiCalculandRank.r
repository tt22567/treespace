#Joe Keating code function of quartet distance
#Xiumei Lu code the rank

#-------------------# 
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls() )

#---------------#
# LOAD PACKAGES #
#---------------#
library(ggplot2)
library(phytools)
library(Quartet)

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

#csv directory
csv_dir <- gsub( pattern = "Scripts",
                 replacement = "csv/",
                 x = full_path_to_file )

#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/",
                 x = full_path_to_file )

# find and replace characters in your path
setwd( wd )

#-------#
# TASKS #
#-------#

#----------------------------------------------------------------------#
# Calculate normolized quartet distance between center tree and others #                           #
#----------------------------------------------------------------------#

# Get the file paths for the gene tree files
criteria_tree_files <- list.files(pattern = ".treefile")

# Assume reference_tree is the reference tree
reference_tree <- phytools::read.newick(file = criteria_tree_files[4])

# Create an empty list to store each gene tree
trees <- vector(mode = "list", length = length(criteria_tree_files))

# Read each tree file into the list
for (i in 1:length(criteria_tree_files)) {
  trees[[i]] <- phytools::read.newick(file = criteria_tree_files[i])
}

# Change the class to "multiPhylo" for subsequent analyses
class(trees) <- "multiPhylo"

# Create an empty vector to store the results
normalized_quartet_diffs <- numeric(length(criteria_tree_files))

# Loop through the gene tree list
for (i in seq_along(trees)) {
  # Check if the labels of the gene tree match those of the reference tree, retain common labels if not
  common_tips <- intersect(trees[[i]]$tip.label, reference_tree$tip.label)
  
  if(length(common_tips) > 0) {
    trees[[i]] <- keep.tip(trees[[i]], common_tips)
    ref_tree_common <- keep.tip(reference_tree, common_tips)
    
    # Calculate the quartet distance between each gene tree and the reference tree
    Qdist <- QuartetStatus(trees = trees[[i]], cf = ref_tree_common)
    
    # Store the normalized quartet distance in the results vector
    normalized_quartet_diffs[i] <- Qdist[, 4] / Qdist[, 2]
  } else {
    # If there are no common labels, set the normalized quartet difference to NA
    normalized_quartet_diffs[i] <- NA
  }
}

# View the results
print(normalized_quartet_diffs)

# Define the output file path
output_file_path <- file.path(csv_dir, "center_criteria_normalized_quartet_diffs.csv")

# Create a data frame to save as CSV
output <- data.frame(criteria_tree_files = criteria_tree_files, normalized_quartet_diffs = normalized_quartet_diffs)

# Write the data frame to a CSV file
write.csv(output, file = output_file_path, row.names = FALSE)


#------------#
#   Ranks    #
#------------#

# Order the data frame by the quartet distance
output <- output[order(-output$normalized_quartet_diffs), ]

# Remove extension
output$criteria_tree_files <- gsub("\\.treefile", "", output$criteria_tree_files)

# Set the factor levels for the 'criteria_tree_files' column based on the quartet distance order
output$criteria_tree_files <- factor(output$criteria_tree_files, levels = output$criteria_tree_files[order(output$normalized_quartet_diffs)])

# Create the bar plot for quartet distances
plot <- ggplot( data = output, aes(x = criteria_tree_files, y = normalized_quartet_diffs)) +
  geom_bar(stat = "identity", fill = "#82AFF9", width = 0.9, #cornflowerblue #82AFF9
           color = "pink", size = 0.0) +  
  labs(x = "Criteria Tree", y = "Quartet Distance") +
  ggtitle("Quartet Distance Rank") +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "gray", linetype = "dashed"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),  # 不显示次要网格线
        panel.border = element_rect(color = "black", fill = NA, size = 1.5),  # 调整边框的颜色和大小
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.title = element_text(size = 10)) +
  coord_cartesian(ylim = c(0, max(output$normalized_quartet_diffs) * 1))
  

# Save the plot as a TIFF file
output_file = paste0(fig_dir, "QT_center_criteria_RANK.tiff")
ggsave(output_file, plot = plot, width = 8, height = 6, units = "in")

# Reset graphics device
dev.off()
while (!is.null(dev.list())) dev.off()


#-----------------#
#   Ranks-color   #
#-----------------#

# Order the data frame by the quartet distance
output <- output[order(-output$normalized_quartet_diffs), ]

# Remove extension
output$criteria_tree_files <- gsub("\\.treefile", "", output$criteria_tree_files)

# Set the factor levels for the 'criteria_tree_files' column based on the quartet distance order
output$criteria_tree_files <- factor(output$criteria_tree_files, levels = output$criteria_tree_files[order(output$normalized_quartet_diffs)])

# Create a color gradient
output$color <- colorRampPalette(c("lightblue", "blue"))(nrow(output))

# Create the bar plot with each bar having a different color
plot <- ggplot(data = output, aes(x = criteria_tree_files, y = normalized_quartet_diffs, fill = criteria_tree_files)) +
  geom_bar(stat = "identity", width = 0.9, color = "pink", size = 0.0) +  
  scale_fill_manual(values = output$color) + 
  labs(x = "Criteria Tree", y = "Quartet Distance") +
  ggtitle("Quartet Distance Rank") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(0, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    legend.position = "none"  # Remove the legend
  )

# Save the plot as a TIFF file
output_file = paste0(fig_dir, "QT_center_criteria_RANK1.tiff")
ggsave(output_file, plot = plot, width = 8, height = 6, units = "in")

# Reset graphics device
dev.off()
while (!is.null(dev.list())) dev.off()


#------------------------#
#   Ranks-default color  #
#------------------------#

plot <- ggplot(data = output, aes(x = criteria_tree_files, y = normalized_quartet_diffs, fill = criteria_tree_files)) +
  geom_bar(stat = "identity", width = 0.9, color = "pink", size = 0.0) +  
  labs(x = "Criteria", y = "Quartet Distance") +
  ggtitle("Tree Distance Rank") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length = unit(0, "cm"),
    panel.background = element_blank(),
    plot.background = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    legend.position = "none"  # Remove the legend
  )

# Save the plot as a TIFF file
output_file = paste0(fig_dir, "QT_center_criteria_RANK2.tiff")
ggsave(output_file, plot = plot, width = 8, height = 6, units = "in")

# Reset graphics device
dev.off()
while (!is.null(dev.list())) dev.off()
