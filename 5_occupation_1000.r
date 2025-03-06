#Xiumei Lu code

#find the center 1000 genes
#plot bins of each property

#-------------------# 
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls() )

#---------------#
# LOAD PACKAGES #
#---------------#
library(ape)
library(ggplot2)
library(cowplot)  #plot_grid
library(viridis)

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#

# Load package for automatically finding the path to this R script
library( rstudioapi ) 

# Get the path to current open R script
path_to_file <- rstudioapi::getActiveDocumentContext()$path
full_path_to_file <- base::dirname( path_to_file )

# find and replace characters in your path
wd <- gsub( pattern = "Scripts",
            replacement = "csv",
            x = full_path_to_file )

#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/occupation/",
                 x = full_path_to_file )

#csv directory
csv_dir <- gsub( pattern = "Scripts",
                 replacement = "csv/criteria_1000/",
                 x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )

setwd( wd )

#-------#
# TASKS #
#-------#

#---------------------------------------------------#
#           choose central 1000 trees               #
#             use distance matrix                   #
#---------------------------------------------------#

#read the distance matrix
input <- paste0(wd, "/res_tree_distance_BMGE_Normalise_noNA_2041.csv")
dist_mat <- read.csv( file = input, row.names = 1)

# Calculate the sum of each row
row_sums <- rowSums(dist_mat)

# Sort the row sums in descending order
sorted_row_sums <- sort(row_sums, decreasing = FALSE)

# Extract the top 1000 rows
top_1000_row_sums <- sorted_row_sums[1:1000]

# Print the row names of the top 1000 rows
name_1000 <- names(top_1000_row_sums)
print(name_1000)

# Save the names of the 1000 closest points 
output <- paste0(wd,"/closest1000.csv")
#write(name_1000, file = output)  

#---------------------------------------------------#
#          choose periphery 1000 trees              #
#             use distance matrix                   #
#---------------------------------------------------#

#read the distance matrix
input <- paste0(wd, "/res_tree_distance_BMGE_Normalise_noNA_2041.csv")
dist_mat <- read.csv( file = input, row.names = 1)

# Calculate the sum of each row
row_sums <- rowSums(dist_mat)

# Sort the row sums in descending order
sorted_row_sums_de <- sort(row_sums, decreasing = TRUE)

# Extract the top 1000 rows
top_1000_row_sums_de <- sorted_row_sums[1:1000]

# Print the row names of the top 1000 rows
name_1000_de <- names(top_1000_row_sums)
print(name_1000_de)

# Save the names of the 1000 closest points 
output <- paste0(wd,"/periphery1000.csv")
#write(name_1000, file = output)  

#-------------------------------------------#
# plot closest 1000 genes on the tree space #
#-------------------------------------------#

#load distance matrix
input <- paste0(Rdata_dir,'1_TreeDistance_NMDS.Rdata')
load(input)

#load nmds analysis result
nmds2

centre_1000 <- name_1000

center_indices <- which(rownames(nmds2$points) %in% centre_1000)

# Create a dataframe for the highlighted points
highlighted_df <- data.frame(
  MDS1 = nmds2$points[center_indices, 1],
  MDS2 = nmds2$points[center_indices, 2],
  label = "Center_1000"
)

plot <- ggplot(data = nmds2$points, aes(MDS1, MDS2)) +
  geom_point(size = 0.5, alpha = 1.0) +
  geom_point(data = highlighted_df, aes(MDS1, MDS2, color = label), size = 0.5, show.legend = TRUE) +
  # Manually specify the colors and labels for the legend in a plot
  scale_color_manual(values = "#AD0B08", labels = c("Center_1000")) +
  labs(color = "") +
  theme_gray() +  # default grey theme
  #theme_classic() +  # classic empty theme
  theme(
    legend.position = "top",
    plot.title = element_text(size = 10),
    axis.text = element_text(size = 10),  # Preserve axis tick labels for small plots
    legend.text = element_text(size = 10),  # Adjust the font size of the legend text
    axis.title = element_text(size = 10)   # Set axis title size
  ) +
  xlab("NMDS1") +
  ylab("NMDS2")   # Set axis labels

# Save
output <- paste0(fig_dir, "closest1000.tiff")
ggsave(plot, filename = output, width = 5, height = 5, units = "in", dpi = 300)


#------------------------------------------------#
# how different properties capture central genes #
#------------------------------------------------#

centre_1000 <- name_1000

center_indices <- which(rownames(nmds2$points) %in% centre_1000)

# Read the CSV file containing the gene names to be mapped
properties <- read.csv("properties_final_2041.csv", header = TRUE)

# Columns to be sorted in descending order
prop_positive_desc <- c("NS","AL","PPIS","NVS","LHF","ABSS","TTL","Treeness","TreenessdeRCV","PC2","ER")

# Columns to be sorted in ascending order
prop_negative_asc <- c("PG","RCV","LHU","ER","ALBS","Saturation") 
  
# Create an empty list to store the plots
plot_list <- list()

i <- 1

# Loop through each group of genes (sorted in descending order)
for (column in prop_positive_desc) {
  
  # Sort the dataframe by the current column in descending order
  sorted_properties <- properties[order(-properties[[column]]), ]
  
  # Select the top 1000 genes
  top_1000_genes <- sorted_properties[1:1000, ]
  
  # Extract the "tree" column for the top 1000 genes
  top_1000_trees <- top_1000_genes$rownames.tree
  
  #save
  output <- paste0(csv_dir,column,"_MAX_1000.csv")
  #write.csv(top_1000_trees, file = output)
  
  # Find the positions of the current group of genes in the PCA results
  top1000_indices <- which(rownames(nmds2$points) %in% top_1000_trees)

  # Find the shared genes between the centre 1000 genes and the top 1000 genes for the current column
  shared_trees <- centre_1000[centre_1000 %in% top_1000_trees]
  
  # Find the positions of the shared genes in the PCA results
  sharetree_indices <- which(rownames(nmds2$points) %in% shared_trees)
  
  # Calculate the number of shared genes
  occup <- length(shared_trees)
  
  # Property name for the legend
  file_id <- paste(column, "MAX_1000", sep = "_")
  
  # Define the color and label vectors
  labels <- c("Center_1000", file_id, paste('Shared', occup, sep = "_" ))
  colors <- c("gold", "blue", "#FF0000")  # Colors for the plot
  
  # Create the plot and save it
  current_plot <- ggplot(data = nmds2$points, aes(MDS1, MDS2)) +
    geom_point(color = 'black', size = 1.5, alpha = 1.0) +
    geom_point(data = nmds2$points[center_indices,], aes(MDS1, MDS2, color = labels[1]), size = 1.0, show.legend = TRUE) +
    geom_point(data = nmds2$points[top1000_indices,], aes(MDS1, MDS2, color = labels[2]), size = 1.0, show.legend = TRUE) +
    geom_point(data = nmds2$points[sharetree_indices,], aes(MDS1, MDS2, color = labels[3]), size = 1.0, show.legend = TRUE) +
    scale_color_manual(labels = labels, values = colors) +
    labs(color = " ") +  # Set the title of the legend
    theme(legend.position = "top",
          #panel.background = element_rect(fill = "grey"),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),  # Preserve the axis tick labels for small plots
          legend.text = element_text(size = 18),  # Adjust the font size of the legend text
          axis.title = element_text(size = 18)) + # Set the axis title size
    xlab("NMDS1") +
    ylab("NMDS2")  # Set the axis labels
  
  # Add the plot to the list
  plot_list[[i]] <- current_plot
  i <- i + 1
}

# Loop through each group of genes (sorted in ascending order)
for (column in prop_negative_asc) {
  
  # Sort the dataframe by the current column in ascending order and select the top 1000 rows
  sorted_properties <- properties[order(properties[[column]]), ]
  
  # Select the top 1000 genes
  top_1000_genes <- sorted_properties[1:1000, ]
  
  # Extract the "tree" column for the top 1000 genes
  top_1000_trees <- top_1000_genes$rownames.tree
  
  #save
  output <- paste0(csv_dir,column,"_MIN_1000.csv")
  #write.csv(top_1000_trees, file = output)
  
  # Find the positions of the current group of genes in the PCA results
  gene_indices_group <- which(rownames(nmds2$points) %in% top_1000_trees)
  
  # Find the shared genes between the centre 1000 genes and the top 1000 genes for the current column
  share_trees <- centre_1000[centre_1000 %in% top_1000_trees]
  
  # Find the positions of the shared genes in the PCA results
  sharetree_indices <- which(rownames(nmds2$points) %in% share_trees)
  
  # Calculate the number of shared genes to display in the top right corner
  occup <- length(share_trees)
  
  # Property name for the legend
  file_id <- paste(column, "MIN_1000", sep = "_")
  
  # Define the color and label vectors
  labels <- c("Center_1000", file_id, paste('Shared', occup, sep = "_"))
  colors <- c("gold", "blue", "#FF0000")  # Colors for the plot
  
  # Create the plot and save it
  current_plot <- ggplot(data = nmds2$points, aes(MDS1, MDS2)) +
    geom_point(color = 'black', size = 1.5, alpha = 1.0) +
    geom_point(data = nmds2$points[center_indices,], aes(MDS1, MDS2, color = labels[1]), size = 1.0, show.legend = TRUE) +
    geom_point(data = nmds2$points[top1000_indices,], aes(MDS1, MDS2, color = labels[2]), size = 1.0, show.legend = TRUE) +
    geom_point(data = nmds2$points[sharetree_indices,], aes(MDS1, MDS2, color = labels[3]), size = 1.0, show.legend = TRUE) +
    scale_color_manual(labels = labels, values = colors) +
    labs(color = " ") +
    theme(legend.position = "top",
          #panel.background = element_rect(fill = "grey"),
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),  # Preserve the axis tick labels for small plots
          legend.text = element_text(size = 18),  # Adjust the font size of the legend text
          axis.title = element_text(size = 18)) + # Set the axis title size
    xlab("NMDS1") +
    ylab("NMDS2")  # Set the axis labels
  
  # Add the plot to the list
  plot_list[[i]] <- current_plot
  i <- i + 1
}


#select the 1000 genes with medium evolutionary rate

# Calculate the median value
median_ER <- median(properties$ER)

# Compute the absolute differences from the median
properties$abs_diff_from_median <- abs(properties$ER - median_ER)

# Sort the data frame by the absolute differences
sorted_properties <- properties[order(properties$abs_diff_from_median), ]

# Select the top 1000 values closest to the median
values_median_1000 <- sorted_properties[1:1000, ]

# Extract the name for the selected trees
name_median_1000 <- values_median_1000$rownames.tree

#save
output <- paste0(csv_dir,"ER_IM_1000.csv")
#write.csv(name_median_1000, file = output)

# Find the indices of these trees in the PCA results
indices_median_1000 <- which(rownames(nmds2$points) %in% name_median_1000)

#find the shared trees
share_trees <- centre_1000[centre_1000 %in% name_median_1000]

# Extract indices of shared genes in the PCA analysis results
sharetree_indices <- which(rownames(nmds2$points) %in% share_trees)

# Number of shared genes to be displayed in the top right corner
occup <- length(share_trees)

# Prepare property name for legend
file_id <- paste( "ER_IM_1000", sep = "_" )

# Define color vector and label vector
labels <- c("Center_1000",file_id, paste('Shared', occup, sep = "_" ))
colors <- c("gold", "blue", "#FF0000")   ##E95146

# Construct and save the plot
current_plot <- ggplot(data = nmds2$points, aes(MDS1, MDS2)) +
  geom_point(color = 'black',size = 1.5, alpha = 1.0) +
  geom_point(data = nmds2$points[center_indices,], aes(MDS1, MDS2, color = labels[1]), size = 1.0, show.legend = TRUE) +
  geom_point(data = nmds2$points[indices_median_1000,], aes(MDS1, MDS2, color = labels[2]), size = 1.0, show.legend = TRUE) +
  geom_point(data = nmds2$points[sharetree_indices,], aes(MDS1, MDS2, color = labels[3]), size = 1.0, show.legend = TRUE) +
  #manually specify the colors and labels for the legend in a plot
  scale_color_manual(labels = labels,values = colors) +
  labs(color = " ") +
  theme(legend.position = "top",
        #panel.background = element_rect(fill = "grey"),
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 18),  # Preserve the axis tick labels for small plots
        legend.text = element_text(size = 18),  # Adjust the font size of the legend text
        axis.title = element_text(size = 18)) + # Set the axis title size
  xlab("NMDS1") +
  ylab("NMDS2")  # Set the axis labels

# Add the current plot to the list
plot_list[[i]] <- current_plot
i <- i + 1

# Combine all small plots into one large plot
combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = "h")                                                                

# Save the combined large plot
output <- paste0(fig_dir,"occupation_1000_NMDS.tiff")
ggsave(output, plot = combined_plot, width = 25, height = 38, dpi = 300)

# Save the workspace to a file
output <- paste0(Rdata_dir, "5_occupation_1000_NMDS.RData") 
save.image(file = output)
