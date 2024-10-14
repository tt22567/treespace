#Xiumei Lu code

#find the closest 500 genes
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

setwd( wd )

#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/occupation/",
                 x = full_path_to_file )

#csv directory
#csv_dir <- gsub( pattern = "Scripts",
                 #replacement = "csv/criteria_500/",
                 #x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )

#-------#
# TASKS #
#-------#

#---------------------------------------------------#
#    choose central 500 trees in multidimension    #
#             use distance matrix                   #
#---------------------------------------------------#

#read the distance matrix
input <- paste0(wd,"/res_tree_distance_BMGE_Normalise_noNA_2041.csv")
dist_mat <- read.csv( file = input, row.names = 1)

# Calculate the sum of each row
row_sums <- rowSums(dist_mat)

# Sort the row sums in descending order
sorted_row_sums <- sort(row_sums, decreasing = FALSE)

# Extract the top 500 rows
top_500_row_sums <- sorted_row_sums[1:500]

# Print the row names of the top 500 rows
name_500 <- names(top_500_row_sums)
print(name_500)

# Save the names of the 500 closest points 
output <- paste0(wd, "/closest500.csv")
write(name_500, file = output)  

#-------------------------------------------#
# plot closest 500 genes on the tree space #
#-------------------------------------------#

#load distance matrix
input <- paste0(Rdata_dir,'1_TreeDistance.Rdata')
load(input)

centre_500 <- name_500

center_indices <- which(rownames(pca_res$x) %in% centre_500)

# Create a dataframe for the highlighted points
highlighted_df <- data.frame(
  PC1 = pca_res$x[center_indices, 1],
  PC2 = pca_res$x[center_indices, 2],
  label = "closest_500"
)

plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
  geom_point(size = 0.5, alpha = 1.0) +
  geom_point(data = highlighted_df, aes(PC1, PC2, color = label), size = 0.5, show.legend = TRUE) +
  # Manually specify the colors and labels for the legend in a plot
  scale_color_manual(values = "red", labels = c("closest_500")) +
  labs(color = "") +
  theme_gray() +  # default grey theme
  #theme_classic() +  # classic empty theme
  theme(
    legend.position = "top",
    plot.title = element_text(size = 6),
    axis.text = element_text(size = 6),  # Preserve axis tick labels for small plots
    legend.text = element_text(size = 6),  # Adjust the font size of the legend text
    axis.title = element_text(size = 6)   # Set axis title size
  ) +
  xlab("PC1 (73.7%)") +
  ylab("PC2 (3.2%)")   # Set axis labels

# Save
output <- paste0(fig_dir, "closest500.tiff")
#ggsave(plot, filename = output, width = 6, height = 3, units = "in", dpi = 300)

#-------------------------------------------------#
# how different properties capture central genes  #
#-------------------------------------------------#

Center_500 <- name_500

center_indices <- which(rownames(pca_res$x) %in% Center_500)

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
  
  # Select the top 500 genes
  top_500_genes <- sorted_properties[1:500, ]
  
  # Extract the "tree" column for the top 500 genes
  top_500_trees <- top_500_genes$rownames.tree
  
  #save
  #output <- paste0(csv_dir,column,"_MAX_500.csv")
  #write.csv(top_500_trees, file = output)
  
  # Find the positions of the current group of genes in the PCA results
  top500_indices <- which(rownames(pca_res$x) %in% top_500_trees)

  # Find the shared genes between the Center 500 genes and the top 500 genes for the current column
  shared_trees <- Center_500[Center_500 %in% top_500_trees]
  
  # Find the positions of the shared genes in the PCA results
  sharetree_indices <- which(rownames(pca_res$x) %in% shared_trees)
  
  # Calculate the number of shared genes
  occup <- length(shared_trees)
  
  # Property name for the legend
  file_id <- paste(column, "MAX_500", sep = "_")
  
  # Define the color and label vectors
  labels <- c("Center_500", file_id, paste('Shared', occup, sep = "_" ))
  colors <- c("gold", "blue", "red")  # Colors for the plot
  
  # Create the plot and save it
  current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
    geom_point(color = 'white', size = 1.0, alpha = 1.0) +
    geom_point(data = pca_res$x[center_indices,], aes(PC1, PC2, color = labels[1]), size = 1.0, show.legend = TRUE) +
    geom_point(data = pca_res$x[top500_indices,], aes(PC1, PC2, color = labels[2]), size = 1.0, show.legend = TRUE) +
    geom_point(data = pca_res$x[sharetree_indices,], aes(PC1, PC2, color = labels[3]), size = 1.0, show.legend = TRUE) +
    scale_color_manual(labels = labels, values = colors) +
    labs(color = " ") +  # Set the title of the legend
    theme(legend.position = "top",
          panel.background = element_rect(fill = "black"),
          plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),  # Preserve the axis tick labels for small plots
          legend.text = element_text(size = 18),  # Adjust the font size of the legend text
          axis.title = element_text(size = 18),  # Set the axis title size
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("PC1 (73.7%)") +
    ylab("PC2 (3.2%)")  # Set the axis labels
  
  # Add the plot to the list
  plot_list[[i]] <- current_plot
  i <- i + 1
}

# Loop through each group of genes (sorted in ascending order)
for (column in prop_negative_asc) {
  
  # Sort the dataframe by the current column in ascending order and select the top 500 rows
  sorted_properties <- properties[order(properties[[column]]), ]
  
  # Select the top 500 genes
  top_500_genes <- sorted_properties[1:500, ]
  
  # Extract the "tree" column for the top 500 genes
  top_500_trees <- top_500_genes$rownames.tree
  
  #save
  #output <- paste0(csv_dir,column,"_MIN_500.csv")
  #write.csv(top_500_trees, file = output)
  
  # Find the positions of the current group of genes in the PCA results
  gene_indices_group <- which(rownames(pca_res$x) %in% top_500_trees)
  
  # Find the shared genes between the Center 500 genes and the top 500 genes for the current column
  share_trees <- Center_500[Center_500 %in% top_500_trees]
  
  # Find the positions of the shared genes in the PCA results
  sharetree_indices <- which(rownames(pca_res$x) %in% share_trees)
  
  # Calculate the number of shared genes to display in the top right corner
  occup <- length(share_trees)
  
  # Property name for the legend
  file_id <- paste(column, "MIN_500", sep = "_")
  
  # Define the color and label vectors
  labels <- c("Center_500", file_id, paste('Shared', occup, sep = "_"))
  colors <- c("gold", "blue", "red")  # Colors for the plot
  
  # Create the plot and save it
  current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
    geom_point(color = 'grey', size = 1.0, alpha = 1.0) +
    geom_point(data = pca_res$x[center_indices,], aes(PC1, PC2, color = labels[1]), size = 1.0, show.legend = TRUE) +
    geom_point(data = pca_res$x[top500_indices,], aes(PC1, PC2, color = labels[2]), size = 1.0, show.legend = TRUE) +
    geom_point(data = pca_res$x[sharetree_indices,], aes(PC1, PC2, color = labels[3]), size = 1.0, show.legend = TRUE) +
    scale_color_manual(labels = labels, values = colors) +
    labs(color = " ") +
    theme(legend.position = "top",
          panel.background = element_rect(fill = "black"),
          plot.title = element_text(size = 18),
          axis.text = element_text(size = 18),  # Preserve the axis tick labels for small plots
          legend.text = element_text(size = 18),  # Adjust the font size of the legend text
          axis.title = element_text(size = 18),  # Set the axis title size
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("PC1 (73.7%)") +
    ylab("PC2 (3.2%)")  # Set the axis labels
  
  # Add the plot to the list
  plot_list[[i]] <- current_plot
  i <- i + 1
}


#select the 500 genes with medium evolutionary rate

# Calculate the median value
median_ER <- median(properties$ER)

# Compute the absolute differences from the median
properties$abs_diff_from_median <- abs(properties$ER - median_ER)

# Sort the data frame by the absolute differences
sorted_properties <- properties[order(properties$abs_diff_from_median), ]

# Select the top 500 values closest to the median
values_median_500 <- sorted_properties[1:500, ]

# Extract the name for the selected trees
name_median_500 <- values_median_500$rownames.tree

#save
#output <- paste0(csv_dir,"ER_IM_500.csv")
#write.csv(name_median_500, file = output)

# Find the indices of these trees in the PCA results
indices_median_500 <- which(rownames(pca_res$x) %in% name_median_500)

#find the shared trees
share_trees <- Center_500[Center_500 %in% name_median_500]

# Extract indices of shared genes in the PCA analysis results
sharetree_indices <- which(rownames(pca_res$x) %in% share_trees)

# Number of shared genes to be displayed in the top right corner
occup <- length(share_trees)

# Prepare property name for legend
file_id <- paste( "ER_IM_500", sep = "_" )

# Define color vector and label vector
labels <- c("Center_500",file_id, paste('Shared', occup, sep = "_" ))
colors <- c("gold", "blue", "red")   ##E95146

# Construct and save the plot
current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
  geom_point(color = 'grey',size = 1.0, alpha = 1.0) +
  geom_point(data = pca_res$x[center_indices,], aes(PC1, PC2, color = labels[1]), size = 1.0, show.legend = TRUE) +
  geom_point(data = pca_res$x[indices_median_500,], aes(PC1, PC2, color = labels[2]), size = 1.0, show.legend = TRUE) +
  geom_point(data = pca_res$x[sharetree_indices,], aes(PC1, PC2, color = labels[3]), size = 1.0, show.legend = TRUE) +
  #manually specify the colors and labels for the legend in a plot
  scale_color_manual(labels = labels,values = colors) +
  labs(color = " ") +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "black"),
        plot.title = element_text(size = 18),
        axis.text = element_text(size = 18),  # Preserve axis tick labels for small plots
        legend.text = element_text(size = 18),  # Adjust the font size of legend text
        axis.title = element_text(size = 18),   # Set axis title size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("PC1 (73.7%)") +
  ylab("PC2 (3.2%)")  # Set the axis labels

# Add the current plot to the list
plot_list[[i]] <- current_plot
i <- i + 1

# Combine all small plots into one large plot
combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = "h")                                                                

# Save the combined large plot
output <- paste0(fig_dir,"occupation_500.tiff")
ggsave(plot = combined_plot, output, width = 25, height = 23, dpi = 300)

# Save the workspace to a file
output <- paste0(Rdata_dir, "5_occupation_500.RData") 
save.image(file = output)






#different color
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

setwd( wd )

#fig directory
fig_dir <- gsub( pattern = "Scripts",
                 replacement = "Figs/occupation/",
                 x = full_path_to_file )

#csv directory
csv_dir <- gsub( pattern = "Scripts",
                 replacement = "csv/criteria_500/",
                 x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )

#-------#
# TASKS #
#-------#

#---------------------------------------------------#
#    choose central 500 trees in multidimension    #
#---------------------------------------------------#

#load distance matrix
input <- paste0(Rdata_dir,'2_TreeDistance.Rdata')
load(input)

#distance matrix
res

#PCA analysis result for dimension reduction
pca_res

#coordination
pca_res$x
head(pca_res$x[1:2, 1:2])

# Calculate the distance of each point to others in multidimension
dist_mat <- stats::dist(pca_res$x, method = "euclidean")  
head(as.matrix(dist_mat)[1:2, 1:2])  

# Calculate the sum of the closest 500 distances for each point, excluding the point itself
D <- apply(as.matrix(dist_mat), 1, function(x) sum(sort(x)[2:501]))
head(D)

# Save the sum of distances to a CSV file
output <- paste0(csv_dir, "dist_mat_closest500.csv")
write.csv(D, file=output)  

# Find the index of the point with the smallest sum of the closest 500 distances
index_pca <- which.min(D)  

# Find the indices of the 500 points closest to this point
index_500 <- order(as.matrix(dist_mat)[index_pca,])[2:501]

# Get the names of its 500 closest points
name_500 <- rownames(pca_res$x)[index_500]  

# Save the names of the 500 closest points 
output <- paste0(csv_dir, "Center500_1.csv")
write(name_500, file = output)  

#-------------------------------------------------#
# how different properties capture central genes  #
#-------------------------------------------------#

Center_500 <- name_500

center_indices <- which(rownames(pca_res$x) %in% Center_500)

# Read the CSV file containing the gene names to be mapped
properties <- read.csv("properties_final_2040.csv", header = TRUE)

# Columns to be sorted in descending order
prop_positive_desc <- c("NS","AL","PPIS","NVS","LHF","ABSS","TTL","Treeness","TreenessdeRCV","PC_two","ER")

# Columns to be sorted in ascending order
prop_negative_asc <- c("PG","RCV","LHU","ER","ALBS","Saturation") 

# Create an empty list to store the plots
plot_list <- list()

i <- 1

# Loop through each group of genes (sorted in descending order)
for (column in prop_positive_desc) {
  
  # Sort the dataframe by the current column in descending order
  sorted_properties <- properties[order(-properties[[column]]), ]
  
  # Select the top 500 genes
  top_500_genes <- sorted_properties[1:500, ]
  
  # Extract the "tree" column for the top 500 genes
  top_500_trees <- top_500_genes$rownames.tree
  
  #save
  #output <- paste0(csv_dir,column,"_MAX_500.csv")
  #write.csv(top_500_trees, file = output)
  
  # Find the positions of the current group of genes in the PCA results
  top500_indices <- which(rownames(pca_res$x) %in% top_500_trees)
  
  # Find the shared genes between the Center 500 genes and the top 500 genes for the current column
  shared_trees <- Center_500[Center_500 %in% top_500_trees]
  
  # Find the positions of the shared genes in the PCA results
  sharetree_indices <- which(rownames(pca_res$x) %in% shared_trees)
  
  # Calculate the number of shared genes
  occup <- length(shared_trees)
  
  # Property name for the legend
  file_id <- paste(column, "MAX_500", sep = "_")
  
  # Define the color and label vectors
  labels <- c("Center_500", file_id, paste('Shared', occup, sep = "_" ))
  colors <- c("gold", "blue", "red")  # Colors for the plot
  
  # Create the plot and save it
  current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
    geom_point(color = 'white', size = 0.3, alpha = 1.0) +
    geom_point(data = pca_res$x[center_indices,], aes(PC1, PC2, color = labels[1]), size = 0.3, show.legend = TRUE) +
    geom_point(data = pca_res$x[top500_indices,], aes(PC1, PC2, color = labels[2]), size = 0.3, show.legend = TRUE) +
    geom_point(data = pca_res$x[sharetree_indices,], aes(PC1, PC2, color = labels[3]), size = 0.3, show.legend = TRUE) +
    scale_color_manual(labels = labels, values = colors) +
    labs(color = " ") +  # Set the title of the legend
    theme(legend.position = "top",
          panel.background = element_rect(fill = "black"),
          plot.title = element_text(size = 10),
          axis.text = element_text(size = 10),  # Preserve the axis tick labels for small plots
          legend.text = element_text(size = 10),  # Adjust the font size of the legend text
          axis.title = element_text(size = 10),  # Set the axis title size
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("PC1 (73.7%)") +
    ylab("PC2 (3.2%)")  # Set the axis labels
  
  # Add the plot to the list
  plot_list[[i]] <- current_plot
  i <- i + 1
}

# Loop through each group of genes (sorted in ascending order)
for (column in prop_negative_asc) {
  
  # Sort the dataframe by the current column in ascending order and select the top 500 rows
  sorted_properties <- properties[order(properties[[column]]), ]
  
  # Select the top 500 genes
  top_500_genes <- sorted_properties[1:500, ]
  
  # Extract the "tree" column for the top 500 genes
  top_500_trees <- top_500_genes$rownames.tree
  
  #save
  #output <- paste0(csv_dir,column,"_MIN_500.csv")
  #write.csv(top_500_trees, file = output)
  
  # Find the positions of the current group of genes in the PCA results
  gene_indices_group <- which(rownames(pca_res$x) %in% top_500_trees)
  
  # Find the shared genes between the Center 500 genes and the top 500 genes for the current column
  share_trees <- Center_500[Center_500 %in% top_500_trees]
  
  # Find the positions of the shared genes in the PCA results
  sharetree_indices <- which(rownames(pca_res$x) %in% share_trees)
  
  # Calculate the number of shared genes to display in the top right corner
  occup <- length(share_trees)
  
  # Property name for the legend
  file_id <- paste(column, "MIN_500", sep = "_")
  
  # Define the color and label vectors
  labels <- c("Center_500", file_id, paste('Shared', occup, sep = "_"))
  colors <- c("gold", "blue", "red")  # Colors for the plot
  
  # Create the plot and save it
  current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
    geom_point(color = 'black', size = 0.3, alpha = 1.0) +
    geom_point(data = pca_res$x[center_indices,], aes(PC1, PC2, color = labels[1]), size = 0.3, show.legend = TRUE) +
    geom_point(data = pca_res$x[top500_indices,], aes(PC1, PC2, color = labels[2]), size = 0.3, show.legend = TRUE) +
    geom_point(data = pca_res$x[sharetree_indices,], aes(PC1, PC2, color = labels[3]), size = 0.3, show.legend = TRUE) +
    scale_color_manual(labels = labels, values = colors) +
    labs(color = " ") +
    theme(legend.position = "top",
          panel.background = element_rect(fill = "grey"),
          plot.title = element_text(size = 10),
          axis.text = element_text(size = 10),  # Preserve the axis tick labels for small plots
          legend.text = element_text(size = 10),  # Adjust the font size of the legend text
          axis.title = element_text(size = 10),  # Set the axis title size
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    xlab("PC1 (73.7%)") +
    ylab("PC2 (3.2%)")  # Set the axis labels
  
  # Add the plot to the list
  plot_list[[i]] <- current_plot
  i <- i + 1
}


#select the 500 genes with medium evolutionary rate

# Calculate the median value
median_ER <- median(properties$ER)

# Compute the absolute differences from the median
properties$abs_diff_from_median <- abs(properties$ER - median_ER)

# Sort the data frame by the absolute differences
sorted_properties <- properties[order(properties$abs_diff_from_median), ]

# Select the top 500 values closest to the median
values_median_500 <- sorted_properties[1:500, ]

# Extract the name for the selected trees
name_median_500 <- values_median_500$rownames.tree

#save
#output <- paste0(csv_dir,"ER_IM_500.csv")
#write.csv(name_median_500, file = output)

# Find the indices of these trees in the PCA results
indices_median_500 <- which(rownames(pca_res$x) %in% name_median_500)

#find the shared trees
share_trees <- Center_500[Center_500 %in% name_median_500]

# Extract indices of shared genes in the PCA analysis results
sharetree_indices <- which(rownames(pca_res$x) %in% share_trees)

# Number of shared genes to be displayed in the top right corner
occup <- length(share_trees)

# Prepare property name for legend
file_id <- paste( "ER_IM_500", sep = "_" )

# Define color vector and label vector
labels <- c("Center_500",file_id, paste('Shared', occup, sep = "_" ))
colors <- c("gold", "blue", "red")   ##E95146

# Construct and save the plot
current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
  geom_point(color = 'black',size = 0.3, alpha = 1.0) +
  geom_point(data = pca_res$x[center_indices,], aes(PC1, PC2, color = labels[1]), size = 0.3, show.legend = TRUE) +
  geom_point(data = pca_res$x[indices_median_500,], aes(PC1, PC2, color = labels[2]), size = 0.3, show.legend = TRUE) +
  geom_point(data = pca_res$x[sharetree_indices,], aes(PC1, PC2, color = labels[3]), size = 0.3, show.legend = TRUE) +
  #manually specify the colors and labels for the legend in a plot
  scale_color_manual(labels = labels,values = colors) +
  labs(color = " ") +
  theme(legend.position = "top",
        panel.background = element_rect(fill = "grey"),
        plot.title = element_text(size = 10),
        axis.text = element_text(size = 10),  # Preserve axis tick labels for small plots
        legend.text = element_text(size = 10),  # Adjust the font size of legend text
        axis.title = element_text(size = 10),   # Set axis title size
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlab("PC1 (73.7%)") +
  ylab("PC2 (3.2%)")  # Set the axis labels

# Add the current plot to the list
plot_list[[i]] <- current_plot
i <- i + 1

# Combine all small plots into one large plot
combined_plot <- plot_grid(plotlist = plot_list, ncol = 3, align = "h")                                                                

# Save the combined large plot
output <- paste0(fig_dir,"occupation_500.tiff")
ggsave(plot = combined_plot, output, width = 25, height = 23, dpi = 300)

# Save the workspace to a file
output <- paste0(Rdata_dir, "5_occupation_500.RData") 
save.image(file = output)
