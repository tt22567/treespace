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
library(cowplot)

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
                 replacement = "Figs/why1000/",
                 x = full_path_to_file )


#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )

#-------#
# TASKS #
#-------#

#--------------------------------#
# why choose 1000 gene trees     #
#--------------------------------#

#load distance matrix
input <- paste0(Rdata_dir,'1_TreeDistance.Rdata')
load(input)

#distance matrix
res
length(rownames(res))
head(rownames(res))

#PCA analysis result for dimension reduction
pca_res

# read property matrix
properties <- read.csv("properties_final_2041.csv", header = TRUE, row.names = 'rownames.tree')

#classify based on attributes
prop_positive_desc <- c("NS","AL","PPIS","NVS","LHF","ABSS","TTL","Treeness","TreenessdeRCV","PC2")
prop_negative_asc <- c("PG","RCV","LHU","ER","ALBS","Saturation") 

# Define tree numbers to plot in each plot
max_values <- c(1800, 1600, 1400, 1200, 1000, 800, 600, 400, 200)

# Loop through each target column and each max value
for (target_col in prop_positive_desc) {
  
  # Create an empty list to store the generated plots
  plots_list <- list()
  
  # Get target column values
  target_prop <- properties[[target_col]]
  
  # Initialize a counter for the plots_list
  plot_index <- 1
  
  # Loop through each max value
  for (max_value in max_values) {
    
    # Get top indices for current max value
    top_indices <- order(target_prop, decreasing = TRUE)[1:max_value]
    
    # Get genes fully
    genes_fully <- rownames(properties)[top_indices]
    
    # Get genes fully indices
    genes_fully_indices <- which(rownames(pca_res$x) %in% genes_fully)
    
    # Generate labels for the current group
    max_label <- paste("MAX_", max_value)
    
    # Generate a dataframe with markers
    plot_df <- data.frame(PC1 = pca_res$x[genes_fully_indices, 1],
                          PC2 = pca_res$x[genes_fully_indices, 2],
                          Group = factor(max_label, levels = max_label))
    
    # Create current plot
    current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
      geom_point(size = 1.0, alpha = 1.0) +
      geom_point(data = plot_df, aes(PC1, PC2, color = Group), size = 1.0) +
      scale_color_manual(values = "blue", labels = max_label) +
      labs(color = target_col) +
      theme(legend.position = "top",
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 15),  # Preserve axis tick labels for small plots
            legend.text = element_text(size = 15),  # Adjust the font size of the legend text
            legend.title = element_text(size = 15),  # Adjust the font size of the legend title
            axis.title = element_text(size = 15)) +   # Set axis title size
      xlab("PC1 (73.7%)") +
      ylab("PC2 (3.2%)")   # Set axis labels
    
    # Add the current plot to the list
    plots_list[[plot_index]] <- current_plot
    
    # Increment the plot index
    plot_index <- plot_index + 1
  }
  
  # Combine all small plots into one large plot
  combined_plot <- plot_grid(plotlist = plots_list, ncol = 3, align = "h")
  
  # Save the combined large plot
  output = paste0(fig_dir, target_col, "_why1000.tiff")
  ggsave(plot = combined_plot, filename = output, width = 18, height = 9, units = "in", dpi = 300)
}


# Loop through each target column and each max value
for (target_col in prop_negative_asc) {
  
  # Create an empty list to store the generated plots
  plots_list <- list()
  
  # Get target column values
  target_prop <- properties[[target_col]]
  
  # Initialize a counter for the plots_list
  plot_index <- 1
  
  # Loop through each max value
  for (max_value in max_values) {
    
    # Get top indices for current max value (sorted in ascending order)
    top_indices <- order(target_prop, decreasing = FALSE)[1:max_value]
    
    # Get genes fully
    genes_fully <- rownames(properties)[top_indices]
    
    # Get genes fully indices
    genes_fully_indices <- which(rownames(pca_res$x) %in% genes_fully)
    
    # Generate labels for the current group
    max_label <- paste("MIN_", max_value)
    
    # Generate a dataframe with markers
    plot_df <- data.frame(PC1 = pca_res$x[genes_fully_indices, 1],
                          PC2 = pca_res$x[genes_fully_indices, 2],
                          Group = factor(max_label, levels = max_label))
    
    # Create current plot
    current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
      geom_point(size = 1.0, alpha = 1.0) +
      geom_point(data = plot_df, aes(PC1, PC2, color = Group), size = 1.0) +
      scale_color_manual(values = "blue", labels = max_label) +
      labs(color = target_col) +
      theme(legend.position = "top",
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 15),  # Preserve axis tick labels for small plots
            legend.text = element_text(size = 15),  # Adjust the font size of the legend text
            legend.title = element_text(size = 15),  # Adjust the font size of the legend title
            axis.title = element_text(size = 15)) +   # Set axis title size
      xlab("PC1 (73.7%)") +
      ylab("PC2 (3.2%)")   # Set axis labels
    
    # Add the current plot to the list
    plots_list[[plot_index]] <- current_plot
    
    # Increment the plot index
    plot_index <- plot_index + 1
  }
  
  # Combine all small plots into one large plot
  combined_plot <- plot_grid(plotlist = plots_list, ncol = 3, align = "h")
  
  # Save the combined large plot
  output = paste0(fig_dir, target_col, "_why1000.tiff")
  ggsave(plot = combined_plot, filename = output, width = 18, height = 9, units = "in", dpi = 300)
}
