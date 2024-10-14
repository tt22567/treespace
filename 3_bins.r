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
library(ape)
library(ggplot2)
library(cowplot)  #plot_grid

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
            replacement = "csv",
            x = full_path_to_file )

# find and replace characters in your path
setwd( wd )

#save directory
fig_dir <- gsub( pattern = "Scripts",
            replacement = "Figs/bins/",
            x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                  replacement = "Rdata/",
                  x = full_path_to_file )
  
#-------#
# TASKS #
#-------#

#---------------------------#
##   bins for properties   ##
##         single          ##
#---------------------------#

#load distance matrix
input <- paste0(Rdata_dir,'1_TreeDistance.Rdata')
load(input)

#distance matrix
res
length(rownames(res))
head(rownames(res))

#PCA analysis result for dimension reduction
pca_res

#read property matrix
properties <- read.csv("properties_final_2041.csv", header = TRUE, row.names = 'rownames.tree')
length(rownames(properties))
head(rownames(properties))

prop_rownames <- rownames(properties)
prop_colnames <- colnames(properties)

#properties_col <- colnames(properties[, c(2:16, 21)])
  
# Check if the row names are the same
if (all.equal(sort(rownames(res)), sort(rownames(properties)))) {
  print("Row names are the same")
} else {
  print("Row names are different")
}
#"Row names are the same"

#classify based on attributes
prop_positive_desc <- c("NS","AL","PPIS","NVS","LHF","ABSS","TTL","Treeness","TreenessdeRCV","PC2")
prop_negative_asc <- c("PG","RCV","LHU","ER","ALBS","Saturation") 

# Loop through each set of genes (in ascending order)
for (colname in prop_positive_desc) {
  
  # Sort the dataframe in ascending order by colname
  properties <- properties[order(properties[[colname]]), ]
  
  # Divide the genes into 6 groups, with the remainder added to the last group
  n <- nrow(properties)
  group_size <- floor(n / 6)   # round down, e.g., 2.3 to 2
  remainder <- n %% 6   # remainder
  
  # Create an empty list to store the generated plots
  plots_list <- list()
  
  # Loop to process each gene group
  for (i in 1:6) {
    start_index <- (i - 1) * group_size + 1
    end_index <- start_index + group_size - 1
    
    if (i == 6) {
      end_index <- start_index + group_size + remainder - 1
    }
    
    # Get the range of gene values for the current gene group
    tree_group <- rownames(properties)[start_index:end_index]
    
    # Extract the positions of genes in the PCA analysis results
    genes_fully_indices <- which(rownames(pca_res$x) %in% tree_group)
    
    # Get the maximum and minimum values for the current gene group
    max_value <- max(properties[[colname]][start_index:end_index])
    min_value <- min(properties[[colname]][start_index:end_index])
    
    # Generate labels for the current group
    labels <- paste(min_value, "-", max_value)
    
    # Generate a dataframe with markers
    plot_df <- data.frame(PC1 = pca_res$x[genes_fully_indices, 1],
                          PC2 = pca_res$x[genes_fully_indices, 2],
                          value_range = factor(labels, levels = labels))
    
    # Generate a scatter plot, marking the maximum and minimum values
    current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
      geom_point(size = 1.0, alpha = 1.0) +
      geom_point(data = plot_df, aes(PC1, PC2, color = value_range), size = 1.0, show.legend = TRUE) +
      # Manually specify the colors and labels for the legend in a plot
      scale_color_manual(values = "blue", labels = labels) +
      labs(color = colname) +
      theme(legend.position = "top",
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 15),  # Preserve axis tick labels for small plots
            legend.text = element_text(size = 15),  # Adjust the font size of the legend text
            legend.title = element_text(size = 15),  # Adjust the font size of the legend title
            axis.title = element_text(size = 15)) +   # Set axis title size
      xlab("PC1 (73.7%)") +
      ylab("PC2 (3.2%)")   # Set axis labels
    
    # Add the current plot to the list
    plots_list[[i]] <- current_plot
  }
  
  # Combine all small plots into one large plot
  combined_plot <- plot_grid(plotlist = plots_list, ncol = 3, align = "h")
  
  # Save the combined large plot
  output = paste0(fig_dir, colname, "_bins.tiff")
  ggsave(plot = combined_plot, filename = output, width = 18, height = 6, units = "in", dpi = 300)
}


# Loop through each set of genes (in descending order)
for (colname in prop_negative_asc) {
  
  # Sort the dataframe in descending order by colname
  properties <- properties[order(properties[[colname]], decreasing = TRUE), ]
  
  # Divide the genes into 6 groups, with the remainder added to the last group
  n <- nrow(properties)
  group_size <- floor(n / 6)   # round down, e.g., 2.3 to 2
  remainder <- n %% 6   # remainder
  
  # Create an empty list to store the generated plots
  plots_list <- list()
  
  # Loop to process each gene group
  for (i in 1:6) {
    start_index <- (i - 1) * group_size + 1
    end_index <- start_index + group_size - 1
    
    if (i == 6) {
      end_index <- start_index + group_size + remainder - 1
    }
    
    # Get the range of gene values for the current gene group
    tree_group <- rownames(properties)[start_index:end_index]
    
    # Extract the positions of genes in the PCA analysis results
    genes_fully_indices <- which(rownames(pca_res$x) %in% tree_group)
    
    # Get the maximum and minimum values for the current gene group
    max_value <- max(properties[[colname]][start_index:end_index])
    min_value <- min(properties[[colname]][start_index:end_index])
    
    # Generate labels for the current group
    labels <- paste(min_value, "-", max_value)
    
    # Generate a dataframe with markers
    plot_df <- data.frame(PC1 = pca_res$x[genes_fully_indices, 1],
                          PC2 = pca_res$x[genes_fully_indices, 2],
                          value_range = factor(labels, levels = labels))
    
    # Generate a scatter plot, marking the maximum and minimum values
    current_plot <- ggplot(data = pca_res$x, aes(PC1, PC2)) +
      geom_point(size = 1.0, alpha = 1.0) +
      geom_point(data = plot_df, aes(PC1, PC2, color = value_range), size = 1.0, show.legend = TRUE) +
      # Manually specify the colors and labels for the legend in a plot
      scale_color_manual(values = "blue", labels = labels) +
      labs(color = colname) +
      theme(legend.position = "top",
            plot.title = element_text(size = 15),
            axis.text = element_text(size = 15),  # Preserve axis tick labels for small plots
            legend.text = element_text(size = 15),  # Adjust the font size of the legend text
            legend.title = element_text(size = 15),  # Adjust the font size of the legend title
            axis.title = element_text(size = 15)) +   # Set axis title size
      xlab("PC1 (73.7%)") +
      ylab("PC2 (3.2%)")   # Set axis labels
    
    # Add the current plot to the list
    plots_list[[i]] <- current_plot
  }
  
  # Combine all small plots into one large plot
  combined_plot <- plot_grid(plotlist = plots_list, ncol = 3, align = "h")
  
  # Save the combined large plot
  output = paste0(fig_dir, colname, "_bins.tiff")
  ggsave(plot = combined_plot, filename = output, width = 18, height = 6, units = "in", dpi = 300)
}

# Save the data
output <- paste0(Rdata_dir, "3_bins.RData")
save.image(file = output)
