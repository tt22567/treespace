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
library(colorspace) 
library(dplyr)

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
                 replacement = "Figs/",
                 x = full_path_to_file )

#Rdata directory
Rdata_dir <- gsub( pattern = "Scripts",
                   replacement = "Rdata/",
                   x = full_path_to_file )

#-------#
# TASKS #
#-------#

#---------------------------------------------------#
#             choose central X trees                #
#             convex hull the groups                #
#---------------------------------------------------#

input <- paste0(Rdata_dir, "15_TreeDistance_2062_NMDS.RData")
load(input)
res
nmds2

#read the distance matrix
dist_mat <- res

# Calculate the sum of each row
row_sums <- rowSums(dist_mat)
min(row_sums)
max(row_sums)

# Sort the row sums in descending order
sorted_row_sums <- sort(row_sums, decreasing = FALSE)

# Median tree
median_tree <- names(sorted_row_sums[1:1])
gene_indices_median_tree <- which(rownames(nmds2$points) %in% median_tree)

# Extract the top n rows
top_100_row_sums <- sorted_row_sums[1:100]
top_500_row_sums <- sorted_row_sums[1:500]
top_1000_row_sums <- sorted_row_sums[1:1000]
top_1500_row_sums <- sorted_row_sums[1:1500]
top_1800_row_sums <- sorted_row_sums[1:1800]
#top_2062_row_sums <- sorted_row_sums[1:2062]

# Print the row names of the top n rows
Center_100 <- names(top_100_row_sums)
Center_500 <- names(top_500_row_sums)
Center_1000 <- names(top_1000_row_sums)
Center_1500 <- names(top_1500_row_sums)
Center_1800 <- names(top_1800_row_sums)
#Center_2062 <- names(top_2062_row_sums)

# Find their position
gene_indices_Center_100 <- which(rownames(nmds2$points) %in% Center_100)
gene_indices_Center_500 <- which(rownames(nmds2$points) %in% Center_500)
gene_indices_Center_1000 <- which(rownames(nmds2$points) %in% Center_1000)
gene_indices_Center_1500 <- which(rownames(nmds2$points) %in% Center_1500)
gene_indices_Center_1800 <- which(rownames(nmds2$points) %in% Center_1800)
#gene_indices_Center_2062 <- which(rownames(nmds2$points) %in% Center_2062)

#-------------------------------------------------------------#
#       convex hull the groups and highlight two trees        #
#-------------------------------------------------------------#

gene_indices_closest_1000 <- which(rownames(nmds2$points) %in% 'Closest_1000.treefile')
gene_indices_Astral_1000 <- which(rownames(nmds2$points) %in% 'Astral_output.treefile')

# Function to get convex hull
get_convex_hull <- function(df) {
  df <- as.data.frame(df)  # Ensure df is a data frame
  df[chull(df$MDS1, df$MDS2), ]
}

# Ensure nmds2$points is a data frame
nmds2$points <- as.data.frame(nmds2$points)

# Calculate convex hulls for each group
#hull_2062 <- get_convex_hull(nmds2$points[gene_indices_Center_2062, ])
hull_1800 <- get_convex_hull(nmds2$points[gene_indices_Center_1800, ])
hull_1500 <- get_convex_hull(nmds2$points[gene_indices_Center_1500, ])
hull_1000 <- get_convex_hull(nmds2$points[gene_indices_Center_1000, ])
hull_500 <- get_convex_hull(nmds2$points[gene_indices_Center_500, ])
hull_100 <- get_convex_hull(nmds2$points[gene_indices_Center_100, ])

# Combine all hulls into a single data frame with a group identifier
hull_df <- dplyr::bind_rows(
  #mutate(hull_2062, group = "Center_2062"),
  mutate(hull_1800, group = "Center_1800"),
  mutate(hull_1500, group = "Center_1500"),
  mutate(hull_1000, group = "Center_1000"),
  mutate(hull_500, group = "Center_500"),
  mutate(hull_100, group = "Center_100")
)

# Calculate the size of the points; the smaller the row sums, the smaller the points
sizes <- (row_sums + 1)  # Add a constant 1 to avoid extreme values due to zero

# Determine the minimum and maximum values of sizes
min_size <- min(sizes)  
max_size <- max(sizes)  

# Ensure nmds2$points is a data frame and contains the size information
nmds2$points$Sum_Distance <- sizes

# Define colors for density points and convex hull
colors <- colorspace::sequential_hcl(n = 5, palette = "ag_Sunset") #viridis, ag_Sunset

# Plotting
plot <- ggplot(data = nmds2$points, aes(MDS1, MDS2)) +
  theme_gray() +  # default grey theme
  geom_point(aes(size = Sum_Distance), alpha = 1.0, color = "black") +  # Map point size based on Sum_Distance
  scale_size_continuous(range = c(2.0, 10.0),  #5 times
                        breaks = c(min_size, max_size), 
                        labels = c(min_size, max_size)) +  # Adjust size range and legend
  
  # Add hulls and group-specific points
  geom_polygon(data = hull_df, aes(x = MDS1, y = MDS2, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  geom_point(data = nmds2$points[gene_indices_Center_1800, ], aes(MDS1, MDS2, color = "Center_1800", size = Sum_Distance)) +
  geom_point(data = nmds2$points[gene_indices_Center_1500, ], aes(MDS1, MDS2, color = "Center_1500", size = Sum_Distance)) +
  geom_point(data = nmds2$points[gene_indices_Center_1000, ], aes(MDS1, MDS2, color = "Center_1000", size = Sum_Distance)) +
  geom_point(data = nmds2$points[gene_indices_Center_500, ], aes(MDS1, MDS2, color = "Center_500", size = Sum_Distance)) +
  geom_point(data = nmds2$points[gene_indices_Center_100, ], aes(MDS1, MDS2, color = "Center_100", size = Sum_Distance)) + 
  geom_point(data = nmds2$points[gene_indices_closest_1000, ], aes(MDS1, MDS2, shape = "Center_1000"), color = "red", size = 8, inherit.aes = FALSE) +
  geom_point(data = nmds2$points[gene_indices_Astral_1000, ], aes(MDS1, MDS2, shape = "Astral"), color = "blue", size = 8, inherit.aes = FALSE) +
  geom_point(data = nmds2$points[gene_indices_median_tree, ], aes(MDS1, MDS2, shape = "Median"), color = "purple", size = 8, inherit.aes = FALSE) +

  
  scale_color_manual(values = setNames(c(colors, "red", "blue",'purple'), c(
    "Center_1800", "Center_1500", 
    "Center_1000", "Center_500", "Center_100", "closest_1000", "Astral",'Median'
  )), breaks = c("Center_100", "Center_500", "Center_1000", "Center_1500", "Center_1800"), guide = guide_legend(title = "Group")) +
  
  scale_fill_manual(values = setNames(colors, c(
    "Center_1800", "Center_1500", 
    "Center_1000", "Center_500", "Center_100"
  )), breaks = c("Center_100", "Center_500", "Center_1000", "Center_1500", "Center_1800"), guide = guide_legend(title = "Group")) +
  
  scale_shape_manual(values = c("Center_1000" = 8, "Astral" = 7, 'Median' = 3), guide = guide_legend(title = "Tree"))+
  theme(
    legend.position = "right",
    plot.title = element_text(size = 28),
    axis.text = element_text(size = 28),  # Preserve axis tick labels for small plots
    legend.text = element_text(size = 28),  # Adjust the font size of the legend text
    legend.title = element_text(size = 28),  # Adjust the font size of the legend title
    axis.title = element_text(size = 28)  # Set axis title size
  ) +
  xlab("MDS1") +
  ylab("MDS2") +  # Set axis labels
  labs(fill = "Group", color = "Group", shape = "Tree", size = "Sum_Distance")  # Add titles to legends

print(plot)

# Save
output <- paste0(fig_dir, "Fig1_TreespacewithConvexHull_NMDS1.tiff")
ggsave(plot, file = output, width = 18, height = 14, dpi = 300)

#---------------------------------------------------#
#      Trees you might want to show on the plot     #
#---------------------------------------------------#

# Define gene indices for all tree files
gene_indices_closest_1000 <- which(rownames(nmds2$points) %in% 'Closest_1000.treefile')
gene_indices_Astral_1000 <- which(rownames(nmds2$points) %in% 'Astral_output.treefile')
gene_indices_AL_MAX_1000 <- which(rownames(nmds2$points) %in% 'AL_MAX_1000.treefile')
gene_indices_PG_MIN_1000 <- which(rownames(nmds2$points) %in% 'PG_MIN_1000.treefile')
gene_indices_MDS2_MAX_1000 <- which(rownames(nmds2$points) %in% 'MDS2_MAX_1000.treefile')
gene_indices_Randomness_MAX_1000 <- which(rownames(nmds2$points) %in% 'Randomness_MAX_1000.treefile')
gene_indices_PPIS_MAX_1000 <- which(rownames(nmds2$points) %in% 'PPIS_MAX_1000.treefile')
gene_indices_RCV_MIN_1000 <- which(rownames(nmds2$points) %in% 'RCV_MIN_1000.treefile')
gene_indices_LHF_MAX_1000 <- which(rownames(nmds2$points) %in% 'LHF_MAX_1000.treefile')
gene_indices_LHU_MIN_1000 <- which(rownames(nmds2$points) %in% 'LHU_MIN_1000.treefile')
gene_indices_ABSS_MAX_1000 <- which(rownames(nmds2$points) %in% 'ABSS_MAX_1000.treefile')
gene_indices_ER_IM_1000 <- which(rownames(nmds2$points) %in% 'ER_IM_1000.treefile')
gene_indices_ER_MAX_1000 <- which(rownames(nmds2$points) %in% 'ER_MAX_1000.treefile')
gene_indices_ER_MIN_1000 <- which(rownames(nmds2$points) %in% 'ER_MIN_1000.treefile')
gene_indices_ALBS_MIN_1000 <- which(rownames(nmds2$points) %in% 'ALBS_MIN_1000.treefile')
gene_indices_TTL_MAX_1000 <- which(rownames(nmds2$points) %in% 'TTL_MAX_1000.treefile')
gene_indices_Treeness_MAX_1000 <- which(rownames(nmds2$points) %in% 'Treeness_MAX_1000.treefile')
gene_indices_Saturation_MIN_1000 <- which(rownames(nmds2$points) %in% 'Saturation_MIN_1000.treefile')
gene_indices_TreenessdeRCV_MAX_1000 <- which(rownames(nmds2$points) %in% 'TreenessdeRCV_MAX_1000.treefile')
gene_indices_NVS_MAX_1000 <- which(rownames(nmds2$points) %in% 'NVS_MAX_1000.treefile')
gene_indices_NS_MAX_1000 <- which(rownames(nmds2$points) %in% 'NS_MAX_1000.treefile')

#---------------------------------------------------#
#                    Zoom in                        #
#        Define colors, sizes and shapes            #
#---------------------------------------------------#

# Define colors for density points
#colors <- c('#082a54','#a00000','#298c8c','#ea801c','#a559aa')
colors <- colorspace::sequential_hcl(n = 5, palette = "ag_Sunset") #viridis, ag_Sunset

# Define colors for tree
#cols <- colorRampPalette(c("blue", "yellow"))(17)
cols <- colorspace::sequential_hcl(n = 19, palette = "Plasma")

# Define shape mapping
tree_shapes <- c(
  "Closest_1000" = 8, 
  "Astral_output" = 7, 
  "AL_MAX_1000" = 2,
  "PG_MIN_1000" = 3,
  "MDS2_MAX_1000" = 4,
  "Randomness_MAX_1000" = 5,
  "PPIS_MAX_1000" = 6,
  "RCV_MIN_1000" = 7,
  "LHF_MAX_1000" = 8,
  "LHU_MIN_1000" = 9,
  "ABSS_MAX_1000" = 10,
  "ER_IM_1000" = 11,
  "ER_MAX_1000" = 12,
  "ER_MIN_1000" = 13,
  "ALBS_MIN_1000" =14,
  "TTL_MAX_1000" = 15,
  "Treeness_MAX_1000" = 2,
  "Saturation_MIN_1000" = 8,
  "TreenessdeRCV_MAX_1000" = 9,
  'NVS_MAX_1000'= 15,
  'NS_MAX_1000' = 8
)

#---------------------------------------------------#
#                    Zoom in                        #
#               Plot tree space                     #
#---------------------------------------------------#

# Define the PCA plot with the original points
plot <- ggplot(data = nmds2$points, aes(MDS1, MDS2)) +
  theme_gray() +  # default grey theme
  geom_point(aes(size = Sum_Distance), alpha = 1.0, color = "black") +  # Map point size based on Distance
  scale_size_continuous(range = c(1, 5), 
                        breaks = c(min_size, max_size), 
                        labels = c(min_size, max_size)) +  # Adjust size range and legend

  # Add the specific gene indices points
  geom_point(data = nmds2$points[gene_indices_Center_1800, ], aes(MDS1, MDS2, size = Sum_Distance, color = "Center_1800"), alpha = 1.0) +
  geom_point(data = nmds2$points[gene_indices_Center_1500, ], aes(MDS1, MDS2, size = Sum_Distance, color = "Center_1500"), alpha = 1.0) +
  geom_point(data = nmds2$points[gene_indices_Center_1000, ], aes(MDS1, MDS2, size = Sum_Distance, color = "Center_1000"), alpha = 1.0) +
  geom_point(data = nmds2$points[gene_indices_Center_500, ], aes(MDS1, MDS2, size = Sum_Distance, color = "Center_500"), alpha = 1.0) +
  geom_point(data = nmds2$points[gene_indices_Center_100, ], aes(MDS1, MDS2, size = Sum_Distance, color = "Center_100"), alpha = 1.0) + 
  geom_point(data = nmds2$points[gene_indices_closest_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "Closest_1000"), color = "red") +
  geom_point(data = nmds2$points[gene_indices_Astral_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "Astral_output"), color = "blue") +
  geom_point(data = nmds2$points[gene_indices_AL_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "AL_MAX_1000"), color = cols[1]) +
  geom_point(data = nmds2$points[gene_indices_PG_MIN_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "PG_MIN_1000"), color = cols[2]) +
  geom_point(data = nmds2$points[gene_indices_MDS2_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "MDS2_MAX_1000"), color = cols[3]) +
  geom_point(data = nmds2$points[gene_indices_Randomness_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "Randomness_MAX_1000"), color = cols[4]) +
  geom_point(data = nmds2$points[gene_indices_PPIS_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "PPIS_MAX_1000"), color = cols[5]) +
  geom_point(data = nmds2$points[gene_indices_RCV_MIN_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "RCV_MIN_1000"), color = cols[6]) +
  geom_point(data = nmds2$points[gene_indices_LHF_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "LHF_MAX_1000"), color = cols[7]) +
  geom_point(data = nmds2$points[gene_indices_LHU_MIN_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "LHU_MIN_1000"), color = cols[8]) +
  geom_point(data = nmds2$points[gene_indices_ABSS_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "ABSS_MAX_1000"), color = cols[9]) +
  geom_point(data = nmds2$points[gene_indices_ER_IM_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "ER_IM_1000"), color = cols[10]) +
  geom_point(data = nmds2$points[gene_indices_ER_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "ER_MAX_1000"), color = cols[11]) +
  geom_point(data = nmds2$points[gene_indices_ER_MIN_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "ER_MIN_1000"), color = cols[12]) +
  geom_point(data = nmds2$points[gene_indices_ALBS_MIN_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "ALBS_MIN_1000"), color = cols[13]) +
  geom_point(data = nmds2$points[gene_indices_TTL_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "TTL_MAX_1000"), color = cols[14]) +
  geom_point(data = nmds2$points[gene_indices_Treeness_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "Treeness_MAX_1000"), color = cols[15]) +
  geom_point(data = nmds2$points[gene_indices_Saturation_MIN_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "Saturation_MIN_1000"), color = cols[16]) +
  geom_point(data = nmds2$points[gene_indices_TreenessdeRCV_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "TreenessdeRCV_MAX_1000"), color = cols[17]) +
  geom_point(data = nmds2$points[gene_indices_NVS_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "NVS_MAX_1000"), color = cols[18]) +
  geom_point(data = nmds2$points[gene_indices_NS_MAX_1000, ], aes(MDS1, MDS2, size = Sum_Distance, shape = "NS_MAX_1000"), color = cols[19]) +
  
  scale_color_manual(values = setNames(colors, c(
    #"Center_2062", 
    "Center_1800", "Center_1500", 
    "Center_1000", "Center_500", "Center_100"
  )), guide = guide_legend(title = "Density")) +
  scale_fill_manual(values = setNames(colors, c(
    #"Center_2062", 
    "Center_1800", "Center_1500", 
    "Center_1000", "Center_500", "Center_100"
  )), guide = guide_legend(title = "Density")) +
  scale_shape_manual(values = tree_shapes, guide = guide_legend(title = "Tree", ncol = 2)) +
  guides(size = guide_legend(override.aes = list(color = "black"))) +  # Set color to white in the Distance label
  theme(legend.position = "right",
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))+
  xlab("MDS1") +
  ylab("MDS2") +  # Set axis labels
  labs(fill = "Group", color = "Group", shape = "Tree")

# Print the plot
print(plot)

# Zoom in on a specific region, e.g., x from -10 to 10 and y from -5 to 5
p <- plot + coord_cartesian(xlim = c(-0.012, -0.0075), ylim = c(-0.005, 0.0025))

# Print the plot
print(p)

# Save
output_pdf <- paste0(fig_dir, "Fig1_zoomin_NMDS.pdf")
ggsave(p, file = output_pdf, width = 15, height = 6)

# Close
dev.off()

