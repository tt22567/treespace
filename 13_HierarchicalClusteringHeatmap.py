
##Xiumei Lu code
#https://seaborn.pydata.org/generated/seaborn.clustermap.html

#powershell
#pip install 

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy

tree_file = "C:\\Users\\tt22567\\OneDrive - University of Bristol\\Documents\\R\\CalculateDistance\\TreespaceAnalysis\\csv\\res_tree_distance_BMGE_Normolise_20.csv"
path_file = "C:\\Users\\tt22567\\OneDrive - University of Bristol\\Documents\\R\\CalculateDistance\\TreespaceAnalysis\\Figs\\"

# Read CSV
df = pd.read_csv(tree_file, index_col=0)

# Define plot size
plt.figure(figsize=(10, 10))

# Generate cluster map
# Color maps: 'plasma' (blue-yellow), 'viridis' (green-purple-blue-yellow), 'inferno' (red-black-yellow), 'vlag' (pink red and blue)
# Methods include 'single', 'complete', 'average', 'ward'
# Use df or df.corr() depending on your data

g = sns.clustermap(df, pivot_kws=None, method='average', metric='euclidean', cmap="vlag",
                   z_score=None, standard_scale=None, figsize=(10, 10), cbar_kws=None, 
                   row_cluster=True, col_cluster=True, row_linkage=None, col_linkage=None,  # Cluster rows and columns
                   row_colors=None, col_colors=None, 
                   mask=None, dendrogram_ratio=0.3, colors_ratio=0.03, # Dendrogram ratio
                   cbar_pos=(0.03, 0.06, 0.02, 0.25),  # Adjust color bar position (left, bottom, width, height)
                   tree_kws={
                       #'linestyles':'dashed', # Line style
                       'colors':'steelblue',  # Line color
                       'linewidths':1         # Line width
                   })

# Plot
#plt.show()

# Save
plt.savefig(path_file + "heatmap.tiff", format='tiff', dpi=300)
