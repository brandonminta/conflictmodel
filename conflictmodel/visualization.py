# -------------------------------------------------------------------------
# Created on 10/03/2023 19:21
# @author: Brandon Minta
# -------------------------------------------------------------------------
#
#   This module search for coalitions and neighbors in a 2d array
#
# --------------------------------------------------------------------------
    # Create a figure and axis
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
    
def plot_commitment(data):
    grid_size = len(data)
    fig, ax = plt.subplots(figsize=(grid_size // 2, grid_size // 2))  # Adjust the division factor to control the size

    # Create the heatmap
    sns.heatmap(data, cmap='Greys', annot=True, fmt='.1f', square=True, ax=ax)

    # Add labels and title (optional)
    plt.xlabel("Actors")
    plt.ylabel("Actors")
    plt.title("Commitment matrix")

    # Display the plot
    plt.show()

    
def plot_grid(grid):
    labels = np.empty_like(grid, dtype=object)
    #labels[:] = np.nan
    labels[labels == None] = ''
    # Find indices of zeros in the array
    zero_indices = np.argwhere(grid == 0)

    # Loop through the zero indices and set the corresponding tick labels to '0 to N'
    n = 0
    for idx in zero_indices:
        x, y = idx[0], idx[1]
        labels[x][y] = n
        n +=1
        
    # Create the heatmap with annotations
    #sns.heatmap(agent_numerization, cmap='tab10', annot=grid, fmt='', cbar=False, square=True)
    sns.heatmap(grid, cmap='Greys', annot=labels, fmt='', cbar=False, square=True, linewidths=1.5, linecolor='black')

    # Remove numerical ticks on x and y axes
    plt.xticks([])
    plt.yticks([])

    # Add title (optional)
    plt.title("Actor's Distribution")

    # Display the plot
    plt.show()
