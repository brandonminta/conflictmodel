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
from scipy.stats import linregress
from matplotlib.colors import LinearSegmentedColormap

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


def plot_actions(grid):
    labels = np.empty_like(grid, dtype=object)
    labels[labels == None] = ''

    # Find indices of zeros in the array
    zero_indices = np.argwhere(grid != 1)

    # Loop through the zero indices and set the corresponding tick labels to '0 to N'
    n = 0
    for idx in zero_indices:
        x, y = idx[0], idx[1]
        labels[x][y] = n
        n += 1

    # Define custom color maps
    cmap_blue = LinearSegmentedColormap.from_list("Blues", ['#FFFFFF', '#9999FF', '#6666FF', '#3333FF'])
    cmap_red = LinearSegmentedColormap.from_list("Reds", ['#FFFFFF', '#FF9999', '#FF6666', '#FF3333'])
   
    # Create the heatmap with annotations and custom colors for each value
    sns.heatmap(grid, cmap=cmap_blue, annot=labels, fmt='', cbar=False, square=True, linewidths=1.5, linecolor='black',
                mask=np.isin(grid, [2, 4]), cbar_kws={'ticks': []})

    sns.heatmap(grid, cmap=cmap_red, annot=labels, fmt='', cbar=False, square=True, linewidths=1.5, linecolor='black',
                mask=np.isin(grid, [3,6]), cbar_kws={'ticks': []})
    # Color cells with value 1 in black
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            if grid[i, j] == 1:
                plt.fill([j, j + 1, j + 1, j], [i, i, i + 1, i + 1], color='black')


    # Remove numerical ticks on x and y axes
    plt.xticks([])
    plt.yticks([])

    # Add title (optional)
    plt.title("Actor's Distribution")

    # Display the plot
    plt.show()


def distibution_plot(x,y, style, N, cycle, log_value = True, ccdf = False):
    plt.figure(figsize=(30, 10), facecolor='white')  # Adjust the values (width, height) as desired

    if ccdf == True:
        cdf = np.cumsum(y)
        ccdf_value = 1 - cdf
        y = ccdf_value
            # Plot the probability distribution on a log-log scale
        if log_value == True:
            plt.loglog(x, y, marker='o', linestyle='',color='#800080',  markersize=8)    
        else: 
            plt.plot(x, y, marker='o', linestyle='', color='#800080',  markersize=8)
            
                # Perform log-log transformation
        log_sorted_values = np.log(x)
        log_prob_distribution = np.log(y)

        # Perform linear regression to estimate the slope
        slope, intercept, r_value, p_value, std_err = linregress(log_sorted_values, log_prob_distribution)

        # Create the line that represents the best fit
        line = slope * log_sorted_values + intercept

        #Plot the line on the log-log plot
        plt.loglog(x, np.exp(line), color='red', label=f'Slope: {slope:.2f}')
        
        if style == 'd':
            plt.xlabel(r"Participants - $\tau$", fontsize=25)
            plt.ylabel(r"Probability $P(\tau)$", fontsize=25)
            plt.title("Complementary Cumulative Distribution Function (CCDF) Defence coalition Size", fontsize=25)
        elif style == 'a':
            plt.xlabel(r"Participants - $\alpha$", fontsize=25)
            plt.ylabel(r"Probability $P(\alpha)$", fontsize=25)
            plt.title("Complementary Cumulative Distribution Function (CCDF)Attacking coalition Size", fontsize=25)
        elif style == 't':
            plt.xlabel(r"Participants - $x$", fontsize=25)
            plt.ylabel(r"$P(X<x)$", fontsize=25)
            plt.title("Complementary Cumulative Distribution Function (CCDF) Conflict Size", fontsize=25)
        elif style =='p':
                # Set the labels and title of the plot
            plt.xlabel("x", fontsize=25)
            plt.ylabel("P(X<x)", fontsize=25)
            plt.title("Complementary Cumulative Distribution Function (CCDF) Peace Intervals", fontsize=25)
    # Plot the probability distribution on a log-log scale
    
    elif ccdf == False:
        if log_value == True:
            plt.loglog(x, y, marker='o', linestyle='', color='blue')    
        else: 
            plt.plot(x, y, marker='o', linestyle='', color='blue')

        log_sorted_values = np.log(x)
        log_prob_distribution = np.log(y)

        # Perform linear regression to estimate the slope
        slope, intercept, r_value, p_value, std_err = linregress(log_sorted_values, log_prob_distribution)

        # Create the line that represents the best fit
        line = slope * log_sorted_values + intercept

        #Plot the line on the log-log plot
        plt.loglog(x, np.exp(line), color='red', label=f'Slope: {slope:.2f}')

        # Increase the font size of the axis labels and tick labels
        if style == 'd':
            plt.xlabel(r"Participants - $\tau$", fontsize=25)
            plt.ylabel(r"Probability $P(\tau)$", fontsize=25)
            plt.title("Probability Distribution - Defence coalition Size", fontsize=25)
        elif style == 'a':
            plt.xlabel(r"Participants - $\alpha$", fontsize=25)
            plt.ylabel(r"Probability $P(\alpha)$", fontsize=25)
            plt.title("Probability Distribution - Attacking coalition Size", fontsize=25)
        elif style == 't':
            plt.xlabel(r"Participants - $x$", fontsize=25)
            plt.ylabel(r"Probability $P(x)$", fontsize=25)
            plt.title("Probability Distribution - Conflict Size", fontsize=25)
        elif style =='p':
                # Set the labels and title of the plot
            plt.xlabel("Peace Intervals - x", fontsize=25)
            plt.ylabel("P(x)", fontsize=25)
            plt.title("Probability Frequency of Peace Intervals", fontsize=25)

    # Add parameter information as a text annotation
    params_info = f"N: {N}\nλ: {N//cycle}"
    plt.text(0.5, 0.95, params_info, fontsize=18, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.8))


    # Increase the font size of the tick labels
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=19)

    plt.legend(fontsize=24)

    plt.grid(True)
    plt.show()
    
def bi_plot(x,y1,y2):
    fig, ax1 = plt.subplots(figsize=(27, 10), facecolor='white')
    ax2 = ax1.twinx()

    ax1.set_ylabel("Activity", fontsize=30)
    ax2.set_ylabel("Total Resources", fontsize=30)
    ax1.set_xlabel("Year", fontsize=30)

    ax1.tick_params(axis='both', labelsize=21)
    ax2.tick_params(axis='both', labelsize=21)

    ax1.plot(x, y1, linestyle=" ", marker="+", color="#FF4C4C", label='Activity')
    ax2.plot(x, y2, linestyle="--", color="#0055FF", label='Total Resources')

    # Create the legend
    ax1.legend(loc='upper left', fontsize=16)
    ax2.legend(loc='upper right', fontsize=16)

    # Add parameter information as a text annotation
    params_info = f"N: {N}\nλ: {N//3}"
    plt.text(0.5, 0.95, params_info, transform=fig.transFigure, fontsize=21, ha='center')

    plt.grid(True)
    plt.show()
    plt.close()