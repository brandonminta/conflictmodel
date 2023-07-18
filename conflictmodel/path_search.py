# -------------------------------------------------------------------------
# Created on 10/03/2023 19:21
# @author: Brandon Minta
# -------------------------------------------------------------------------
#
#   This module search for coalitions and neighbors in a 2d array
#
# --------------------------------------------------------------------------

import numpy as np

def vicinal(positions, i, L): 
    """
    This function extract the 4 neighbors of actor i.
    
    Parameters
    -----------
    positions: dictionary
        Positional indexes for each actor
    i: int
        Main actor
    L: int
        N = LxL
        
    return    
    ------
    neighbors: list
        List of neighbors
        
    """
    # Initialize the list of neighbors
    neighbors = []
    row, col = positions[i]
    
    # Check and append neighboring positions if they exist in the dictionary
    for indexes in [((row - 1) % L, col), ((row + 1) % L, col),
                 (row, (col - 1) % L), (row, (col + 1) % L)]:
        for key, value in positions.items():
            if value == indexes:
                neighbors.append(key)

    return neighbors

def filt(topography, positions, neighbors, coalition_list, value): 
    """
    This function filters the neighbors
    
    Parameters
    -----------
    topography: array
        Matrix of the topographical space  
    positions: dictionary
        actors indexes 
    neighbors: list
        List of the neighbors
    coalition_list: list
        Coalition list
    value: int
        value identifying a certain coalition

        
    return    
    ------
    vicinal_list_filtered: list
        List of the remained neighbors
    c: list
        Updated coalition list
        
    """
    neighbors_filtered  = []
    for i in neighbors:
        row, col = positions[i]
        if (topography[row, col] == value) and (i  not in coalition_list):
            coalition_list.append(i)
            neighbors_filtered.append(i)
    
    return neighbors_filtered, coalition_list

def group(i, topography, value, positions):
    """
    This function computes the coalition of i
    
    Parameters
    -----------
    i: int
        Actor to whom we will generate its coalition
    topography: array
        Matrix of the topographical space 
    val: int
        Value identifying a certain coalition    
    return    
    ------
    coalition: list
        List of i's coalition

    """
    L = topography.shape[0]
    neighbors = vicinal(positions, i, L)
    coalition_list = [i]
    
    #filterd neighbors and updated coalitions
    neighbors_filtered, coalition_list = filt(topography, positions, neighbors, coalition_list, value) 
    while len(neighbors_filtered)>0:
        temporal_neighbors = []
        #print(neighbors_filtered)
        for j in neighbors_filtered:
            local_neighbors = vicinal(positions, j, L)
         #   print(j, local_neighbors)
            new_local, coalition_list = filt(topography, positions, local_neighbors, coalition_list, value)
          #  print(new_local, coalition_list)
            if len(new_local)>0:
                temporal_neighbors = temporal_neighbors + new_local
        neighbors_filtered = temporal_neighbors
    return coalition_list 