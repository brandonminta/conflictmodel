# -------------------------------------------------------------------------
# Created on 10/03/2023 19:21
# @author: Brandon Minta
# -------------------------------------------------------------------------
#
#   This module search for coalitions and neighbors in a 2d array
#
# --------------------------------------------------------------------------
import numpy as np
  
def group_resources(agent, coalition_arr, capital_arr, loyalty_mtx):                                                                                   
    """
    The function calculates the total resources of a coalition
    
    Parameters
    -----------
    agent: int
        The actor of interest 
    coalition_arr: array
        Coalition of the agent
    capital_arr: array
        Capital of each member
    loyalty_mtx: array
        Commitment matrix
    
    return
    ------
    money: float
        The total resources of the coalition
        
    """
    money = 0
    for i in coalition_arr:
        money += 0.1*loyalty_mtx[i][agent] * capital_arr[i]        
    return money   

def vulnerability(r_i, r_j):
    """
    The function computes the vulnerability of agent j with respect to agent i
    
    Parameters
    -----------
    r_i: float
        Resources of i's coalition
    r_j: float
        Resources of j's coalition
        
    Returns
    -------
    float
        The vulnerability of j with respect to i
    """
    return (r_i - r_j) / r_i if r_i > 0 and r_j > 0 else 0

def groups_matrix(attacker, target, indices, loyalty_mtx, grid) :
    """
    This function generates an array that represents commitments based on two groups
    
    Parameters
    -----------
    matrix: array matrix (NxN) 
        Commitment matrix
    grid: array matrix (LxL) 
        Geographical matrix
    indices: dictionary
            indices of elements
    attacker: int
        Selected attacker
    target: int
        Selected target
        
    return    
    ------
    topology: array LxL
        Segregated groups
        
    """
    topology =  np.copy(grid)
    for key in range(len(indices)):
        row, col = indices[key]
        if loyalty_mtx[key][attacker] != loyalty_mtx[key][target] :    # Not equal loyalty
            if (loyalty_mtx[key][attacker] > loyalty_mtx[key][target]) :
                topology[row][col] = 2               # attackers
            elif (loyalty_mtx[key][attacker] < loyalty_mtx[key][target]):
                topology[row][col] = 3               # defenders
            
    row_a, col_a = indices[attacker]
    row_i, col_i = indices[target]
    
    topology[row_a][col_a] = 2
    topology[row_i][col_i] = 3
    
    return topology

def vicinal(agent, L, positions_arr): 
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
    row, col = positions_arr[agent]
    
    # Check and append neighboring positions if they exist in the dictionary
    for indexes in [((row - 1) % L, col), ((row + 1) % L, col),
                 (row, (col - 1) % L), (row, (col + 1) % L)]:
        if np.any(np.all(positions_arr == indexes, axis=1)):
            index = np.where(np.all(positions_arr == indexes, axis=1))[0][0]
            neighbors.append(index)

    return neighbors

def filt(value, positions, neighbors, coalition_list, topography): 
    """
    This function filters the neighbors
    
    Parameters
    -----------
    topography: array
        Matrix of the topographical space  
    positions: array
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

def group(agent, topography, value, positions):
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
    neighbors = vicinal(agent, L, positions)
    coalition_list = [agent]
    #filterd neighbors and updated coalitions
    neighbors_filtered, coalition_list = filt(value, positions, neighbors, coalition_list, topography) 
    while len(neighbors_filtered)>0:
        temporal_neighbors = []
        #print(neighbors_filtered)
        for j in neighbors_filtered:
            local_neighbors = vicinal(j, L, positions)
         #   print(j, local_neighbors)
            new_local, coalition_list = filt(value, positions, local_neighbors, coalition_list, topography)
          #  print(new_local, coalition_list)
            if len(new_local)>0:
                temporal_neighbors = temporal_neighbors + new_local
        neighbors_filtered = temporal_neighbors
    return np.array(coalition_list, dtype=int)
    
def finding_candidates(attacker,capital, actors_pos, loyalty_mtx, grid): 
    """
    This function filters the neighbors
    
    Parameters
    -----------
    M: array (NxN)
        Commitment matrix
    grid: array (LxL)
        Grid with connected valid elements (0's) and inaccessible elements (1's). 
    actors_pos: dictionary
        actors indexes 
    attaacker: int
        Demander
    key: int
        Main target
    capital: dictionary
        actor's resources

        
    return    
    ------
    attacker_alley: list
        attacker coalition
    target_alley: list
        demander helpers
    """
    current_status = [-1, -1, -1, 0]       # target, target_resources, attacker_resources, susceptibility
    L = grid.shape[0]
    target_alley_new = np.array([-1])
    attacker_alley_new = np.array([-1])
    for key in range(len(actors_pos)):
        if key == attacker:
            continue
        topology = groups_matrix(attacker, key, actors_pos, loyalty_mtx, grid)
        #The attacker can only attack if there is a path to connnect with a target  
        attacker_alley = group(attacker, topology, 2, actors_pos)
        target_neighbors = vicinal(key, L, actors_pos)

        if not any(item in target_neighbors for item in attacker_alley):
            continue

        target_alley = group(key, topology, 3, actors_pos)

        w_att = group_resources(attacker, attacker_alley, capital, loyalty_mtx)
        w_def = group_resources(key, target_alley, capital, loyalty_mtx)

        susceptibility_target = vulnerability(w_att, w_def) * min(250,capital[key])
        if susceptibility_target > current_status[-1]:
            current_status = [key, w_def, w_att, susceptibility_target]
            attacker_alley_new = attacker_alley
            target_alley_new = target_alley
    return current_status, attacker_alley_new, target_alley_new
