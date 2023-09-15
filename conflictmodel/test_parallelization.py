from mpi4py import MPI
import numpy as np
import time
import random

def broadcast_constants():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        start_time = MPI.Wtime()
        # Set a fixed seed value for reproducibility
        seed_value = 42  # You can choose any integer value you like

        # Seed the random number generator
        random.seed(seed_value)
        np.random.seed(seed_value)
        # Define your constants here
        attacker = 423
        L = 30
        N = L * L
        capital = np.zeros(N + 2)
        for i in range(N):
            capital[i] = round(random.randrange(300000, 500000, 1), 1)
        grid = np.zeros((L, L), dtype=int)
        actors_pos = np.argwhere(grid == 0)
        loyalty_mtx = np.eye(N)
        for i in range(N):
            for j in range(i+1, N):
                random_value = round(np.random.uniform(0.1, 1.0), 1)
                loyalty_mtx[i, j] = loyalty_mtx[j, i] = random_value
    
        # Broadcast constants to all processes
        constants = (attacker, capital, actors_pos, loyalty_mtx, grid)
        end_time = MPI.Wtime()
        print("Constants initialization time: " + str(end_time-start_time))
    else:
        constants = None

    # Broadcast the constants from rank 0 to all other processes
    attacker, capital, actors_pos, loyalty_mtx, grid = comm.bcast(constants, root=0)

    return attacker, capital, actors_pos, loyalty_mtx, grid

def group_resources(agent, coalition_arr, capital_arr, loyalty_mtx):                                                                                   
    money = 0
    for i in coalition_arr:
        money += loyalty_mtx[i][agent] * capital_arr[i]        
    return money   

def vulnerability(r_i, r_j):
    return (r_i - r_j) / r_i if r_i > 0 and r_j > 0 else 0

def groups_matrix(attacker, target, indices, loyalty_mtx, grid) :
    topology =  np.copy(grid)
    for key in range(len(indices)):
        row, col = indices[key]
        if abs(loyalty_mtx[key][attacker] - loyalty_mtx[key][target]) > 1e-9 :    # Not equal loyalty
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

    neighbors_filtered  = []
    for i in neighbors:
        row, col = positions[i]
        if (topography[row, col] == value) and (i  not in coalition_list):
            coalition_list.append(i)
            neighbors_filtered.append(i)
    
    return neighbors_filtered, coalition_list

def group(agent, topography, value, positions):

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

    
def finding_candidates(attacker, capital, actors_pos, loyalty_mtx, grid, land_combat=True):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    current_status = [-1, -1, -1, 0]  # target, target_resources, attacker_resources, susceptibility
    L = grid.shape[0]
    target_alley_new = np.array([-1])
    attacker_alley_new = np.array([-1])

    # Determine the range of keys to process for each MPI process
    keys_per_process = len(actors_pos) // size
    start_key = rank * keys_per_process
    end_key = start_key + keys_per_process
    start_time = MPI.Wtime()
    for key in range(start_key, end_key):
        if key == attacker:
            continue

        topology = groups_matrix(attacker, key, actors_pos, loyalty_mtx, grid)

        if land_combat == True:
            attacker_alley = group(attacker, topology, 2, actors_pos)
            target_neighbors = vicinal(key, L, actors_pos)

            if not any(item in target_neighbors for item in attacker_alley):
                continue

            target_alley = group(key, topology, 3, actors_pos)

        w_att = group_resources(attacker, attacker_alley, capital, loyalty_mtx)
        w_def = group_resources(key, target_alley, capital, loyalty_mtx)

        susceptibility_target = vulnerability(w_att, w_def) * min(250, capital[key])
        if susceptibility_target > current_status[-1]:
            current_status = [key, w_def, w_att, susceptibility_target]
            attacker_alley_new = attacker_alley
            target_alley_new = target_alley
    end_time = MPI.Wtime()
    if rank == 0:
        print("Finding Candidates per rank: " + str(end_time-start_time))
    # Communicate results to find the maximum susceptibility_target across all processes
    start_time = MPI.Wtime()
    max_susceptibility = comm.allreduce(current_status[-1], op=MPI.MAX)

    # Determine which rank has the maximum susceptibility_target
    rank_with_max_susceptibility = comm.allreduce(rank if current_status[-1] == max_susceptibility else -1, op=MPI.MAX)

    # Broadcast the final results to all processes
    current_status = comm.bcast(current_status, root=rank_with_max_susceptibility)
    attacker_alley_new = comm.bcast(attacker_alley_new, root=rank_with_max_susceptibility)
    target_alley_new = comm.bcast(target_alley_new, root=rank_with_max_susceptibility)
    end_time = MPI.Wtime()
    if rank == 0:
        print("Best candidate among rank - time: " + str(end_time-start_time))
        return current_status, attacker_alley_new, target_alley_new
    
if __name__ == "__main__":
    # Define your inputs here
    attacker, capital, actors_pos, loyalty_mtx, grid = broadcast_constants()
  
    #comm = MPI.COMM_WORLD
    result = finding_candidates(attacker, capital, actors_pos, loyalty_mtx, grid)
    #inter(actors_pos)
    print(result)
    
