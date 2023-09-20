#!/usr/bin/env python


# Created on January - 2023
# @author: Brandon Minta
# -------------------------------------------------------------------------
#
#   A simulation of a socio-dynamical model of conflicts based on 
#   "pay or else" for setting simple interaction  rules among agents.  
#
# --------------------------------------------------------------------------

import argparse
import csv
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import h5py
from tqdm import tqdm
import multiprocessing as mp

import random
import candidate_finder as finder

#seed_value = 42
#random.seed(seed_value)
#np.random.seed(seed_value)

def parse_arguments():
    parser = argparse.ArgumentParser(description='Run the simulation script.')
    parser.add_argument('-L', type=int, default=5, help='Value for L')
    parser.add_argument('--years', type=int, default=1000, help='Number of iterations')
    parser.add_argument('--density', type=float, default=0.0, help='Density value')
    parser.add_argument("--r", type=int, default=20, help="Value of productivity/income")
    parser.add_argument("--q", type=int, default=250, help="Value of tribute")
    parser.add_argument('--ncpu', type=int, default=4, help='Number of processes')
    parser.add_argument("--output_dir", required=True, help="Output directory")
    return parser.parse_args()
#Constants     
args =  parse_arguments()
c = 0.1                         # commitment fluctuation 
k = 0.25                        # Destructiveness 
demands = 3                     # Demands per Year Cycle (1/3 of N)
period_steps = 25               # Period for data collection
r = args.r                      # Production
q = args.q                      # Tribute

#Functions
def saving_data(rank, output_dir, loyalty_list, tribute_list, wealth_matrix, info_matrix):
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    subdirectory_name = f"run_{rank}"  # Create the subdirectory name
    subdirectory_path = os.path.join(output_dir, subdirectory_name)
    # Create the subdirectory if it doesn't exist (with exist_ok=True, it won't raise an error if it already exists)
    os.makedirs(subdirectory_path, exist_ok=True)

    # Save data to HDF5 files
    with h5py.File(os.path.join(subdirectory_path, f"loyalty_output_{rank}.h5"), 'w') as hf:
        for i, array in enumerate(loyalty_list):
            hf.create_dataset(f'loyalty_matrix_{i}', data=array)

    with h5py.File(os.path.join(subdirectory_path, f"tribute_output_{rank}.h5"), 'w') as hf:
        for i, array in enumerate(tribute_list):
            hf.create_dataset(f'tribute_matrix_{i}', data=array)

    with h5py.File(os.path.join(subdirectory_path, f"wealth_output_{rank}.h5"), 'w') as hf:
        hf.create_dataset('wealth', data=wealth_matrix)

    with h5py.File(os.path.join(subdirectory_path, f"states_output_{rank}.h5"), 'w') as hf:
        hf.create_dataset('info', data=info_matrix)

    # Save data to CSV file
    #if groups_list != None:
    #    with open(os.path.join(subdirectory_path, f"groups_output_${rank}.csv"), 'w', newline='') as csvfile:
    #        fieldnames = ['iteration', 'att', 'def']
    #        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    #        writer.writeheader()
    #        for i, iteration_data in enumerate(groups_list):
    #            writer.writerow({'iteration': i+1, 'att': iteration_data['att'], 'def': iteration_data['def']})
   
def loss(r_i, r_j, contribution):         
    """
    The function calculates agent j's loss of resources, caused by conflict with i
    
    Parameters
    -----------
    r_i: float
        Resources of i's coalition  
    r_j: float
        Resources of j's coalition  
    contribution: float
        j's contribution of resources to the coalition
        
    return
    -------
    float:
        Resources loss of j member of it's coalition
        
    """
    return k*r_i*(contribution)/r_j

def response(attacker, state, target_alley, attacker_alley, capital, loyalty_mtx, tribute_mtx):
    target = state[0]
    w_def = state[1]
    w_att = state[2]
    
    damage_by_attacker = min(k*w_att, w_def)   #can't cause more damage 
    damage_by_defender = min(k*w_def, w_att)
    loyalty = np.copy(loyalty_mtx[attacker][target])
    if min(q,capital[target]) > (damage_by_attacker*capital[target]/w_def ):    #Consider only target's damage
        for i in target_alley:
            offering = loyalty_mtx[i][target] * capital[i]
            contribution_loss = damage_by_attacker*offering/w_def
            capital[i] -= round(contribution_loss,1)
            for m in target_alley:
                if 1. - c >= loyalty_mtx[i][m] >= 0.:
                    loyalty_mtx[i][m] = loyalty_mtx[i][m] + c
            for n in attacker_alley:
                if 1. >= loyalty_mtx[i][n] >= c:
                    loyalty_mtx[i][n] = loyalty_mtx[i][n] - c
                    loyalty_mtx[n][i] = loyalty_mtx[n][i] - c

        for j in attacker_alley:
            offering = loyalty_mtx[j][attacker] * capital[j]
            contribution_loss = damage_by_defender*offering/w_att
            capital[j] -= round(contribution_loss,1) 
            for l in attacker_alley:
                if 1. - c >= loyalty_mtx[j][l] >= 0.:
                    loyalty_mtx[j][l] = loyalty_mtx[j][l] + c
        return 1, loyalty  # Conflict: True ; Loyalty between attacker and target
    else:
        money = round(min(q,capital[target]),1)
        capital[target], capital[attacker] = capital[target] - money, capital[attacker] + money
        tribute_mtx[target][attacker] += 1
        if 1. - c >= loyalty_mtx[target][attacker] >= 0. :
            loyalty_mtx[target][attacker] += c
            loyalty_mtx[attacker][target] += c
        return 0, loyalty   # Conflict: False ; Loyalty between attacker and target

class GridGenerator:
    def __init__(self, L):
        """
        Initialize the AccessibilityGrid class.

        Parameters
        ----------
        L: int
            Size of the topological grid.
        """
        self.L = L
        self.grid = np.zeros((L, L), dtype=int)

    def generate_grid(self, density=None):
        """
        Generate an array with a connected grid of valid elements (0's)
        surrounded by inaccessible elements (1's) using a specified density of obstacles.

        Parameters
        ----------
        density: float (optional)
            Density of obstacles on the topological grid.

        Returns
        -------
        topology: array LxL
            Grid with connected valid elements (0's) and inaccessible elements (1's).
        """
        if density is None:
            return self.grid
        
        self.grid = np.ones((self.L, self.L), dtype=int)
        num_obstacles = int(density * self.L * self.L)
        valid_elements = self.L * self.L - num_obstacles

        # Select a random cell as the origin for expanding valid elements
        origin = np.random.randint(self.L), np.random.randint(self.L)

        # Expand valid elements from the origin
        self._expand_valid_elements(origin, valid_elements)

        return self.grid

    def _expand_valid_elements(self, origin, valid_elements):
        """
        Expand valid elements (0's) from the given origin in the grid until reaching the specified count
        using a Monte Carlo method.

        Parameters
        ----------
        origin: tuple
            Coordinates (row, col) of the origin cell.
        valid_elements: int
            Number of valid elements (0's) to reach.

        Returns
        -------
        int
            Number of valid elements added.
        """
        count = 0

        while count < valid_elements:
            # Randomly select a neighbor cell
            neighbors = [
                [(origin[0] - 1) % self.L, origin[1]],    # Up
                [(origin[0] + 1) % self.L, origin[1]],    # Down
                [origin[0], (origin[1] - 1) % self.L],    # Left
                [origin[0], (origin[1] + 1) % self.L]     # Right
            ]
            neighbor_row, neighbor_col = neighbors[np.random.randint(4)]

            if self.grid[neighbor_row, neighbor_col] == 1:
                # Add the neighbor cell as a valid element
                self.grid[neighbor_row, neighbor_col] = 0
                count += 1

            # Update the origin to the neighbor cell
            origin = (neighbor_row, neighbor_col)

class Simulation:
    def __init__(self, L, years, density=None):
        self.L = L
        self.years = years                                          # Years to be run 
        self.density = density                                      # Density of obstacles
        access_grid = GridGenerator(self.L)
        self.grid = access_grid.generate_grid(self.density)         # Actors Distribution in the grid
        self.actors_pos = np.argwhere(self.grid == 0)               # Positions of actors in the grid
        self.N = len(self.actors_pos)
        self.capital = np.zeros(self.N + 2)                         # Wealth + W_a + W_d
        self.loyalty_mtx = np.identity(self.N, dtype=np.float64) 
        self.tribute_mtx = np.zeros((self.N, self.N), dtype=int)
        
    def simulate_activation(self):
        attacker = random.randrange(0, self.N)
        current_status, attacker_alley,target_alley  = finder.finding_candidates(attacker, self.capital, self.actors_pos, self.loyalty_mtx, self.grid)
        if current_status[0] != -1:    #Target
            decision, loyalty = response(attacker, current_status, target_alley, attacker_alley, self.capital, self.loyalty_mtx, self.tribute_mtx)
            return decision, loyalty, target_alley, attacker_alley, attacker, current_status #Decision,loyalty,T alley, A alley, attecker, (target, wt,wa)
        else:         
            return -1, -1, target_alley, attacker_alley, attacker, current_status,      #1 conflict, 0 tribute, -1 not capable 

    def run_simulation(self, pbar = None):  
        #Get the total number of iterations
        total_iterations = self.years * (self.N // demands)
        #Initialize the resources for each actor
        for i in range(self.N):
            self.capital[i] = round(random.randrange(300,500,1), 1)
      
        loyalty_list, tribute_list = [], []  # Lists for periodic loyalty, tribute matrices, and groups formation
        wealth_matrix = np.zeros((total_iterations+1, self.N+2), dtype=np.float32)  # N + 2 columns: actors wealths and W_d, W_a  
        wealth_matrix[0] = self.capital
        info_matrix = np.zeros((total_iterations+1, 6), dtype=int)  # 5 columns: decision, loyalty*10, att, def, size_coal_a, size_coal_t
        #Create the charging bar if pbar is provided
        if pbar:
            pbar.reset(total=total_iterations)
        # Loop for each simulation run 
        iterator, period = 1, 0
        for year in range(self.years):
            if year == period:
                loyalty_list.append(np.copy(self.loyalty_mtx))
                tribute_list.append(np.copy(self.tribute_mtx))
                period += period_steps
                
            for k in range( self.N // demands):
                decision, loyalty, target_alley, attacker_alley, attacker, current_status = self.simulate_activation()
                selected_target, self.capital[-2:] = current_status[0], current_status[1:3]  # Assigning values from current_status
                info_matrix[iterator] = np.array([decision, int(loyalty*10), attacker, selected_target, len(target_alley), len(attacker_alley)])            
                #iteration_data = {'att': attacker_alley.tolist(), 'def': target_alley.tolist()}
                #groups_list.append(iteration_data)
                wealth_matrix[iterator] = self.capital
                iterator += 1
                if pbar:
                    pbar.update(1)
            self.capital[:-2] += r      
        return loyalty_list, tribute_list, wealth_matrix, info_matrix

def run(rank, output_dir, L, years, density):
    simulation = Simulation(L, years, density)
    with tqdm(total=total_iterations, desc="Simulation Progress", unit="iteration") as pbar:
        loyalty_list, tribute_list, wealth_matrix, info_matrix = simulation.run_simulation(pbar=pbar) 
        saving_data(rank, output_dir, loyalty_list, tribute_list, wealth_matrix, info_matrix)
    
if __name__ == "__main__":
    args = parse_arguments()
    L = args.L
    years = args.years
    density = args.density
    output_dir = args.output_dir
    ncpu = args.ncpu
    
    N = L * L - int(density * L * L)
    total_iterations = years * (N // demands)
    
    # Get number of laptop CPUs
    n_cpu = mp.cpu_count()
    # Call Pool
    pool = mp.Pool(processes=ncpu)
    # Create a list of tuples containing all combinations of XM and B values
    parameters = [(rank, output_dir,L, years, density) for rank in range(1,ncpu+1)]
    # Call run for all parameter tuples using pool.map
    pool.starmap(run, parameters)
    # Close the pool
    pool.close()
