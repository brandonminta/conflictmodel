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
import os

import numpy as np
import math
import h5py
from tqdm import tqdm
import multiprocessing as mp
import networkx as nx

import random
# --------------------------------------------------------------------------
def parse_arguments():
    parser = argparse.ArgumentParser(description='Run the simulation script.')
    parser.add_argument('-N', type=int, default=100, help='Number of actors')
    parser.add_argument('--years', type=int, default=1000, help='Number of iterations')
    parser.add_argument("--network_type", choices=["watts_strogatz", "grid_2d", "complete", "wheel"], default="grid_2d", help="Type of network")
    parser.add_argument('--ncpu', type=int, default=1, help='Number of processes')
    parser.add_argument("--output_dir", required=True, help="Output directory")
    return parser.parse_args()
#Constants     
args =  parse_arguments()
c = 1                           # commitment fluctuation 
k = 0.25                        # Destructiveness 
demands = 3                     # Demands per Year Cycle (1/3 of N)
period_steps = 1               # Period for data collection
r = 20                          # harvest
q = 250                         # Tribute

# Functions --------------------------------------------------------------------------
def saving_data(rank, output_dir, loyalty_list, tribute_list, data_matrix):
    subdirectory_name = f"run_{rank}"  # Create the subdirectory name
    subdirectory_path = os.path.join(output_dir, subdirectory_name)
    # Create the subdirectory if it doesn't exist (with exist_ok=True, it won't raise an error if it already exists)
    os.makedirs(subdirectory_path, exist_ok=True)
    # Save data to HDF5 files
    with h5py.File(os.path.join(subdirectory_path, f"loyalty_output_{rank}.h5"), 'w') as hf:
        for i, array in enumerate(loyalty_list):
            hf.create_dataset(f'loyalty_matrix_{i}', data = array)

    with h5py.File(os.path.join(subdirectory_path, f"tribute_output_{rank}.h5"), 'w') as hf:
        for i, array in enumerate(tribute_list):
            hf.create_dataset(f'tribute_matrix_{i}', data = array)

    with h5py.File(os.path.join(subdirectory_path, f"simulation_output_{rank}.h5"), 'w') as hf:
        hf.create_dataset('simulation_data', data = data_matrix)

def create_graph(num_nodes, graph_type='grid', **kwargs):
    """
    Create and return a NetworkX graph based on the specified number of nodes and graph type.
    
    Parameters:
    num_nodes (int): Number of nodes in the graph.
    graph_type (str): Type of graph ('grid', 'scale_free', or 'watts_strogatz').
    **kwargs: Additional keyword arguments specific to the chosen graph type.
    
    Returns:
    nx.Graph: NetworkX graph based on the specified parameters.
    """
    if graph_type == 'grid_2d':
        factors = []
        for i in range(1, int(math.sqrt(num_nodes)) + 1):
            if num_nodes % i == 0:
                factors.append((i, num_nodes // i))
        m, n = min(factors, key=lambda x: abs(x[0] - x[1]))
        G = nx.grid_2d_graph(m, n, periodic=True)
        mapping = {(i, j): i * n + j for i in range(m) for j in range(n)}
        G = nx.relabel_nodes(G, mapping)
    elif graph_type == 'watts_strogatz':
        k = kwargs.get('k', 4)  # Default value for k is set to 4
        p = kwargs.get('p', 0.1)  # Default value for p is set to 0.1
        G = nx.watts_strogatz_graph(num_nodes, k, p)
    elif graph_type == 'complete':
        G = nx.complete_graph(num_nodes)
    elif graph_type == 'wheel':
        G = nx.wheel_graph(num_nodes)
    return G
     
def loss(coAWealth, coBWealth, contribution):         
    """
    The function calculates loss of resources from an actor of coalition B, caused by conflict with coalition A
    """
    return k*coAWealth*(contribution)/coBWealth

def vulnerability(r_i, r_j):
    """
    The function computes the vulnerability of agent j with respect to agent i
    """
    return (r_i - r_j) / r_i if r_i > 0 and r_j > 0 else 0

def group_resources(agent, coalition_arr, capital_arr, loyalty_mtx):                                                                                   
    """
    The function calculates the total resources of a coalition
    
    Parameters
    -----------
    agent: int
        The group leader
    coalition_arr: array
        Coalition 
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

def candidate_connection(G,loyalty_mtx, attacker, target):
    """
    This function checks is there is a path between Attacker and target
    """
    # Create a subgraph to remove nodes
    H = G.copy()
    nodes_to_isolate = []
    for node in range(len(loyalty_mtx)):
        conditionA = loyalty_mtx[node, attacker]
        conditionB = loyalty_mtx[node, target] 
        if (conditionA <= conditionB ) and node not in [attacker, target]:
        #if (conditionA <= conditionB or  conditionA < 5) and node not in [attacker, target]:  #
            nodes_to_isolate.append(node)
    H.remove_nodes_from(nodes_to_isolate)
    #Check if there is a path between nodes
    if nx.has_path(H, attacker, target):
        path_length = nx.shortest_path_length(H, source=attacker, target=target)
        return path_length
    else: 
        return None

def group(G, a, b, loyalty_mtx):
    """
    This function find the spatial connected coalition of 'a'
    """
    # Create an initial grouping where each node is its own group
    group_a = {node for node in G.nodes() if (loyalty_mtx[node, a] > loyalty_mtx[node, b] and node != b) or (node == a)}
    #group_a = {node for node in G.nodes() if (loyalty_mtx[node, a] > loyalty_mtx[node, b] and node != b and loyalty_mtx[node, a] >=0.5) or (node == a)}
    subgraph = G.subgraph(group_a) 
    coalition = list(nx.descendants(subgraph, a))+ [a]
    return np.array(coalition, dtype=int)  

def candidate_exposure(G, attacker, target, capital, loyalty_mtx, path_length):
    """
    This function calculates the susceptibility 
    """
    attacker_alley = group(G, attacker, target , loyalty_mtx)
    target_alley = group(G, target, attacker, loyalty_mtx)
    w_att = group_resources(attacker, attacker_alley, capital, loyalty_mtx)
    w_def = group_resources(target, target_alley, capital, loyalty_mtx)
    #susceptibility = vulnerability(w_att, w_def) * min(q,capital[target])/(path_length**(1/4))      #Variation Path length penalization
    susceptibility = vulnerability(w_att, w_def) * min(q,capital[target])
    return susceptibility, attacker_alley, target_alley, w_att, w_def

def candidate_selection(G, loyalty_mtx, attacker, capital):
    """
    This function selects the best target possible if any
    """
    neighbors_attacker = list(G.neighbors(attacker))
    
    def get_optimal_params(target):
        """
        Helper function to calculate susceptibility and other parameters for a target node.
        """
        path_length = candidate_connection(G, loyalty_mtx, attacker, target)
        # Exclude not profitable targets and not spatially connected. 
        if capital[target] < 0.1 or path_length is None: 
        #if capital[target] < 0.1: 
            return {
            "susceptibility": 0,
            "attacker_alley": None,
            "target_alley": None,
            "w_att": None,
            "w_def": None,
            "target": None,
            "attacker":None, 
            "path_len":None
        }
        # Calculate parameters using candidate_exposure function
        susceptibility, attacker_alley, target_alley, w_att, w_def = candidate_exposure(G, attacker, target, capital, loyalty_mtx, path_length)
        return {
            "susceptibility": susceptibility,
            "attacker_alley": attacker_alley,
            "target_alley": target_alley,
            "w_att": w_att,
            "w_def": w_def,
            "target": target,
            "attacker":attacker,
            "path_len":path_length
        }
    
    # Initialize optimal parameters with the first neighbor of the attacker
    optimal = get_optimal_params(neighbors_attacker[0])   
    # Iterate through other neighbors of the attacker and update optimal parameters if a better target is found
    for target in neighbors_attacker[1:]:
        params = get_optimal_params(target)
        # Check if the target is valid and has higher susceptibility
        if params and params["susceptibility"] > optimal["susceptibility"]:
            optimal.update(params)  
            
    # Iterate through all nodes that are not neighbors of the attacker and update optimal parameters if a better target is found
    for target in [node for node in range(len(loyalty_mtx)) if node not in neighbors_attacker and node != attacker]:
        # Check if target's neighbors are not commited with attacker (no possible path)
        if all(loyalty_mtx[node][attacker] <= loyalty_mtx[node][target] for node in list(G.neighbors(target))):
            continue
        params = get_optimal_params(target)
        # Update optimal parameters if a better target is found
        if params and params["susceptibility"] > optimal["susceptibility"]:
            optimal.update(params)
    
    # Return the optimal parameters
    return optimal

def response(optimal, capital, loyalty_mtx, tribute_mtx):
    attacker = optimal["attacker"]
    target = optimal["target"]
    w_def = optimal["w_def"]
    w_att = optimal["w_att"]
    target_alley = optimal["target_alley"] 
    attacker_alley = optimal["attacker_alley"]     
    damage_by_attacker = min(k*w_att, w_def)   #can't cause more damage 
    damage_by_defender = min(k*w_def, w_att)
    loyalty = np.copy(loyalty_mtx[attacker][target])
    if min(q,capital[target]) > (damage_by_attacker*capital[target]/w_def):    #Consider only target's damage
        for i in target_alley:
            offering = 0.1*loyalty_mtx[i][target] * capital[i]
            contribution_loss = damage_by_attacker*offering/w_def
            capital[i] -= contribution_loss
            for m in target_alley:
                if 10 - c >= loyalty_mtx[i][m] >= 0:
                    loyalty_mtx[i][m] = loyalty_mtx[i][m] + c
            for n in attacker_alley:
                if 10 >= loyalty_mtx[i][n] >= c:
                    loyalty_mtx[i][n] = loyalty_mtx[i][n] - c
                    loyalty_mtx[n][i] = loyalty_mtx[n][i] - c

        for j in attacker_alley:
            offering = 0.1*loyalty_mtx[j][attacker] * capital[j]
            contribution_loss = damage_by_defender*offering/w_att
            capital[j] -= contribution_loss
            for l in attacker_alley:
                if 10 - c >= loyalty_mtx[j][l] >= 0:
                    loyalty_mtx[j][l] = loyalty_mtx[j][l] + c
        return 1, loyalty  # Conflict: True ; Loyalty between attacker and target
    else:
        money = min(q,capital[target])
        capital[target], capital[attacker] = capital[target] - money, capital[attacker] + money
        tribute_mtx[target][attacker] += 1
        if 10 - c >= loyalty_mtx[target][attacker] >= 0 :
            loyalty_mtx[target][attacker] += c
            loyalty_mtx[attacker][target] += c
        return 0, loyalty   # Conflict: False ; Loyalty between attacker and target

# Simulation --------------------------------------------------------------------------
class Simulation:
    def __init__(self, N, years, G, rank):
        self.rank = rank
        self.N = N                                              # Number of actors
        self.years = years                                      # Years 
        self.G = G                                              # Network
        self.capital = np.zeros(self.N)                         # Wealth 
        self.loyalty_mtx = np.identity(self.N, dtype= int)*10 
        self.tribute_mtx = np.zeros((self.N, self.N), dtype=int)
        
    def simulate_activation(self):
        attacker = random.randrange(0, self.N)
        optimal = candidate_selection(self.G, self.loyalty_mtx, attacker, self.capital)
        if optimal["susceptibility"] > 0:    #Target
            decision, loyalty = response(optimal, self.capital, self.loyalty_mtx, self.tribute_mtx)
            return decision, loyalty, optimal 
        else:         
            return None, None, optimal      

    def run_simulation(self, pbar = None):  
        #Get the total number of iterations
        total_iterations = self.years * (self.N // demands)
        #Initialize the resources for each actor
        for i in range(self.N):
            self.capital[i] = random.randrange(300,500,1)
        loyalty_list, tribute_list = [], []  # Lists for periodic loyalty, tribute matrices
        data_matrix = np.zeros((total_iterations+1, self.N + 10), dtype=np.float32)  
        data_matrix[0,10:] = self.capital
        #Create the charging bar if pbar is provided
        if pbar:
            pbar.reset(total=total_iterations)
        # Loop for each simulation run 
        iterator, period = 1, period_steps
        for year in range(self.years):
            if self.rank == 1:
                loyalty_list.append(np.copy(self.loyalty_mtx))
                tribute_list.append(np.copy(self.tribute_mtx))
            loyalty_prev = np.copy(self.loyalty_mtx)
            for k in range( self.N // demands):
                decision, loyalty, optimal = self.simulate_activation()
                if decision is not None:
                    # Perform element-wise comparison and count the number of differing elements
                    differences = loyalty_prev != self.loyalty_mtx
                    num_differences = np.sum(differences)
                    info = np.array([decision, optimal["attacker"], optimal["target"], loyalty, 
                            len(optimal["target_alley"] ), len(optimal["attacker_alley"]), optimal["w_def"], 
                            optimal["w_att"], optimal["path_len"],num_differences ])
                else:
                    info = np.full(10, None)
                data_matrix[iterator,:10], data_matrix[iterator,10:] = info, self.capital
                iterator += 1
                if pbar:
                    pbar.update(1)
            self.capital += r
           
            # Calculate the increment amount (5% of each element's value)    #Proportional to what they have
            #increment = 0.05 * self.capital
            #self.capital += increment
            
        loyalty_list.append(np.copy(self.loyalty_mtx))
        tribute_list.append(np.copy(self.tribute_mtx))
        return loyalty_list, tribute_list, data_matrix

def run(rank, output_dir, N, years, G):
    simulation = Simulation(N, years, G, rank)
    with tqdm(total=total_iterations, desc="Simulation Progress", unit="iteration") as pbar:
        loyalty_list, tribute_list, data_matrix = simulation.run_simulation(pbar=pbar) 
        saving_data(rank, output_dir, loyalty_list, tribute_list, data_matrix)
    
if __name__ == "__main__":
    args = parse_arguments()
    N = args.N
    years = args.years
    network_type = args.network_type
    output_dir = args.output_dir
    ncpu = args.ncpu
    
    network = create_graph(N, graph_type = network_type)
    total_iterations = years * (N // demands)
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Get number of laptop CPUs
    n_cpu = mp.cpu_count()
    # Call Pool
    pool = mp.Pool(processes=ncpu)
    # Create a list of tuples containing all combinations of XM and B values
    parameters = [(rank, output_dir, N, years, network) for rank in range(1,ncpu+1)]
    # Call run for all parameter tuples using pool.map
    pool.starmap(run, parameters)
    # Close the pool
    pool.close()