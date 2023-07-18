# -------------------------------------------------------------------------
# Created on January - 2023
# @author: Brandon Minta
# -------------------------------------------------------------------------
#
#   A simulation of a socio-dynamical model of conflicts based on 
#   "pay or else" for setting simple interaction  rules among agents.  
#
# --------------------------------------------------------------------------


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
import path_search as pth


#Constants 
c = 0.1                         # commitment fluctuation 
l_wealth, u_wealth = 300, 500   # Wealth limits distribution
k = 0.25                        # Tribute 
r = 20                          # Resource injection
q = 250                         #Tribute

#Functions
def resources(actor, coalition, matrix, capital):    
    """
    The function calculates the total resources of a coalition
    
    Parameters
    -----------
    actor: int
        The actor of interest 
    coalition: array
        Coalition of the actor
    matrix: array
        Commitment matrix
    capital: dictionary
        Capital of each member
    
    return
    ------
    money: float
        The total resources of the coalition
        
    """
    money = 0
    for i in coalition:
        money += matrix[i][actor] * capital[i]        
    return money   

def wealth():
    # This function assigns a random value from upper and lower wealth limits
    return random.randrange(l_wealth,u_wealth,1) 
    
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
  
def min_pay(r_j):   
    # This function computes the minimum ability to pay of agent j
    return q if r_j > q else r_j 

def activity(M0, M1):
    """
    This function computes the change in the commitment matrix
    used by K.Kaneko.
    
    Parameters
    -----------
    M0: array matrix (nxn) 
        Commitment matrix at time t
    M1: array matrix (nxn) 
        Commitment matrix at time t+1 
        
    return    
    ------
    float
        The activity of the commitment matrix
        
    """
    num = M0.shape[0]
    new = np.abs(M1 - M0)
    return np.sum(new) / ((num - 1) ** 2)

def groups_matrix(matrix, grid, indices, attacker,  target) :
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
    for j in range(len(indices)):
        row, col = indices[j]
        if matrix[j][attacker] > matrix[j][target]:
            topology[row][col] = 2               #attackers
        elif matrix[j][attacker] < matrix[j][target]:
            topology[row][col] = 3               #defenders
            
    row_a, col_a = indices[attacker]
    row_i, col_i = indices[target]
    
    topology[row_a][col_a] = 2
    topology[row_i][col_i] = 3
    
    return topology

class AccessibilityGrid:
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

class DictionaryUpdater:
    increment = r
    @classmethod
    def update_dic(cls, resources_dic):
        # This method increases the values of the resources of every agent
        for key in resources_dic:
            resources_dic[key] += cls.increment

class Simulation:
    def __init__(self, L, years, density=None):
        self.L = L
        self.years = years
        self.density = density
        self.data = pd.DataFrame(columns=['Attackers', 'Defenders', 'Activity', 'Status', 'TCL', 'ACL', 'Active', 'Target'])
        self.matrix_list = []
        access_grid = AccessibilityGrid(self.L)
        self.grid = access_grid.generate_grid(self.density)
        self.actors_pos =self.get_actor_positions()
        self.N = len(self.actors_pos)
        self.capital ={i: round(wealth(), 3) for i in range(self.N)}
        self.M = np.identity(self.N) 
        self.old_M = np.copy(self.M)

    def get_actor_positions(self):
        actors_pos = {}
        for idx, (row, col) in enumerate(np.argwhere(self.grid == 0)):
            actor_index = idx
            actors_pos[actor_index] = (row, col)
        return actors_pos

    def simulate_activation(self):
        class Status:
            def __init__(self, target, target_alley, target_resources, attacker_alley, attacker_resources, susceptibility):
                self.target = target
                self.target_alley = target_alley
                self.target_resources = target_resources
                self.attacker_alley = attacker_alley
                self.attacker_resources = attacker_resources
                self.susceptibility = susceptibility

        def candidates_2d(attacker, capital, M, grid, actors_pos):
            current_status = Status(np.nan, np.nan, np.nan, np.nan, np.nan, 0)
            L = grid.shape[0]

            for i in range(len(self.actors_pos)):
                if i == attacker:
                    continue
                topology = groups_matrix(M, grid, actors_pos, attacker, i)
                """if land_combat = False
                attacker_alley = []
                target_alley = []

                for key, (row, col) in actors_pos.items():
                    if topology[row][col] == 2:
                        attacker_alley.append(key)
                    elif topology[row][col] == 3:
                        target_alley.append(key)

                attacker_resources = resources(attacker, attacker_alley, M, capital)
                target_resources = resources(i, target_alley, M, capital)

                susceptibility_target = vulnerability(attacker_resources, target_resources) * min_pay(capital[i])
                if susceptibility_target > current_status.susceptibility:
                    current_status = Status(i, target_alley, target_resources, attacker_alley, attacker_resources,
                                            susceptibility_target)
                 """                           
                #else if land_combat = True
                attacker_alley = pth.group(attacker, topology, 2, actors_pos)
                target_neighbors = pth.vicinal(actors_pos,i, L)

                if not any(item in target_neighbors for item in attacker_alley):
                    continue

                target_alley = pth.group(i, topology, 3, actors_pos)

                attacker_resources = resources(attacker, attacker_alley, M, capital)
                target_resources = resources(i, target_alley, M, capital)

                susceptibility_target = vulnerability(attacker_resources, target_resources) * min_pay(capital[i])
                if susceptibility_target > current_status.susceptibility:
                    current_status = Status(i, target_alley, target_resources, attacker_alley, attacker_resources,
                                            susceptibility_target)
            return current_status

        def response(attacker, state, capital, M):
            target = state.target
            target_alley = state.target_alley
            target_resources = state.target_resources
            attacker_alley = state.attacker_alley
            attacker_resources = state.attacker_resources
            target_loss, active_loss = 0, 0

            if min_pay(capital[target]) > loss(attacker_resources, target_resources, capital[target]):
                for i in target_alley:
                    offering = M[i][target] * capital[i]
                    contribution_loss = loss(attacker_resources, target_resources, offering)
                    capital[i] -= contribution_loss
                    target_loss += contribution_loss
                    for m in target_alley:
                        if 1 - c >= M[i][m] >= 0:
                            M[i][m] = M[i][m] + c
                    for n in attacker_alley:
                        if 1 >= M[i][n] >= c:
                            M[i][n] = M[i][n] - c
                            M[n][i] = M[n][i] - c

                for j in attacker_alley:
                    offering = M[j][attacker] * capital[j]
                    contribution_loss = loss(target_resources, attacker_resources, offering)
                    capital[j] -= contribution_loss
                    active_loss += contribution_loss
                    for l in attacker_alley:
                        if 1 - c >= M[j][l] >= 0:
                            M[j][l] = M[j][l] + c

                return 1, target_loss, active_loss

            else:
                money = min_pay(capital[target])
                capital[target], capital[attacker] = capital[target] - money, capital[attacker] + money
                if 1 - c >= M[target][attacker] >= 0:
                    M[target][attacker] += c
                    M[attacker][target] += c

                return 0, money, 0


        attacker = random.randrange(0, len(self.actors_pos))
        prospect = candidates_2d(attacker, self.capital, self.M, self.grid, self.actors_pos)
        tau, alpha = prospect.target_alley, prospect.attacker_alley
        if prospect.susceptibility > 0:
            decision, t_loss, a_loss = response(attacker, prospect, self.capital, self.M)
            return decision, tau, alpha, t_loss, a_loss, attacker, prospect.target
        else:
            return 0, tau, alpha, 0, 0, attacker, prospect.target


    def run_simulation(self):
        for _ in range(self.years):
            for _ in range(self.N // 3):
                status, tau, alpha, t_loss, a_loss, active, target = self.simulate_activation()
                resources = {k: round(v, 3) for k, v in self.capital.items()}
                iteration_data = {
                    **resources,
                    'Attackers': alpha,
                    'Defenders': tau,
                    'Activity': activity(self.old_M, self.M),
                    'Status': status,
                    'TCL': round(t_loss, 3),
                    'ACL': round(a_loss, 3),
                    'Active': active,
                    'Target': target
                }
                self.data = pd.concat([self.data, pd.DataFrame([iteration_data])], ignore_index=True)

                self.old_M = np.round(np.copy(self.M), 1)

            DictionaryUpdater.update_dic(self.capital)

        return self.data, self.grid, self.M
