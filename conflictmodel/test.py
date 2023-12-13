import numpy as np
import matplotlib.pyplot as plt
import time
import math
import multiprocessing as mp


global c
c = 0.1

#Functions
def create_arrays(N):
    matrix = np.random.uniform(300, 501, size=N)
    commitment = np.eye(N)
    return matrix, commitment
    
def von_neumann_neighbors(site, lattice_size, order=1):
    row, col = site
    indices = np.arange(-order, order + 1)
    neighbors = np.array(np.meshgrid(indices, indices)).T.reshape(-1, 2)
    neighbors = neighbors[np.abs(neighbors).sum(axis=1) == order]
    neighbors = (neighbors + site) % lattice_size

    return neighbors

def vulnerability(r_i, r_j):
    return (r_i - r_j) / r_i if r_i > 0 and r_j > 0 else 0

def response(capital, commitment, target, active, w_att, w_def, target_alley, attacker_alley, k, q, c=0.1):  
    damage_by_attacker = min(k*w_att, w_def)                                   #can't cause more damage 
    damage_by_defender = min(k*w_def, w_att)
    activity_val = 0
    if min(q,capital[target]) > (damage_by_attacker*capital[target]/w_def):    #Consider only target's damage
        for i in target_alley:
            offering = commitment[i][target] * capital[i]
            contribution_loss = damage_by_attacker*offering/w_def
            capital[i] -= contribution_loss
            for m in target_alley:
                if 1. - c >= commitment[i][m] >= 0.:
                    commitment[i][m] = commitment[i][m] + c
                    activity_val += 1
            for n in attacker_alley:
                if 1. >= commitment[i][n] >= c:
                    commitment[i][n] = commitment[i][n] - c
                    commitment[n][i] = commitment[n][i] - c
                    activity_val += 2

        for j in attacker_alley:
            offering = commitment[j][active] * capital[j]
            contribution_loss = damage_by_defender*offering/w_att
            capital[j] -= contribution_loss
            for l in attacker_alley:
                if 1. - c >= commitment[j][l] >= 0.:
                    commitment[j][l] = commitment[j][l] + c
                    activity_val += 1
    else:
        money = min(q,capital[target])
        capital[target], capital[active] = capital[target] - money, capital[active] + money
        if 1. - c >= commitment[target][active] >= 0. :
            commitment[target][active] += c
            commitment[active][target] += c
            activity_val += 2
            
    return capital, commitment, activity_val 


def lattice_evolution(capital, commitment, decay_rate, k, q, c=0.1):
    #Random ambitious leader
    active = np.random.randint(0, N)
    #Target selection with probability 
    order = np.random.choice(orders_arr, p=probabilities)
    neighbors_target = von_neumann_neighbors([active//L,active%L], L, order=order)
    target_idx =  neighbors_target[np.random.choice(neighbors_target.shape[0])]
    target = L*target_idx[0]+target_idx[1]

    #Coalition formation
    target_alley = np.where(commitment[:, active] < commitment[target, :])[0]
    attacker_alley = np.where(commitment[:, active] > commitment[target, :])[0]
    if not target_alley.size>0:
        target_alley = np.array([target])
    if not attacker_alley.size>0:
        attacker_alley = np.array([active]) 
    #Resources
    w_att = np.sum(capital[attacker_alley]*commitment[attacker_alley,active])
    w_def = np.sum(capital[target_alley]*commitment[target_alley,target])
    #Susceptibility metric
    susceptibility = vulnerability(w_att, w_def) * min(q,capital[target])
    #Demand
    if susceptibility <= 0:
        return capital, commitment, 0
    capital, commitment, activity_val = response(capital, commitment, target, active,w_att,w_def, target_alley, attacker_alley, k,q, c)
    return capital, commitment, activity_val

def gini(arr):
    N = len(arr)
    numerator = np.sum(np.abs(np.subtract.outer(arr, arr)))
    denominator = 2 * N * np.sum(arr)
    index = numerator / denominator
    return index

def parallel_simulation(args):
    i, j, k, r = args
    gini_values = 0
    activity_values = 0

    for _ in range(simulations):
        capital, commitment = create_arrays(N)
        gini_indx = np.zeros((steps - transient))
        activity = np.zeros((steps - transient))

        for m in range(steps):
            capital, commitment, activity_val = lattice_evolution(capital, commitment, decay_rate, k, q)
            if m >= transient:
                gini_indx[m - transient] = gini(capital)
                activity[m - transient] = activity_val
            if m % (N // 3) == 0:
                capital += r

        gini_values += np.mean(gini_indx)
        activity_values += np.mean(activity)

    average_gini = gini_values / simulations
    average_activity = activity_values / simulations

    return i, j, average_gini, average_activity


N = 100
L = int(np.sqrt(N))
simulations = 20
transient = 110000
steps = 120000
decay_rate = 0.1  # Make sure to define decay_rate

orders_arr = np.arange(1, L // 2 + 1)
probabilities = np.exp(-decay_rate * orders_arr)
probabilities /= np.sum(probabilities)

k_num = 25
r_num = 25
k_values = np.linspace(0, 1, k_num, endpoint=True)
r_values = np.linspace(0, 400, r_num, endpoint=True)
gini_matrix = np.zeros((k_num, r_num))
act_matrix = np.zeros((k_num, r_num))

if __name__ == "__main__":
    start_time = time.time()

    # Use Pool to parallelize the outer loop
    n_cpu = mp.cpu_count()
    pool = Pool(processes=n_cpu)
    args_list = [(i, j, k_values[i], r_values[j]) for i in range(k_num) for j in range(r_num)]
    results = list(tqdm(pool.imap(parallel_simulation, args_list), total=len(args_list)))
    
    pool.close()
    pool.join()

    # Update matrices based on results
    for result in results:
        i, j, average_gini, average_activity = result
        gini_matrix[i, j] = average_gini
        act_matrix[i, j] = average_activity
        
    filename1 = 'Gini-Grid25x25-Sim20-N100.npy'
    filename2 = 'Act-Grid25x25-Sim20-N100.npy'
    # Save the array to the specified filename
    np.save(filename1, gini_matrix)
    np.save(filename2, act_matrix)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Simulation completed in {elapsed_time:.4f} seconds.")
    