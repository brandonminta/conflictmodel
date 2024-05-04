import numpy as np
import matplotlib.pyplot as plt
import time
import math
import multiprocessing as mp
from multiprocessing import Pool
from tqdm import tqdm

global c, N
c = 0.1
N = 100

#Functions
def create_arrays(N):
    matrix = np.random.uniform(300, 501, size=N)
    commitment = np.eye(N)
    return matrix, commitment
    
def target_selection(attacker, distance, L=10):
    x0, y0 = attacker // L, attacker % L  # Convert site number to (x, y) coordinates
    points = []

    # Iterate over possible x-coordinates within the range of the grid
    for dx in range(-distance, distance + 1):
        # Calculate the corresponding y-coordinate based on the Manhattan distance
        dy = distance - abs(dx)
        for dy_sign in [-1, 1]:
            y = (y0 + dy * dy_sign) % L  # Apply periodic boundary conditions
            x = (x0 + dx) % L  # Apply periodic boundary conditions
            points.append(x * L + y)  # Convert (x, y) coordinates back to site number
    arr = np.unique(points)
    arr_filtered = np.delete(arr, np.where(arr == attacker))
    return np.random.choice(arr_filtered)


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


def lattice_evolution(capital, commitment, probabilities, k, q, c=0.1):
    #Random ambitious leader
    active = np.random.randint(0, N)
    #Target selection with probability 
    order = np.random.choice(orders_arr, p=probabilities)
    target = target_selection(active, order, 10)

    #Coalition formation
    target_alley = np.where(commitment[:, active] < commitment[target, :])[0]
    target_alley = np.union1d(target_alley, [target])
    attacker_alley = np.where(commitment[:, active] > commitment[target, :])[0]
    attacker_alley = np.union1d(attacker_alley, [active])
    
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
    i, j, k, r, q, decay_rate = args
    gini_values = 0
    activity_values = 0
    gini_values_std = 0
    activity_values_std = 0
    probabilities = np.exp(-decay_rate * orders_arr)
    probabilities /= np.sum(probabilities)

    for _ in range(simulations):
        capital, commitment = create_arrays(N)
        gini_indx = np.zeros((steps))
        activity = np.zeros((steps))

        for m in range(steps+transient):
            capital, commitment, activity_val = lattice_evolution(capital, commitment, probabilities, k, q)
            if m >= transient:
                gini_indx[m - transient] = gini(capital)
                activity[m - transient] = activity_val
            if m % (N // 3) == 0:
                capital += r

        gini_values += np.mean(gini_indx)
        activity_values += np.mean(activity)
        gini_values_std += np.std(gini_indx)
        activity_values_std += np.std(activity)

    average_gini = gini_values / simulations
    average_activity = activity_values / simulations
    var_gini = gini_values_std / simulations
    var_activity = activity_values_std / simulations

    return i, j, average_gini, average_activity, var_gini, var_activity


L = int(np.sqrt(N))
simulations = 20                                   ####### change ######
transient = 100000
steps = 100000
orders_arr = np.arange(1, 10 + 1)

k = 0.25                                           ####### uncomment ######
k_num = 25
k_values = np.linspace(0, 1, k_num, endpoint=True)

r = 20                                            ####### change ######
r_num = 25
r_values = np.linspace(0, 200, r_num, endpoint=True)

q = 250                                           ####### change ######
q_num = 25
q_values = np.linspace(1, 750, q_num, endpoint=True)

decay_rate = 0.7                               ####### change ######
decay_num = 25
d_values = np.linspace(0, 1, decay_num, endpoint=True)

gini_matrix_av = np.zeros((k_num, r_num))  ####### change ######
act_matrix_av = np.zeros((k_num, r_num))   ####### change ######

gini_matrix_sd = np.zeros((k_num, r_num))   ####### change ######
act_matrix_sd = np.zeros((k_num, r_num))   ####### change ######



if __name__ == "__main__":
    start_time = time.time()

    # Use Pool to parallelize the outer loop
    n_cpu = mp.cpu_count()
    pool = Pool(processes=n_cpu)
    args_list_kr = [(i, j, k_values[i], r_values[j], q, decay_rate) for i in range(k_num) for j in range(r_num)]
    #args_list_rq = [(i, j, k, r_values[i], q_values[j], decay_rate) for i in range(r_num) for j in range(q_num)]  
    #args_list_rd = [(i, j, k, r_values[i], q, d_values[j]) for i in range(r_num) for j in range(decay_num)]  
    #args_list_kd = [(i, j, k_values[i], r, q, d_values[j]) for i in range(k_num) for j in range(decay_num)]  
    results_kr = list(tqdm(pool.imap(parallel_simulation, args_list_kr), total=len(args_list_kr)))
    
    pool.close()
    pool.join()

    # Update matrices based on results
    for result in results_kr:
        i, j, average_gini, average_activity, var_gini, var_activity = result
        gini_matrix_av[i, j] = average_gini
        act_matrix_av[i, j] = average_activity
        gini_matrix_sd[i, j] = var_gini
        act_matrix_sd[i, j] = var_activity
        
    filename1 = 'GiniMean-Grid25x25-kr-Sim20-N100.npy'    ####### change ######
    filename2 = 'ActMean-Grid25x25-kr-Sim20-N100.npy'     ####### change ######
    filename3 = 'GiniStd-Grid25x25-kr-Sim20-N100.npy'    ####### change ######
    filename4 = 'ActStd-Grid25x25-kr-Sim20-N100.npy'      ####### change ######
    # Save the array to the specified filename
    np.save(filename1, gini_matrix_av)
    np.save(filename2, act_matrix_av)
    np.save(filename3, gini_matrix_sd)
    np.save(filename4, act_matrix_sd)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Simulation completed in {elapsed_time:.4f} seconds.")
    
