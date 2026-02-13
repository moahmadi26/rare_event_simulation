import sys
import json
import time
import multiprocessing
from guidedwssa import guidedwssa
from prism_parser import parser
import numpy as np
from utils import suppress_c_output
import math


def main(json_path):
    start_time = time.time()
    #############################################################################################
    num_procs = 16       # number of processors used for parallel execution
    
    # Hyperparameters
    K = 4                # K ensembles of size N are used to estimate the probability of event
    N = 1_000_000       # number of trajectories used in each ensemble
    negative_method = 'C'  # method for resolving negative propensities ('A', 'B', or 'C')
    #############################################################################################
    
    with open(json_path, 'r') as f:
        json_data = json.load(f)
    
    model_path = json_data['model_path']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    t_max = float(json_data['max_time'])
    
    with suppress_c_output():
        model = parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    # Define F vector (projection vector for observable)
    # By default, F projects onto the target variable
    num_species = len(model.get_initial_state())
    F = np.zeros(num_species)
    F[target_index] = 1.0
    
    print(f"Running guided wSSA...")
    print(f"Model: {model_path}")
    print(f"Target: {target_var} = {target_value}")
    print(f"Time horizon: {t_max}")
    print(f"Negative resolution method: {negative_method}")
    print("-" * 50)
    
    start_time = time.time()
       
    # Estimation phase
    print(f"Running estimation with K={K} ensembles of N={N} trajectories each...")
    start_time = time.time()
    
    p_vector = []
    
    # Run K ensembles of size N
    for i in range(K):
        print(f"Ensemble {i+1}/{K}")
        
        N_vec = [N // num_procs 
                 if j != num_procs - 1 
                 else N - ((num_procs - 1) * (N // num_procs)) 
                 for j in range(num_procs)]
        
        tasks = [(model_path, N_vec_j, t_max, target_index, target_value, F, negative_method) 
                 for N_vec_j in N_vec]
        
        with multiprocessing.Pool(processes=num_procs) as pool:
            results = pool.starmap(guidedwssa, tasks)
        
        # Aggregate results for this ensemble
        ensemble_p = 0
        ensemble_var = 0
        for p, var, conf in results:
            ensemble_p += p * (N_vec[0] if len(set(N_vec)) == 1 else N // num_procs)
        
        ensemble_p /= N
        p_vector.append(ensemble_p)
    
    # Calculate final statistics
    p_hat = sum(p_vector) / K
    if K > 1:
        s_2 = sum([(p_vector[i] - p_hat)**2 for i in range(len(p_vector))]) / (K - 1)
        error = math.sqrt(s_2) / math.sqrt(K)
    else:
        error = 0
    
    print(f"Estimation finished in {time.time() - start_time:.2f} seconds.")
    print("-" * 50)
    print(f"Final probability estimate: {p_hat:.6e}")
    print(f"Standard error: {error:.6e}")
    print(f"95% confidence interval: [{p_hat - 1.96*error:.6e}, {p_hat + 1.96*error:.6e}]")

    return {
        "method": "Guided-dwSSA",
        "probability": p_hat,
        "std_error": error,
        "ci_lower": p_hat - 1.96*error,
        "ci_upper": p_hat + 1.96*error,
        "total_time": time.time() - start_time
    }


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_guidedwssa.py <config.json>")
        sys.exit(1)
    
    config_path = sys.argv[1]
    main(config_path)
