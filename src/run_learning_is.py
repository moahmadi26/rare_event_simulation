
import sys
import json
import numpy as np
import math
import time
import multiprocessing
from learning_is import LearningBasedIS, ISConfig
from utils import suppress_c_output
from prism_parser import parser

def main(json_path):
    # Load Config
    with open(json_path, 'r') as f:
        json_data = json.load(f)
        
    model_path = json_data['model_path']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    max_time = float(json_data['max_time'])
    
    # Parse model for index
    with suppress_c_output():
        model = parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    # Heuristics for params if not in JSON
    # Learning params
    dt_learning = json_data.get('dt_learning', max_time / 50.0)
    dt_estimation = json_data.get('dt_estimation', max_time / 100.0)
    num_iterations = json_data.get('num_iterations', 50)
    batch_size = json_data.get('batch_size', 10)
    learning_rate = json_data.get('learning_rate', 0.01)
    
    # Ensemble params
    K = 4
    N = 10000
    num_procs = 4
    
    config = ISConfig(
        target_index=target_index,
        target_value=target_value,
        max_time=max_time,
        dt_learning=dt_learning,
        dt_estimation=dt_estimation,
        num_iterations=num_iterations,
        batch_size=batch_size,
        learning_rate=learning_rate,
        # Default architecture
        hidden_layers=[32, 32]
    )
    
    print(f"Running Learning-Based IS...")
    print(f"Model: {model_path}")
    print(f"Target: {target_var} = {target_value}")
    print(f"Time horizon: {max_time}")
    print("-" * 50)
    
    learner = LearningBasedIS(model_path, config)
    
    start_time = time.time()
    
    # 1. Learning Phase
    print("PHASE 1: Parameter Learning")
    learner.train()
    print(f"Learning finished in {time.time() - start_time:.2f} seconds.")
    print("-" * 50)
    
    # 2. Estimation Phase
    print(f"Running estimation with K={K} ensembles of N={N} trajectories each...")
    start_time = time.time()
    
    p_vector = []
    
    # Extract NN state to pass to workers (avoid pickling the whole learner with stormpy model)
    nn_state = {
        'params': learner.nn.params,
        'mean': learner.nn.mean,
        'std': learner.nn.std,
        'input_dim': learner.nn.input_dim
    }
    
    for i in range(K):
        print(f"Ensemble {i+1}/{K}")
        
        # Split N into chunks for processors
        N_vec = [N // num_procs for _ in range(num_procs)]
        remaining = N % num_procs
        for j in range(remaining):
            N_vec[j] += 1
            
        # Pass minimal data to workers: path, config, nn_state, count
        tasks = [(model_path, config, nn_state, n) for n in N_vec]
        
        with multiprocessing.Pool(processes=num_procs) as pool:
            results = pool.starmap(run_batch, tasks)
            
        # Aggregation
        total_W = sum(r[0] for r in results) # Sum of Weights
        # Mean = Sum(W) / N
        ensemble_p = total_W / N
        p_vector.append(ensemble_p)
        
    # Stats
    p_hat = sum(p_vector) / K
    if K > 1:
        s_2 = sum([(p - p_hat)**2 for p in p_vector]) / (K - 1)
        error = math.sqrt(s_2) / math.sqrt(K)
    else:
        error = 0.0
        
    print(f"Estimation finished in {time.time() - start_time:.2f} seconds.")
    print("-" * 50)
    print(f"Final probability estimate: {p_hat:.6e}")
    print(f"Standard error: {error:.6e}")
    print(f"95% confidence interval: [{p_hat - 1.96*error:.6e}, {p_hat + 1.96*error:.6e}]")

def run_batch(model_path, config, nn_state, n_trajs):
    """Worker function to run a batch of trajectories"""
    # Re-seed random number generator for this process
    np.random.seed()
    
    # Parse model locally (avoid pickling stormpy objects)
    from prism_parser import parser
    from utils import suppress_c_output
    with suppress_c_output():
        model = parser(model_path)
    
    # Reconstruct simulator and NN
    from learning_is import TauLeapIS, NeuralNetworkAnsatz
    
    simulator = TauLeapIS(model, config)
    
    nn = NeuralNetworkAnsatz(nn_state['input_dim'], config.hidden_layers, config.activation)
    nn.params = nn_state['params']
    nn.mean = nn_state['mean']
    nn.std = nn_state['std']
    
    total_W = 0.0
    for _ in range(n_trajs):
        res = simulator.simulate_trajectory(nn, config.dt_estimation, train_mode=False)
        if res['reached']:
            total_W += np.exp(res['log_W']) # Accumulate Likelihood Ratio
            
    return (total_W, )

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python run_learning_is.py <config.json>")
        sys.exit(1)
    
    main(sys.argv[1])
