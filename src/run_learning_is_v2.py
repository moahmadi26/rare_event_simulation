
import sys
import json
import numpy as np
import math
import time
import multiprocessing
import traceback
from learning_is_v2 import LearningBasedISV2, ISConfigV2, SimulatorV2, NeuralNetworkAnsatzV2
from utils import suppress_c_output
from prism_parser import parser

def main(json_path):
    global_start_time = time.time()
    with open(json_path, 'r') as f:
        json_data = json.load(f)
        
    model_path = json_data['model_path']
    target_var = json_data['target_variable']
    target_value = int(json_data['target_value'])
    max_time = float(json_data['max_time'])
    
    # Parse for index
    with suppress_c_output():
        model = parser(model_path)
    target_index = model.species_to_index_dict[target_var]
    
    # Check condition string
    csl = json_data.get('csl_property', '')
    target_condition = 'less_equal'
    if '>=' in csl:
        target_condition = 'greater_equal'
    
    config = ISConfigV2(
        target_index=target_index,
        target_value=target_value,
        max_time=max_time,
        target_condition=target_condition,
        num_iterations=20, # Short training
        batch_size=1_000,
        learning_rate=0.01,
        hidden_layers=[20, 20]
    )
    
    print("-" * 50)
    print("Running Learning-Based IS V2 (Robust)")
    print(f"Model: {model_path}")
    print(f"Target: {target_var} {target_condition} {target_value}")
    
    # Check Initial State
    x0 = model.get_initial_state()
    val0 = x0[target_index]
    print(f"Initial State Value: {val0}")
    
    # Immediate Satisfaction Check
    sat = False
    if target_condition == 'greater_equal' and val0 >= target_value: sat = True
    if target_condition == 'less_equal' and val0 <= target_value: sat = True
    
    if sat:
        print("!!! WARNING: Target satisfied at t=0. Probability = 1.0 !!!")
        
    learner = LearningBasedISV2(model_path, config)
    
    # Phase 1: Train
    print("Starting Learning Phase...")
    start_t = time.time()
    learner.train()
    print(f"Learning done in {time.time()-start_t:.2f}s")
    
    # Phase 2: Estimate
    K = 4
    N = 10_000
    print(f"Starting Estimation (K={K}, N={N})...")
    
    p_vector = []
    
    nn_state = {
        'params': learner.nn.params,
        'mean': learner.nn.mean,
        'std': learner.nn.std,
        'input_dim': learner.nn.input_dim
    }
    
    for k in range(K):
        # Parallelize within ensemble
        num_procs = 4
        N_vec = [N // num_procs] * num_procs
        
        args = [(model_path, config, nn_state, n) for n in N_vec]
        
        with multiprocessing.Pool(num_procs) as pool:
            results = pool.starmap(run_worker, args)
            
        # Results: list of (sum_W, count)
        total_W = sum(r[0] for r in results)
        total_N = sum(r[1] for r in results) # Should be N
        
        p_est = total_W / total_N if total_N > 0 else 0.0
        p_vector.append(p_est)
        print(f"Ensemble {k+1}: {p_est:.4e}")
        
    # Stats
    p_hat = np.mean(p_vector)
    error = np.std(p_vector, ddof=1) / np.sqrt(K) if K > 1 else 0.0
    
    print("-" * 50)
    print(f"Final probability estimate: {p_hat:.6e}")
    print(f"Standard error: {error:.6e}")
    ci_low = 0.0
    ci_high = 0.0
    if p_hat > 0:
        ci_low = p_hat - 1.96*error
        ci_high = p_hat + 1.96*error
        print(f"95% confidence interval: [{ci_low:.6e}, {ci_high:.6e}]")
    else:
        print("95% confidence interval: [0.0, 0.0]")

    return {
        "method": "Learning-Based IS",
        "probability": p_hat,
        "std_error": error,
        "ci_lower": ci_low,
        "ci_upper": ci_high,
        "total_time": time.time() - global_start_time
    }

def run_worker(model_path, config, nn_state, n_trajs):
    np.random.seed()
    with suppress_c_output():
        model = parser(model_path)
    simulator = SimulatorV2(model, config)
    nn = NeuralNetworkAnsatzV2(nn_state['input_dim'], config.hidden_layers)
    nn.params = nn_state['params']
    nn.mean = nn_state['mean']
    nn.std = nn_state['std']
    
    sum_W = 0.0
    count = 0
    
    for _ in range(n_trajs):
        # Use simple estimation (no training gradients)
        res = simulator.simulate_trajectory(nn, 0.0, train_mode=False)
        if res['reached']:
            w = np.exp(res['log_W'])
            sum_W += w
        count += 1
        
    return sum_W, count

if __name__ == '__main__':
    main(sys.argv[1])
