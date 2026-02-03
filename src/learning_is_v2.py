"""
Learning-Based Importance Sampling via Stochastic Optimal Control (v2)
Robust Implementation using Cross-Entropy / Variance Minimization
"""

import numpy as np
import time
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from utils import suppress_c_output
from prism_parser import parser

@dataclass
class ISConfigV2:
    """Configuration for the learning-based IS algorithm"""
    target_index: int
    target_value: float
    max_time: float
    dt_learning: float = 0.1
    dt_estimation: float = 0.1
    num_iterations: int = 50
    batch_size: int = 20
    learning_rate: float = 0.005
    target_condition: str = 'greater_equal'
    hidden_layers: List[int] = None
    
    def __post_init__(self):
        if self.hidden_layers is None:
            self.hidden_layers = [32, 32]

class NeuralNetworkAnsatzV2:
    """
    MLP for V(t, x).
    Input: [t, x_normalized]
    Output: Scalar V
    """
    def __init__(self, input_dim: int, hidden_sizes: List[int], activation: str = 'tanh'):
        self.input_dim = input_dim
        self.params = []
        
        # Init weights
        prev_dim = input_dim
        for size in hidden_sizes:
            limit = np.sqrt(6 / (prev_dim + size))
            W = np.random.uniform(-limit, limit, (prev_dim, size))
            b = np.zeros(size)
            self.params.append({'W': W, 'b': b})
            prev_dim = size
            
        # Last layer
        limit = np.sqrt(6 / (prev_dim + 1))
        W = np.random.uniform(-limit, limit, (prev_dim, 1))
        b = np.zeros(1)
        self.params.append({'W': W, 'b': b})
        
        self.mean = np.zeros(input_dim)
        self.std = np.ones(input_dim)
        
    def forward(self, x_input: np.ndarray) -> Tuple[np.ndarray, List[Dict]]:
        # Normalize
        safe_std = np.where(self.std > 1e-6, self.std, 1.0)
        a = (x_input - self.mean) / safe_std
        
        cache = [{'a': a}] 
        
        for i, layer in enumerate(self.params):
            W, b = layer['W'], layer['b']
            z = np.dot(a, W) + b
            if i < len(self.params) - 1:
                a = np.tanh(z)
            else:
                a = z 
            cache.append({'z': z, 'a': a})
            
        return a, cache

    def get_gradients(self, x_input: np.ndarray, grad_output: np.ndarray, cache: List[Dict]):
        # Backprop
        grads = []
        delta = grad_output
        
        for i in range(len(self.params) - 1, -1, -1):
            input_to_layer = cache[i]['a']
            W = self.params[i]['W']
            
            dW = np.dot(input_to_layer.T, delta)
            db = np.sum(delta, axis=0)
            grads.insert(0, {'W': dW, 'b': db})
            
            delta_prev = np.dot(delta, W.T)
            
            if i > 0:
                # Tanh derivative: 1 - a^2
                activ = cache[i]['a'] # Output of prev layer? No, input to this layer was 'a' from prev
                # Wait: cache[i] is input to layer i. cache[i+1] is output.
                # delta_prev is dL/da_{i-1}.
                # We need dL/dz_{i-1}.
                # a_{i-1} = tanh(z_{i-1}). da/dz = 1 - a^2.
                # So delta = delta_prev * (1 - input_to_layer**2)
                
                # Careful with indices.
                # Layer i takes cache[i]['a'] -> produces cache[i+1]['a'].
                # To backprop through layer i-1's activation:
                # We need deriv of cache[i]['a'] w.r.t cache[i]['z'].
                # But we don't store z of prev layer explicitly in this loop easily?
                # Actually cache[i] stores 'a'. If i>0, that 'a' came from layer i-1.
                # We need z_{i-1} or just use 1 - a_{i-1}^2.
                
                delta_prev = delta_prev * (1.0 - input_to_layer**2)
                
            delta = delta_prev
            
        return grads

class SimulatorV2:
    def __init__(self, model, config: ISConfigV2):
        self.model = model
        self.config = config
        self.num_reactions = len(model.get_reactions_vector())
        self.stoichiometry = np.array(model.get_reactions_vector())
        
    def compute_rates(self, x):
        rates = np.zeros(self.num_reactions)
        x_tuple = tuple(max(0, int(xi)) for xi in x)
        for j in range(self.num_reactions):
            rates[j] = self.model.get_reaction_rate(x_tuple, j)
        return rates

    def simulate_trajectory(self, nn: NeuralNetworkAnsatzV2, dt: float, train_mode=False):
        T = self.config.max_time
        x = np.array(self.model.get_initial_state(), dtype=float)
        t = 0.0
        
        log_W = 0.0
        trajectory = []
        reached = False
        
        while t < T:
            # 1. Compute nominal rates
            a = self.compute_rates(x)
            lambda0 = np.sum(a)
            
            if lambda0 < 1e-12:
                # Absoring state
                break
                
            # 2. Compute Controls (Twisted rates)
            # Need V(x) and V(x+nu)
            
            # Batch inputs for efficiency: current x and all neighbors
            states_to_eval = [np.concatenate(([t], x))]
            for j in range(self.num_reactions):
                # Always eval neighbors to get gradient later? 
                # Or only if a[j] > 0?
                if a[j] > 0:
                    nu = self.stoichiometry[j]
                    states_to_eval.append(np.concatenate(([t], np.maximum(0, x+nu))))
                else:
                    # Dummy
                    states_to_eval.append(np.concatenate(([t], x)))
            
            inputs = np.array(states_to_eval)
            V_vals, _ = nn.forward(inputs)
            V_curr = V_vals[0, 0]
            
            a_tilde = np.zeros_like(a)
            log_ratios = np.zeros_like(a)
            
            for j in range(self.num_reactions):
                if a[j] > 0:
                    V_next = V_vals[j+1, 0]
                    # Control: a_tilde = a * exp(V_c - V_n)
                    diff = V_curr - V_next
                    # Clip for stability
                    diff = np.clip(diff, -10.0, 10.0)
                    a_tilde[j] = a[j] * np.exp(diff)
                    log_ratios[j] = diff
            
            lambda_tilde = np.sum(a_tilde)
            
            # 3. Step (SSA / Tau-Leap). Using Direct Method (SSA) for accuracy in rare events usually.
            # Paper uses Tau-Leap. Let's use SSA for cleaner gradients if T is small?
            # Config has dt. If dt is small step, it's basically SSA.
            # Actually, let's stick to SSA-like step for exact timing if allowed.
            # Or constant time step? The user config asks for `dt`.
            # Let's use Gillespie (SSA) since rare events are sensitive to discrete steps. Target `U<=T`.
            
            # Sampling next reaction time tau
            # In twisted dynamics: tau ~ Exp(lambda_tilde)
            if lambda_tilde < 1e-12:
                break
                
            r1 = np.random.random()
            tau = -np.log(r1) / lambda_tilde
            
            if t + tau > T:
                t = T
                break
            
            t += tau
            
            # Sampling reaction mu
            r2 = np.random.random() * lambda_tilde
            mu = -1
            cum = 0.0
            for j in range(self.num_reactions):
                cum += a_tilde[j]
                if r2 <= cum:
                    mu = j
                    break
            
            # 4. Update Weight
            # log W += log(P_orig / P_new)
            # P_orig_step = a_mu * exp(-lambda0 * tau)
            # P_new_step = a_tilde_mu * exp(-lambda_tilde * tau)
            # Ratio = (a/a_tilde) * exp( (lambda_tilde - lambda0) * tau )
            
            log_step_W = np.log(a[mu]) - np.log(a_tilde[mu]) + (lambda_tilde - lambda0) * tau
            log_W += log_step_W
            
            if train_mode:
                # Store data for gradient
                trajectory.append({
                    't': t - tau,
                    'x': x.copy(),
                    'mu': mu,
                    'a': a,
                    'a_tilde': a_tilde,
                    'lambda_tilde': lambda_tilde,
                    'tau': tau,
                    'log_ratios': log_ratios, # V_c - V_n
                })
            
            # 5. Update State
            x = np.maximum(0, x + self.stoichiometry[mu])
            
            # Check target
            if self.check_target(x):
                reached = True
                break
                
        return {'reached': reached, 'log_W': log_W, 'trajectory': trajectory}

    def check_target(self, x):
        idx = self.config.target_index
        val = x[idx]
        tval = self.config.target_value
        cond = self.config.target_condition
        if cond == 'greater_equal': return val >= tval
        if cond == 'less_equal': return val <= tval
        return val >= tval

class LearningBasedISV2:
    def __init__(self, model_path, config: ISConfigV2):
        with suppress_c_output():
            self.model = parser(model_path)
            
        self.config = config
        self.simulator = SimulatorV2(self.model, config)
        
        x0 = np.array(self.model.get_initial_state())
        input_dim = 1 + len(x0)
        self.nn = NeuralNetworkAnsatzV2(input_dim, config.hidden_layers)
        
        # Norm stats
        self.nn.mean = np.concatenate(([0], x0))
        self.nn.std = np.concatenate(([config.max_time], np.maximum(1.0, x0)))

    def train(self):
        optimizer_params = self.nn.params
        lr = self.config.learning_rate
        
        # Adam params
        m = [{'W': 0, 'b': 0} for _ in optimizer_params]
        v = [{'W': 0, 'b': 0} for _ in optimizer_params]
        beta1, beta2 = 0.9, 0.999
        eps = 1e-8
        
        for iteration in range(self.config.num_iterations):
            batch_grads = []
            valid_trajs = 0
            
            # Collect batch (parallelize if needed, here serial for safety)
            trajectories = []
            for _ in range(self.config.batch_size):
                res = self.simulator.simulate_trajectory(self.nn, 0.0, train_mode=True)
                trajectories.append(res)
            
            # Compute Gradient minimizing Second Moment ? Or Likelihood Ratio?
            # loss = E[ (W * 1_E)^2 ]
            # grad approx = - E[ 1_E * W^2 * grad(log P_tilde) ]
            
            total_grad = [{'W': np.zeros_like(p['W']), 'b': np.zeros_like(p['b'])} for p in optimizer_params]
            
            count_success = 0
            
            for res in trajectories:
                if not res['reached']:
                    continue
                count_success += 1
                
                log_W = res['log_W']
                # Clip W to avoid explosion
                # log_W can be +100 or -100.
                # If W is huge, this trajectory dominates.
                real_W = np.exp(log_W)
                
                weight = -(real_W**2) # Negative because we minimize, but derivation says -E[...].
                # Wait, minimize J -> subtract grad.
                # grad J = - E [ ... ].
                # So update = - lr * grad J = + lr * E [ ... ].
                # Effectively we maximize likelihood of these weighted paths.
                
                traj_grad = self.compute_trajectory_score_grad(res['trajectory'])
                
                # Accumulate
                for l in range(len(total_grad)):
                    total_grad[l]['W'] += weight * traj_grad[l]['W']
                    total_grad[l]['b'] += weight * traj_grad[l]['b']
            
            if count_success > 0:
                # Normalize?
                # Usually 1/N
                scale = 1.0 / self.config.batch_size
                
                # Adam Step
                for i in range(len(optimizer_params)):
                    gw = total_grad[i]['W'] * scale
                    gb = total_grad[i]['b'] * scale
                    
                    # Clip gradients
                    gw = np.clip(gw, -1.0, 1.0)
                    gb = np.clip(gb, -1.0, 1.0)
                    
                    m[i]['W'] = beta1 * m[i]['W'] + (1-beta1) * gw
                    v[i]['W'] = beta2 * v[i]['W'] + (1-beta2) * gw**2
                    m[i]['b'] = beta1 * m[i]['b'] + (1-beta1) * gb
                    v[i]['b'] = beta2 * v[i]['b'] + (1-beta2) * gb**2
                    
                    m_hat_w = m[i]['W'] / (1 - beta1**(iteration+1))
                    v_hat_w = v[i]['W'] / (1 - beta2**(iteration+1))
                    
                    m_hat_b = m[i]['b'] / (1 - beta1**(iteration+1))
                    v_hat_b = v[i]['b'] / (1 - beta2**(iteration+1))
                    
                    optimizer_params[i]['W'] -= lr * m_hat_w / (np.sqrt(v_hat_w) + eps)
                    optimizer_params[i]['b'] -= lr * m_hat_b / (np.sqrt(v_hat_b) + eps)
            
            print(f"Iter {iteration+1}: {count_success}/{self.config.batch_size} hit target.")

    def compute_trajectory_score_grad(self, trajectory):
        """
        Compute sum of grad(log P_tilde_step)
        log P_tilde_step = log a_tilde_mu - lambda_tilde * tau
        
        grad(log a_tilde_mu) = grad(log a + V(x) - V(x+nu)) = grad(V(x)) - grad(V(x+nu))
        grad(lambda_tilde) = sum_j a_tilde_j * (grad(V(x)) - grad(V(x+nu_j)))
        """
        grads_by_layer = [{'W': np.zeros_like(p['W']), 'b': np.zeros_like(p['b'])} for p in self.nn.params]
        
        # We process step by step? 
        # Batched backprop is better.
        
        # Collect all points where we need gradients: x and all neighbors (x+nu)
        # For each step k:
        #   We need grad V(x_k) with coeff:
        #      +1 (from log a_tilde_mu part)
        #      - tau * sum_j a_tilde_j (from lambda_tilde part)
        #   We need grad V(x_k + nu_mu) with coeff:
        #      -1
        #   We need grad V(x_k + nu_j) for all j with coeff:
        #      + tau * a_tilde_j
        
        all_inputs = []
        all_coeffs = []
        
        for step in trajectory:
            t = step['t']
            x = step['x']
            mu = step['mu']
            tau = step['tau']
            a_tilde = step['a_tilde']
            
            # Terms Analysis
            # 1. log a_tilde_mu = ... + V(x) - V(x+nu_mu)
            # 2. - tau * lambda_tilde = - tau * sum_j a_tilde_j [ ... + V(x) - V(x+nu_j) ]
            
            # Coeff for V(x): 1 - tau * sum(a_tilde) = 1 - tau * lambda_tilde
            coeff_x = 1.0 - tau * step['lambda_tilde']
            all_inputs.append(np.concatenate(([t], x)))
            all_coeffs.append(coeff_x)
            
            # Coeff for V(x + nu_mu) from part 1: -1
            nu_mu = self.simulator.stoichiometry[mu]
            all_inputs.append(np.concatenate(([t], np.maximum(0, x+nu_mu))))
            all_coeffs.append(-1.0)
            
            # Coeff for V(x + nu_j) from part 2: + tau * a_tilde_j
            for j in range(len(a_tilde)):
                if a_tilde[j] > 1e-12:
                     nu_j = self.simulator.stoichiometry[j]
                     all_inputs.append(np.concatenate(([t], np.maximum(0, x+nu_j))))
                     all_coeffs.append(tau * a_tilde[j])
        
        if not all_inputs:
            return grads_by_layer
            
        inputs_arr = np.array(all_inputs)
        coeffs_arr = np.array(all_coeffs).reshape(-1, 1)
        
        # Backprop
        _, cache = self.nn.forward(inputs_arr)
        step_grads = self.nn.get_gradients(inputs_arr, coeffs_arr, cache)
        
        # Sum up
        for l in range(len(grads_by_layer)):
            grads_by_layer[l]['W'] += step_grads[l]['W']
            grads_by_layer[l]['b'] += step_grads[l]['b']
            
        return grads_by_layer
