"""
Learning-Based Importance Sampling via Stochastic Optimal Control
Implementation using NumPy-based Neural Network
"""

import numpy as np
import time
from dataclasses import dataclass
from typing import List, Dict, Tuple, Optional
from utils import suppress_c_output
from prism_parser import parser

@dataclass
class ISConfig:
    """Configuration for the learning-based IS algorithm"""
    target_index: int
    target_value: float
    max_time: float
    dt_learning: float
    dt_estimation: float
    num_iterations: int
    batch_size: int
    learning_rate: float
    
    # Optional parameters with defaults
    target_condition: str = 'greater'
    hidden_layers: List[int] = None  # e.g., [20, 20]
    activation: str = 'tanh'  # 'tanh' or 'relu'
    mc_samples_final: int = 100000 
    clip_exp: float = 50.0 

    def __post_init__(self):
        if self.hidden_layers is None:
            self.hidden_layers = [20, 20]

class NeuralNetworkAnsatz:
    """
    Multilayer Perceptron (MLP) to approximate the value function V(t, x).
    V(t, x) ≈ -log(Psi(t, x)) where Psi is the committor function.
    
    We model u(t, x) = exp(-V(t, x)) directly or model V(t, x) and compute u.
    Paper tends to model V or u. Let's model V(t, x).
    
    Architecture:
    Input: [t, x_1, ..., x_d] (normalized)
    Hidden Layers: Tanh activation
    Output: Scalar V (Linear activation, but we might want to enforce non-negativity or bounds)
    """
    
    def __init__(self, input_dim: int, hidden_sizes: List[int], activation: str = 'tanh'):
        self.input_dim = input_dim
        self.hidden_sizes = hidden_sizes
        self.activation_type = activation
        
        # Initialize weights and biases
        # Structure: list of (W, b) tuples
        self.params = []
        
        # Input -> First Hidden
        prev_dim = input_dim
        for size in hidden_sizes:
            # Xavier/Glorot initialization
            limit = np.sqrt(6 / (prev_dim + size))
            W = np.random.uniform(-limit, limit, (prev_dim, size))
            b = np.zeros(size)
            self.params.append({'W': W, 'b': b})
            prev_dim = size
            
        # Last Hidden -> Output (Scalar)
        limit = np.sqrt(6 / (prev_dim + 1))
        W = np.random.uniform(-limit, limit, (prev_dim, 1))
        b = np.zeros(1)
        self.params.append({'W': W, 'b': b})
        
        # Statistics for input normalization (mean and std will be updated)
        self.mean = np.zeros(input_dim)
        self.std = np.ones(input_dim)
        
    def _activation(self, x):
        if self.activation_type == 'tanh':
            return np.tanh(x)
        elif self.activation_type == 'relu':
            return np.maximum(0, x)
        return x

    def _activation_deriv(self, x, output):
        """Derivative of activation function"""
        if self.activation_type == 'tanh':
            return 1.0 - output**2
        elif self.activation_type == 'relu':
            return (x > 0).astype(float)
        return np.ones_like(x)

    def forward(self, x_input: np.ndarray) -> Tuple[np.ndarray, List[Dict]]:
        """
        Forward pass.
        args:
            x_input: shape (batch_size, input_dim)
        returns:
            output: shape (batch_size, 1)
            cache: list of intermediate values for backprop
        """
        # Normalize input
        # Avoid division by zero
        safe_std = np.where(self.std > 1e-6, self.std, 1.0)
        a = (x_input - self.mean) / safe_std
        
        cache = [{'a': a}] # Store normalized input as first activation
        
        for i, layer in enumerate(self.params):
            W, b = layer['W'], layer['b']
            z = np.dot(a, W) + b
            
            # Apply activation for all except last layer
            if i < len(self.params) - 1:
                a = self._activation(z)
            else:
                a = z # Linear output for V
            
            cache.append({'z': z, 'a': a})
            
        return a, cache

    def evaluate(self, t: float, x: np.ndarray, T: float) -> float:
        """Evaluate V(t, x) for a single state"""
        # Input vector: [t, x_1, ..., x_d]
        input_vec = np.concatenate(([t], x)).reshape(1, -1)
        output, _ = self.forward(input_vec)
        return output[0, 0]

    def get_gradients(self, x_input: np.ndarray, grad_output: np.ndarray, cache: List[Dict]):
        """
        Backward pass to get gradients w.r.t parameters AND input x.
        
        args:
            x_input: Shape (batch_size, input_dim)
            grad_output: dL/dV, Shape (batch_size, 1). Typically 1.0 for Jacobian w.r.t V.
            cache: From forward pass
        
        returns:
            param_grads: List of dicts {'W': dW, 'b': db}
            input_grads: dV/dx, Shape (batch_size, input_dim)
        """
        batch_size = x_input.shape[0]
        params_grads = []
        
        # Backpropagate
        # Start from last layer
        delta = grad_output # (batch_size, 1)
        
        # Iterate backwards
        for i in range(len(self.params) - 1, -1, -1):
            prev_a = cache[i]['a'] # Input to this layer
            
            # Gradients for W and b
            dW = np.dot(prev_a.T, delta)
            db = np.sum(delta, axis=0)
            
            params_grads.insert(0, {'W': dW, 'b': db})
            
            # Gradient w.r.t input of this layer
            W = self.params[i]['W']
            delta = np.dot(delta, W.T) # (batch_size, prev_dim)
            
            # If not the first layer (which is input), multiply by activation derivative
            if i > 0:
                z_prev = cache[i]['z'] # This would be z of the PREVIOUS layer... wait.
                # cache[i] stores 'a' which is output of layer i-1. 
                # cache[i+1] stores 'z' and 'a' of layer i.
                # So to backprop through layer i-1 activation, we need z of layer i-1.
                # Actually, my cache structure:
                # cache[0] = input (normalized)
                # cache[1] = layer 0 output
                # ...
                # cache[i] is input to layer i.
                
                # Let's clean this up.
                # cache[i] contains 'a' which is input to layer i.
                # cache[i+1] contains 'z' and 'a' (output) of layer i.
                 
                # Backprop through activation of layer i-1 (unless i=0)
                # The 'delta' computed above is dL/da_{i-1}.
                # We need dL/dz_{i-1} = dL/da_{i-1} * sigma'(z_{i-1})
                
                # Activation of layer i-1 was applied to z_{i-1}, stored in cache[i]['z'] if i>0?
                # cache[0] is just input. logic check:
                # Layer 0: input=cache[0]['a'], output=cache[1]['a']
                # To go from Layer 0 output to Layer 0 input, we need deriv of Layer 0 activation? No.
                # We need deriv of activation of layer i-1.
                
                # Let's redo loop logic carefully.
                pass

        # Re-implementation of clean backprop
        grads = []
        delta = grad_output
        
        # Top layer (Linear)
        # Output = Last 'a'
        # Input to weights = cache[-2]['a']
        
        # Iterate from last layer index L-1 down to 0
        L = len(self.params)
        for i in range(L - 1, -1, -1):
            # i is layer index
            # Layer i takes cache[i]['a'] -> produces cache[i+1]['a']
            
            input_to_layer = cache[i]['a']
            
            # dL/dW = input.T * delta
            dW = np.dot(input_to_layer.T, delta)
            db = np.sum(delta, axis=0)
            grads.insert(0, {'W': dW, 'b': db})
            
            # Propagate delta to previous layer
            # delta_prev = delta * W.T * activation_deriv(z_prev)
            W = self.params[i]['W']
            delta_prev = np.dot(delta, W.T)
            
            if i > 0:
                # Backprop through activation of previous layer (i-1)
                # But wait, input_to_layer is a from previous layer.
                # input_to_layer = activation(z_{i-1})
                # We need z_{i-1}. That is stored in cache[i]['z']? 
                # No, cache[i] only has 'a' for i=0. But for i>0, cache[i] has 'z' and 'a'.
                z_prev = cache[i]['z'] 
                delta_prev = delta_prev * self._activation_deriv(z_prev, input_to_layer)
            
            delta = delta_prev
            
        # The final 'delta' is gradient w.r.t. normalized input
        safe_std = np.where(self.std > 1e-6, self.std, 1.0)
        input_grads = delta / safe_std
        
        return grads, input_grads


class TauLeapIS:
    """Tau-leap simulation with IS control via Neural Network"""
    
    def __init__(self, model, config: ISConfig):
        self.model = model
        self.config = config
        self.num_reactions = len(model.get_reactions_vector())
        self.num_species = len(model.get_initial_state())
        self.stoichiometry = np.array(model.get_reactions_vector())
        
    def compute_propensities(self, x: np.ndarray) -> np.ndarray:
        a = np.zeros(self.num_reactions)
        x_tuple = tuple(max(0, int(xi)) for xi in x)
        for j in range(self.num_reactions):
            a[j] = self.model.get_reaction_rate(x_tuple, j)
        return a
    
    def simulate_trajectory(self, nn: NeuralNetworkAnsatz, dt: float, train_mode: bool = False) -> Dict:
        """
        Simulate trajectory.
        If train_mode=True, we accumulate gradients for SOC objective.
        But for "Learning-Based IS", they often use the Cross-Entropy method or minimizing Variance.
        The paper minimizes variance of the estimator. 
        Usually via pathwise gradients or simply fitting V to learn the committor.
        
        The paper section 2.3 mentions: 
        Minimizing E[(W * 1_E)^2] or similar.
        Actually, they often use a specific loss function.
        
        Here, we implement the "Minimizing Second Moment" approach or simply "Learning the Committor".
        The paper "Learning-Based IS" typically uses the fact that optimal control u* = -grad V.
        
        Let's implement the simulation with IS controls derived from V.
        Control lambda_j = a_j * exp( -0.5 * (V(x + nu_j) - V(x)) ) ? 
        Or Eq 2.16: delta_j = a_j * sqrt( u(x+nu)/u(x) ). if u = exp(-V), then sqrt(ratio) = exp( -0.5 * (V_next - V_curr) ).
        """
        T = self.config.max_time
        x = np.array(self.model.get_initial_state(), dtype=float)
        t = 0.0
        
        log_W = 0.0
        states = [x.copy()]
        
        # Accumulate gradients w.r.t parameters?
        # Implementing full backprop through time (BPTT) or pathwise grad is complex.
        # Often simple policy gradient or "Galilean" updates are used.
        # The implementation in the previous file used a direct gradient of log_W.
        # We will track trajectory for batch learning if in train mode.
        
        trajectory_data = [] # List of tuples (t, x, delta, a, K_j)
        
        steps = 0
        max_steps = int(2 * T / dt) + 100 # Safety
        
        while t < T and steps < max_steps:
            current_dt = min(dt, T - t)
            
            a = self.compute_propensities(x)
            
            # Compute V(t, x) and V(t, x + nu_j) to get controls
            # We need efficient evaluation.
            
            # Current V
            input_curr = np.concatenate(([t], x)).reshape(1, -1)
            V_curr, _ = nn.forward(input_curr)
            V_curr = V_curr[0,0]
            
            delta = np.zeros(self.num_reactions)
            
            for j in range(self.num_reactions):
                if a[j] > 1e-10:
                    x_next = np.maximum(0, x + self.stoichiometry[j])
                    
                    # Next V estimate (approximate at current t for control calc, or t+dt?)
                    # Usually control is instantaneous.
                    input_next = np.concatenate(([t], x_next)).reshape(1, -1)
                    V_next, _ = nn.forward(input_next)
                    V_next = V_next[0,0]
                    
                    # u(x) = exp(-V(x)). 
                    # ratio u_next/u_curr = exp(V_curr - V_next)
                    # delta = a * ratio = a * exp(V_curr - V_next)
                    
                    exponent = (V_curr - V_next)
                    # Clip exponent
                    exponent = np.clip(exponent, -20.0, 20.0) 
                    
                    delta[j] = a[j] * np.exp(exponent)
                else:
                    delta[j] = 0.0
            
            # Sample reactions (Tau-leap)
            # K_j ~ Poisson(delta_j * dt)
            K = np.zeros(self.num_reactions, dtype=int)
            for j in range(self.num_reactions):
                if delta[j] > 0:
                    mean = delta[j] * current_dt
                    # Safety check: Clip mean to avoid "lam value too large" error
                    # Numpy's random.poisson can handle up to large floats but fails on Inf or overflows
                    if mean > 1e9:
                         # For very large mean, Poisson converges to Normal(mean, sqrt(mean))
                         # But let's just clip for stability or take a robust sample
                         # Or just cap it if it's unreasonably high (reaction explosion)
                         mean = 1e9
                    
                    K[j] = np.random.poisson(mean)
            
            # Update Weight
            # log W += sum ( K log(a/delta) + (delta - a)dt )
            step_log_W = 0.0
            for j in range(self.num_reactions):
                if delta[j] > 1e-10:
                    term1 = (delta[j] - a[j]) * current_dt
                    term2 = 0.0
                    if K[j] > 0:
                        if a[j] > 0:
                            term2 = K[j] * np.log(a[j] / delta[j])
                        else:
                            step_log_W = -np.inf # Impossible path
                            break
                    step_log_W += term1 + term2
            
            log_W += step_log_W
            
            if train_mode:
                trajectory_data.append({
                    't': t, 'x': x.copy(), 'a': a, 'delta': delta, 'K': K, 'dt': current_dt
                })
            
            # Update state
            for j in range(self.num_reactions):
                if K[j] > 0:
                    x += K[j] * self.stoichiometry[j]
            x = np.maximum(0, x)
            t += current_dt
            steps += 1
            states.append(x.copy())
            
        reached = self._check_target(states)
        
        return {
            'reached': reached,
            'log_W': log_W,
            'trajectory': trajectory_data,
            'final_state': x
        }
        
    def _check_target(self, states):
        idx = self.config.target_index
        val = self.config.target_value
        cond = self.config.target_condition
        
        for x in states:
            if cond == 'greater' or cond == 'greater_equal':
                if x[idx] >= val: return True
            elif cond == 'less' or cond == 'less_equal':
                if x[idx] <= val: return True
        return False

class LearningBasedIS:
    
    def __init__(self, model_path: str, config: ISConfig):
        with suppress_c_output():
            self.model = parser(model_path)
            
        self.config = config
        self.simulator = TauLeapIS(self.model, config)
        
        # Determine input stats from pilot run
        x0 = np.array(self.model.get_initial_state())
        input_dim = 1 + len(x0) # t + state
        
        self.nn = NeuralNetworkAnsatz(input_dim, config.hidden_layers, config.activation)
        
        # Initialize normalization stats with x0
        self.nn.mean = np.concatenate(([0], x0))
        self.nn.std = np.concatenate(([config.max_time], np.maximum(x0, 1.0)))

    def train(self):
        """
        Train using Cross-Entropy inspired method or similar SOC gradient descent.
        Since we don't have auto-diff, we use a simple specific gradient update derived for this problem
        or a generic gradient-free method (CEM) or approximate gradient.
        
        Given the constraints and lack of torch, implementing exact SOC gradient (pathwise) is hard.
        However, the paper "Learning-based IS" often uses:
        Loss = E[ (W * 1_E)^2 ] -> Minimize second moment.
        
        Alternative: Cross-Entropy Method.
        1. Generate samples with current control.
        2. Select "elite" samples (those that reach target or close to it).
        3. Fit V(t,x) to -log(optimal_weight) or similar.
        
        Let's try a simplified "Weight-weighted" regression loop which is common in CE-IS.
        Objective: Minimize KL divergence between u_theta and u_optimal.
        u_optimal(t,x) approx probability of reaching target.
        
        Iterative algorithm:
        1. Run N trajectories with current NN.
        2. Assign "value" to each state visited based on whether the trajectory reached target.
           Value(x_k) should be related to the "Committor Probability" from x_k.
           Estimator: P(reach | x_k).
        3. Train NN to predict this Value.
        
        Actually, let's implement the standard Policy Gradient style update from the provided `learning_is.py` logic
        but adapted for NN. The previous code used `delta_j` gradients. 
        
        We will use a simpler approach for NumPy: **Cross-Entropy with Elite Samples**.
        This is robust and derivative-free for the trajectory part, only supervised learning for NN.
        
        Algorithm:
        Loop:
          1. Sample N trajectories using current NN.
          2. Identify "Elite" trajectories (reached target). If too few, pick closest ones.
          3. For each state (t, x) in elite trajectories:
             Target Value y = 1.0 (or decaying as we get further from success? No, just 1.0 is simpler for committor).
             Wait, committor is P(reach). 
             Better: For elite trajectories, we want V(t, x) to encourage this path.
             
             Let's blindly follow the gradients of the "Variance" loss if possible? No.
             
             Let's use the "Many-to-One" approach.
             For every successful trajectory i, with weight W_i.
             We want to maximize E[ W * 1_E ]. But we are minimizing Variance.
             The optimal control is related to the zero-variance change of measure.
             
             Let's use a **Supervised Learning** approach on Elite Paths.
             States on successful paths should have high "Committor" value.
             States on failed paths have low value.
             
             Regression Target:
             For trajectory j, let Z_j = 1 if reached else 0.
             If we use just elites, we set Z_j = 1.
             We want V(t, x) approx -log( P(reach | t, x) ).
             
             We can estimate P(reach | t, x) via simulation? No, too expensive.
             
             Let's use the **Cross-Entropy Method (CEM)** on parametrised policy.
             We want to minimize Cross Entropy between current sampling distribution and the optimal zero-variance one.
             The gradient is E_u [ (1_E * W) * grad( log u_theta ) ].
             
             We can compute grad( log u_theta ) with backprop.
             
             Control delta_j = a_j * exp(0.5*(V_curr - V_next)).
             Likelihood ratio involve delta.
             
             Let's stick to a simpler heuristic for stability without Autodiff:
             **Train V(t, x) to approximate the "distance" or "probability" via regression.**
             
             Heuristic:
             1. Gather batch.
             2. Select top rho percentile of trajectories (closest to target).
             3. For these elite trajectories, we want to decrease the "energy" V(t,x).
             4. Define Target V*(t, x) based on "progress".
                E.g. V*(t, x) = -k * (Target - Current). Scaling to 0 at target.
                This effectively learns a biasing potential.
             
             Wait, the user wants us to implement the paper.
             The paper uses "Stochastic Approximation" with "Adam".
             The gradient is derived in Lemma 2.9 (as seen in previous code).
             
             Gradient of J(theta) = Variance.
             Grad J = E [ W^2 * 1_E * grad(log W) ].
             Wait, grad(Variance) involves knowing the constant P.
             Usually we minimize Second Moment E[ (W * 1_E)^2 ].
             Gradient ~ E [ W^2 * 1_E * grad(log W) ]. 
             NOTE: W depends on theta.
             grad(W) = W * grad(log W).
             So, Gradient ~ E [ W * 1_E * W * grad(log W) ].
             Weighted average of score function.
             
             We can implement this!
             
             We need grad(log W) w.r.t parameters theta.
             log W = Sum [ K_j * log(a/delta) + (delta - a)dt ]
                   = Sum [ K_j * (log a - log delta) + delta*dt - a*dt ]
             
             d(log W)/dTheta = Sum [ K_j * (-1/delta * dDelta/dTheta) + dt * dDelta/dTheta ]
                             = Sum [ (dt - K_j/delta) * dDelta/dTheta ]
             
             Recall delta_j = a_j * exp(0.5 * (V_curr - V_next)).
             log delta_j = log a_j + 0.5 * (V_curr - V_next).
             
             dDelta/dTheta = delta * d(log delta)/dTheta 
                           = delta * 0.5 * (dV_curr/dTheta - dV_next/dTheta).
             
             Plug back:
             d(log W)/dTheta = Sum [ (dt * delta - K_j) * 0.5 * (dV_curr/dTheta - dV_next/dTheta) ].
             
             (dV/dTheta) is obtained via backprop on the NN.
        """
        
        optimizer = AdamOptimizer(self.nn.params, lr=self.config.learning_rate)
        
        for iteration in range(self.config.num_iterations):
            batch_loss = []
            
            # Run batch
            # We process one by one or in mini-batches? 
            # NumPy NN can handle batches, but simulation is trajectory-based.
            # We'll run simulation, accumulating "gradients" per trajectory, then average.
            
            # Storage for batch update
            batch_grads = [] # List of param_grad dicts
            
            success_count = 0
            
            for b in range(self.config.batch_size):
                # 1. Simulate
                res = self.simulator.simulate_trajectory(self.nn, self.config.dt_learning, train_mode=True)
                
                # 2. Compute Trajectory Gradient
                if res['reached']:
                    success_count += 1
                    W = res['log_W']
                    # Numerical stability? W can be tiny or huge.
                    # Objective: Minimize E[ (W 1_E)^2 ].
                    # Grad = 2 * (W 1_E) * d(W 1_E)/dTheta = 2 * W^2 * 1_E * d(log W)/dTheta.
                    # We can drop the factor 2.
                    # Weight for this trajectory gradient is W^2.
                    
                    real_W = np.exp(res['log_W'])
                    weight = real_W**2
                    
                    # Accumulate grad(log W)
                    # We need to sum over time steps
                    
                    traj_grad_acc = None # Structure matching self.nn.params
                    
                    # This implies BPTT or similar cost.
                    # Doing this for every trajectory in Python loop will be slow.
                    # Optimization: Accumulate tuples (t, x, coeff) then do one big NN backward pass?
                    # Coeff for step k, reaction j: C_{kj} = 0.5 * (dt * delta - K_j).
                    # Contribution: C_{kj} * (dV(t_k)/dTheta - dV(t_{k+1})/dTheta).
                    
                    # Regroup by V(t_k):
                    # V(t_k) is involved in step k (as V_curr) and step k-1 (as V_next).
                    # Coeff for V(t_k):
                    # From step k: Sum_j C_{kj} * (+1)
                    # From step k-1: Sum_j C_{k-1,j} * (-1)
                    
                    # Let Beta_k = Sum_j (dt*delta_kj - K_kj).
                    # d(log W)/dTheta = Sum_k Beta_k * (dV_k - dV_{k+1})
                    
                    # So we calculate coeff "D_k" for each time step k.
                    
                    trajectory = res['trajectory']
                    
                    inputs = []
                    coeffs = []
                    
                    for k, step in enumerate(trajectory):
                        # Step k has t, x, delta, K, dt
                        # Beta_k = Sum( delta*dt - K )
                        beta_k = np.sum( step['delta'] * step['dt'] - step['K'] )
                        
                        # We need to assign coefficients to V(t_k) and V(t_{k+1})
                        
                        # V_curr (at t_k, x_k). Coeff: +beta_k
                        inputs.append( np.concatenate(([step['t']], step['x'])) )
                        coeffs.append( beta_k )
                        
                        # V_next.
                        # Sum_j [ (dt*delta_j - K_j) * (dV(x)/dTheta - dV(x+nu_j)/dTheta) ]
                        # Terms:
                        # For x: Sum_j (term_j)
                        # For x+nu_j: -term_j
                        
                        term_j = (step['delta'] * step['dt'] - step['K'])
                        
                        # Point x (re-adding contribution from V_next part?)
                        # Actually simpler:
                        # d(log delta)/dTheta = dV_curr - dV_next
                        # d(log W)/dTheta = Sum_j (delta*dt - K) * (dV_curr - dV_next)
                        #                 = Beta_k * dV_curr - Sum_j (delta*dt - K)_j * dV_next_j
                        
                        # So for V_curr (x): Coeff is Beta_k
                        # For V_next (x+nu_j): Coeff is -(delta*dt - K)_j
                        
                        # inputs.append(...) already done for V_curr above?
                        # Let's just collect all terms.
                        
                        # Term for V(x) from this step
                        coeff_x = np.sum(term_j)
                        inputs.append( np.concatenate(([step['t']], step['x'])) )
                        coeffs.append( coeff_x )
                        
                        # Terms for V(x+nu_j)
                        for j in range(self.simulator.num_reactions):
                            if np.abs(term_j[j]) > 1e-12:
                                x_next = np.maximum(0, step['x'] + self.simulator.stoichiometry[j])
                                inputs.append( np.concatenate(([step['t']], x_next)) )
                                coeffs.append( -term_j[j] )

                    # Now batch backprop
                    if not inputs: continue
                    
                    inputs = np.array(inputs)
                    coeffs = np.array(coeffs) * weight # Scale by W^2
                    
                    # Forward pass to populate cache
                    _, cache_list = self.nn.forward(inputs)
                    # We only need cache for backward. 
                    # Actually forward() returns 'a' and 'cache'.
                    # cache contains 'a' and 'z' for all layers.
                    
                    # Backward
                    # grad_output is dLoss/dV per sample.
                    # Here Loss is linear combination of V's.
                    # So grad_output = coeffs.
                    grad_output = coeffs.reshape(-1, 1)
                    
                    param_grads_traj, _ = self.nn.get_gradients(inputs, grad_output, cache_list)
                    
                    # Add to batch accumulator
                    if not batch_grads:
                        # Initialize
                        batch_grads = param_grads_traj
                    else:
                        # Sum parameters
                        for l in range(len(batch_grads)):
                            batch_grads[l]['W'] += param_grads_traj[l]['W']
                            batch_grads[l]['b'] += param_grads_traj[l]['b']
                            
            # End of batch
            if batch_grads:
                # Average? Or Sum?
                # Expectation is average.
                scale = 1.0 / self.config.batch_size
                for l in range(len(batch_grads)):
                    batch_grads[l]['W'] *= scale
                    batch_grads[l]['b'] *= scale
                    
                optimizer.step(batch_grads)
                
            print(f"Iteration {iteration+1}/{self.config.num_iterations}: Success {success_count}/{self.config.batch_size}")

class AdamOptimizer:
    def __init__(self, params, lr=0.01, beta1=0.9, beta2=0.999, epsilon=1e-8):
        self.params = params
        self.lr = lr
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon
        self.m = []
        self.v = []
        self.t = 0
        
        for layer in params:
            self.m.append({'W': np.zeros_like(layer['W']), 'b': np.zeros_like(layer['b'])})
            self.v.append({'W': np.zeros_like(layer['W']), 'b': np.zeros_like(layer['b'])})
            
    def step(self, grads):
        self.t += 1
        for i, (param, grad) in enumerate(zip(self.params, grads)):
            for key in ['W', 'b']:
                g = grad[key]
                self.m[i][key] = self.beta1 * self.m[i][key] + (1-self.beta1) * g
                self.v[i][key] = self.beta2 * self.v[i][key] + (1-self.beta2) * g**2
                
                m_hat = self.m[i][key] / (1 - self.beta1**self.t)
                v_hat = self.v[i][key] / (1 - self.beta2**self.t)
                
                param[key] -= self.lr * m_hat / (np.sqrt(v_hat) + self.epsilon)
