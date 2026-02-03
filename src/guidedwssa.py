import math
import random
import numpy as np
from prism_parser import parser
from utils import get_reaction_rate, is_target, suppress_c_output


def get_propensities(model, state):
    reactions = model.get_reactions_vector()
    a = [get_reaction_rate(state, model, r_idx) 
         for r_idx in range(len(reactions))]
    a_0 = sum(a)
    return a, a_0


def calculate_propensity_matrix(model, state):
    reactions = model.get_reactions_vector()
    S_in = np.zeros((len(state), len(reactions)))
    
    for r_idx in range(len(reactions)):
        for s_idx in range(len(state)):
            if reactions[r_idx][s_idx] < 0:
                S_in[s_idx, r_idx] = abs(reactions[r_idx][s_idx])
    
    return S_in


def expectation(state, F, S, h, dt):
    e = F.T @ (np.array(state) + S @ h * dt)
    return e


def presimulation_check(state, F, S, h, dt, xp, x0):
    e = expectation(state, F, S, h, dt)
    d = np.zeros(len(e))
    FTx0 = F.T @ x0
    for i in range(len(e)):
        d[i] = (e[i] - xp[i]) * (xp[i] - FTx0[i])
    return d


def calculate_predilection(state, h, S, F, dt, xp):
    H = np.diag(h)
    e = expectation(state, F, S, h, dt)
    denominator = F.T @ S @ H @ S.T @ F * dt
    
    # Handle scalar or near-zero denominator
    if np.ndim(denominator) == 0:
        denom_val = denominator
    else:
        denom_val = denominator.item() if denominator.size == 1 else denominator[0, 0]
    
    if np.abs(denom_val) < 1e-10:
        return h
    
    h_tilde = h + ((H @ S.T @ F) / denom_val * (xp - e)).flatten()
    return h_tilde


def resolve_negatives(alpha, method='C'):
    alpha_out = alpha.copy()
    
    if method == 'A':
        if np.any(alpha_out < 0):
            alpha_out = alpha_out + abs(np.min(alpha_out)) + 0.01
        
        for i in range(len(alpha_out)):
            if np.isnan(alpha_out[i]) or alpha_out[i] < 0:
                alpha_out[i] = 1.0
                
    elif method == 'B':
        min_val = np.min(alpha_out)
        if min_val != 0:
            alpha_out = alpha_out / min_val
        
        for i in range(len(alpha_out)):
            if np.isnan(alpha_out[i]):
                alpha_out[i] = 1.0
                
    elif method == 'C':
        for i in range(len(alpha_out)):
            if np.isnan(alpha_out[i]) or alpha_out[i] < 0:
                alpha_out[i] = 1.0
    
    return alpha_out


def guidedwssa(model_path, N, t_max, target_index, target_value, F=None, negative_method='C'):
    with suppress_c_output():
        model = parser(model_path)
    
    reactions = model.get_reactions_vector()
    num_reactions = len(reactions)
    num_species = len(model.get_initial_state())
    
    S = np.zeros((num_species, num_reactions))
    for r_idx in range(num_reactions):
        for s_idx in range(num_species):
            S[s_idx, r_idx] = reactions[r_idx][s_idx]
    
    if F is None:
        F = np.zeros(num_species)
        F[target_index] = 1.0
    else:
        F = np.array(F)
    
    F = F.reshape(-1, 1)
    
    x0 = np.array(model.get_initial_state())
    xp = np.zeros(1)
    xp[0] = target_value
    
    mn = 0
    square_sum = 0
    
    for i in range(N):
        t = 0
        x = np.array(model.get_initial_state())
        w = 1
        
        while t < t_max:
            h, h0 = get_propensities(model, tuple(x))
            h = np.array(h)
            
            if is_target(tuple(x), target_index, target_value):
                mn += w
                square_sum += w**2
                break
            
            delta_t = t_max - t
            d = presimulation_check(x, F, S, h, delta_t, xp, x0)
            
            flag = 0
            for v in range(len(d)):
                if d[v] <= 0:
                    flag = 1
                    break
            
            if flag == 0:
                h_tilde = h
            else:
                h_tilde = calculate_predilection(x, h, S, F, delta_t, xp)
            
            # Safe division to avoid warnings
            with np.errstate(divide='ignore', invalid='ignore'):
                alph = np.where(h > 0, h_tilde / h, 1.0)
            
            alph = resolve_negatives(alph, negative_method)
            h_tilde = alph * h
            h0_tilde = np.sum(h_tilde)
            
            if h0_tilde <= 0:
                break
            
            r1 = random.random()
            r2 = random.random()
            
            tau = -(1.0 / h0_tilde) * math.log(r1)
            
            temp_sum = 0
            j = 0
            while temp_sum <= r2 * h0_tilde and j < len(h_tilde):
                temp_sum += h_tilde[j]
                j += 1
            j -= 1
            
            if j < 0:
                j = 0
            
            w = w * (h[j] / h_tilde[j]) * math.exp((h0_tilde - h0) * tau)
            t += tau
            x = x + S[:, j]
    
    p = mn / N if N > 0 else 0
    var = square_sum / N - p**2 if N > 0 else 0
    SE = (1 / math.sqrt(N)) * math.sqrt(var) if N > 0 and var >= 0 else 0
    zstar = 1.96
    conf = [p - zstar * SE, p + zstar * SE]
    
    return p, var, conf
