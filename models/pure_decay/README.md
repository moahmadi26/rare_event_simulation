# Pure Decay Model

## Description

The Pure Decay model is a simple stochastic process representing exponential decay or death processes. This fundamental model serves as a benchmark for testing rare event simulation methods and understanding basic stochastic dynamics in chemical and biological systems.

## Model Details

### Chemical Species
- **X**: Population or molecule count undergoing decay

### Reactions

Single decay reaction:
```
R₁: X → ∅     (rate: k·X)
```

Where:
- k is the decay rate constant
- The reaction rate is proportional to the current population

### Initial Conditions

Various initial populations are tested:
- Standard: X₀ = 100
- Large population: X₀ = 1000
- Small population: X₀ = 10

## Mathematical Properties

### Deterministic Behavior
The deterministic equation is:
```
dX/dt = -kX
X(t) = X₀·e^(-kt)
```

### Stochastic Properties
- **Extinction time**: Time until X = 0
- **Survival probability**: P(X(t) > 0)
- **Half-life**: t₁/₂ = ln(2)/k
- **Mean lifetime**: τ = 1/k

## Configuration Files

Multiple scenarios available:
- `survival.json`: Survival probability estimation
- `pure_decay_survival_95.json`: 95% survival threshold
- `pure_decay_survival_50.json`: 50% survival (median)
- `pure_decay_slow_survival.json`: Slow decay regime
- `pure_decay_extinction_0.json`: Complete extinction
- `pure_decay_extinction_10.json`: Near extinction (X ≤ 10)
- `pure_decay_extinction_50.json`: Half population remaining

## Usage Example

From the `src/` directory:

```bash
# Run with dwSSA
python run_dwssa.py ../models/pure_decay/survival.json

# Run with dwSSA++
python run_dwssa_plus_plus.py ../models/pure_decay/pure_decay_survival_50.json

# Run with learning-based IS
python run_learning_is.py ../models/pure_decay/pure_decay_extinction_10.json

# Run with guided wSSA
python run_guidedwssa.py ../models/pure_decay/pure_decay_slow_survival.json
```

## Rare Events of Interest

1. **Long survival**: Population persisting beyond expected lifetime
2. **Rapid extinction**: Unusually fast decay to zero
3. **Threshold crossing**: Time to reach specific population levels
4. **Fractional survival**: Maintaining certain percentage of initial population

## Model File

The PRISM model is defined in `pure_decay.sm`

## Applications

Despite its simplicity, pure decay models:

### Physical Systems
- **Radioactive decay**: Nuclear disintegration
- **Fluorescence decay**: Excited state relaxation
- **Chemical degradation**: Molecular breakdown

### Biological Systems
- **Cell death**: Apoptosis and necrosis
- **Drug clearance**: Pharmacokinetics
- **RNA/protein degradation**: Molecular turnover
- **Population extinction**: Ecological dynamics

### Engineering
- **Reliability analysis**: Component failure
- **Queue theory**: Customer departure
- **Signal decay**: Communication systems

## Analytical Solutions

### Survival Probability
```
P(X(t) > n) = Σ(k=n+1 to X₀) (X₀ choose k) · e^(-kt·k) · (1-e^(-kt))^(X₀-k)
```

### First Passage Time
The time to reach threshold n has distribution:
```
f(t|n) = k·(X₀-n)·(X₀ choose n)·e^(-kt·n)·(1-e^(-kt))^(X₀-n-1)
```

### Extinction Time
Expected time to extinction:
```
E[T_ext] = (1/k)·Σ(i=1 to X₀) 1/i ≈ (1/k)·(ln(X₀) + γ)
```
where γ is Euler's constant.

## Stochastic Simulation

### Gillespie Algorithm
1. Calculate propensity: a = k·X
2. Generate time to next event: τ = -ln(r₁)/a
3. Decrement population: X → X-1
4. Repeat until X = 0 or t > t_max

### Importance Sampling
For rare survival events:
- Bias towards slower decay
- Weight trajectories appropriately
- Focus computational effort on rare outcomes

## Validation

This model serves as validation for:
- **Algorithm correctness**: Known analytical solutions
- **Variance reduction**: Quantifiable improvements
- **Convergence rates**: Theoretical predictions
- **Computational efficiency**: Benchmark comparisons

## Extensions

Common variations include:
- **Birth-death**: Adding birth/production terms
- **Multiple species**: Coupled decay processes
- **Time-dependent rates**: k(t)
- **Spatial effects**: Diffusion-decay
- **Environmental noise**: Stochastic rate constants

## Theoretical Importance

Pure decay is fundamental for:
- **Master equation theory**: Exact solutions
- **Large deviation theory**: Rate functions
- **Martingale methods**: Optimal stopping
- **Information theory**: Channel capacity

## References

This model is a cornerstone of stochastic process theory, appearing in classic texts by Van Kampen, Gardiner, and Gillespie on stochastic methods in physics, chemistry, and biology.