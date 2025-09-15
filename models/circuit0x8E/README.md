# Genetic Circuit 0x8E Model

## Description

Genetic Circuit 0x8E is a complex synthetic genetic circuit model consisting of 18 chemical species interacting through 15 reaction channels. This model represents an engineered biological system designed to exhibit specific computational or regulatory behaviors, using Hill functions to capture cooperativity in gene regulation.

## Model Details

### Chemical Species

The circuit contains 18 species (S₁ through S₁₈) representing various molecular components:
- **S₁ - S₈**: Primary signaling molecules (10-fold amplified in reactions)
- **S₉ - S₁₈**: Regulatory proteins/transcription factors

Note: Species S₁, S₃, S₄, S₆, and S₈ are produced in groups of 10 molecules per reaction event (abstracted representation).

### Reactions

The model uses abstracted reactions with Hill function kinetics:

```
R₁:  10S₁ → ∅                    (degradation)
R₂:  10S₃ → ∅                    (degradation)
R₃:  10S₄ → ∅                    (degradation)
R₄:  10S₆ → ∅                    (degradation)
R₅:  10S₈ → ∅                    (degradation)
R₆:  S₉ + S₂ → S₉ + S₂ + 10S₁   (catalytic production)
R₇:  S₁₀ + S₄ → S₁₀ + S₄ + 10S₁ (catalytic production)
R₈:  S₁₁ + S₄ → S₁₁ + S₄ + 10S₃ (catalytic production)
R₉:  S₁₂ + S₇ → S₁₂ + S₇ + 10S₃ (catalytic production)
R₁₀: S₁₄ + S₂ → S₁₄ + S₂ + 10S₄ (catalytic production)
R₁₁: S₁₃ + S₇ → S₁₃ + S₇ + 10S₄ (catalytic production)
R₁₂: S₁₆ + S₅ → S₁₆ + S₅ + 10S₆ (catalytic production)
R₁₃: S₁₅ + S₁ → S₁₅ + S₁ + 10S₆ (catalytic production)
R₁₄: S₁₈ + S₆ → S₁₈ + S₆ + 10S₈ (catalytic production)
R₁₅: S₁₇ + S₃ → S₁₇ + S₃ + 10S₈ (catalytic production)
```

### Initial Conditions

Default initial populations:
```
[S₁, S₂, S₃, S₄, S₅, S₆, S₇, S₈, S₉, S₁₀, S₁₁, S₁₂, S₁₃, S₁₄, S₁₅, S₁₆, S₁₇, S₁₈]
= [70, 60, 70, 0, 0, 70, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]
```

## Circuit Properties

### Design Features
- **Multilayer regulation**: Complex interconnected regulatory loops
- **Signal amplification**: 10-fold production in catalytic reactions
- **Cooperative binding**: Hill functions model molecular cooperativity
- **Modular architecture**: Distinct functional modules

### Dynamic Behaviors
- **State transitions**: Circuit can switch between distinct states
- **Signal processing**: Processes input signals to generate specific outputs
- **Memory effects**: Can maintain states over time
- **Noise filtering**: Robust to molecular fluctuations

## Configuration Files

Available configurations:
- `circuit0x8E.json`: Standard configuration
- `circuit0x8E_temp.json`: Alternative parameter set

## Usage Example

From the `src/` directory:

```bash
# Run with dwSSA
python run_dwssa.py ../models/circuit0x8E/circuit0x8E.json

# Run with dwSSA++
python run_dwssa_plus_plus.py ../models/circuit0x8E/circuit0x8E.json

# Run with learning-based IS
python run_learning_is.py ../models/circuit0x8E/circuit0x8E.json

# Run with guided wSSA
python run_guidedwssa.py ../models/circuit0x8E/circuit0x8E.json
```

## Rare Events of Interest

1. **Complete state transitions**: Moving from initial state [1,0,0] to target state [1,1,1]
2. **Signal extinction**: Complete depletion of specific species
3. **Synchronization events**: Coordinated changes across multiple species
4. **Extreme accumulation**: Unusual build-up of intermediate species
5. **Oscillatory extrema**: Peak values in oscillating components

## Model File

The PRISM model is defined in `Circuit0x8E_100to111_unb.sm`

The model name suggests a transition from state "100" to "111" (binary representation), indicating the circuit is designed to compute or switch between these logical states.

## Applications

Synthetic genetic circuits like 0x8E are used for:
- **Biocomputation**: Implementing logical operations in living cells
- **Biosensing**: Detecting and responding to environmental signals
- **Metabolic engineering**: Controlling production pathways
- **Therapeutic devices**: Smart drug delivery systems
- **Research tools**: Studying network dynamics and emergence

## Mathematical Features

### Hill Functions
The model uses Hill kinetics with parameters:
- Hill coefficients (n): Control steepness of response
- Dissociation constants (K): Set activation thresholds
- Maximum rates (k): Determine reaction speeds

### Network Topology
- **Feed-forward loops**: Rapid signal propagation
- **Feedback regulation**: Stability and homeostasis
- **Cross-talk**: Inter-module communication
- **Cascades**: Sequential activation patterns

## Computational Complexity

This circuit presents significant challenges for simulation:
- Large state space (18 species)
- Complex interdependencies
- Rare transitions between stable states
- Stiff dynamics with multiple timescales

These properties make it an excellent test case for rare event simulation methods.

## References

This model represents a class of synthetic circuits studied in synthetic biology for implementing complex behaviors in engineered organisms.