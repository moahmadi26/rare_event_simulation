# Rare Event Simulation for Stochastic Reaction Networks

This repository implements advanced importance sampling methods for estimating rare event probabilities in stochastic chemical reaction networks (CRNs). The framework provides efficient algorithms to compute probabilities of events that occur with extremely low frequency in biochemical systems.

## Methods Implemented

### 1. **dwSSA** - Dynamic Weighted Stochastic Simulation Algorithm
Cross-entropy based importance sampling method that iteratively learns optimal biasing parameters.

### 2. **dwSSA++** - Enhanced dwSSA  
Extends dwSSA with polynomial leaping acceleration for improved efficiency in large-scale systems.

### 3. **Learning-based IS** - Learning-based Importance Sampling
Uses gradient-based optimization with sigmoid ansatz functions to learn optimal importance distributions.

### 4. **Guided wSSA** - Guided Weighted Stochastic Simulation Algorithm
Implements guided importance sampling with predilection functions to steer trajectories toward rare events.

## Installation

### Prerequisites

- **Python 3.7+** with the following packages:
  - NumPy
  - Multiprocessing (standard library)
  - JSON (standard library)

- **Stormpy** - Python bindings for STORM model checker
  ```bash
  # Install Stormpy following instructions at:
  # https://moves-rwth.github.io/stormpy/installation.html
  ```

- **z3-solver** (Optional - for certain analysis features)
  ```bash
  pip install z3-solver
  ```

### Installation Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/rlwssa.git
   cd rlwssa
   ```

2. Install Python dependencies:
   ```bash
   pip install numpy stormpy
   ```

3. Verify installation:
   ```bash
   cd src
   python -c "import stormpy; import numpy; print('Installation successful!')"
   ```

## Repository Structure

```
rlwssa/
├── src/                           # Source code and main scripts
│   ├── run_dwssa.py              # Run dwSSA method
│   ├── run_dwssa_plus_plus.py    # Run dwSSA++ method
│   ├── run_learning_is.py        # Run learning-based IS
│   ├── run_guidedwssa.py         # Run guided wSSA method
│   ├── dwssa.py                  # dwSSA implementation
│   ├── dwssa_plus_plus.py        # dwSSA++ implementation
│   ├── learning_is.py            # Learning IS implementation
│   ├── guidedwssa.py             # Guided wSSA implementation
│   ├── prism_parser.py           # PRISM model parser
│   └── utils.py                  # Utility functions
└── models/                        # Biochemical reaction network models
    ├── enzymatic_futile_cycle/    # Enzymatic futile cycle model
    ├── motility_regulation/       # Motility regulation network
    ├── yeast_polarization/        # Yeast polarization model
    ├── circuit0x8E/              # Genetic circuit 0x8E
    ├── michaelis_menten/         # Michaelis-Menten kinetics
    └── pure_decay/               # Pure decay process
```

## Running the Framework

The framework accepts CTMC models written in [PRISM modeling language](https://www.prismmodelchecker.org/manual/ThePRISMLanguage/Introduction). Each method requires a JSON configuration file with the following format:

```json
{
    "model_path": "../models/model_name/model.sm",
    "model_name": "descriptive_name",
    "target_variable": "X",
    "target_value": "100",
    "max_time": "10.0",
    "description": "Optional description of the rare event"
}
```

### Parameters:
- `model_path`: Path to the PRISM model file (.sm)
- `model_name`: Arbitrary name for the model
- `target_variable`: The species/variable of interest
- `target_value`: The threshold value defining the rare event
- `max_time`: Time horizon for the simulation

### Running Different Methods

Navigate to the `src/` directory and run:

```bash
# Dynamic weighted SSA (dwSSA)
python run_dwssa.py ../models/enzymatic_futile_cycle/low_s5.json

# dwSSA with polynomial leaping (dwSSA++)
python run_dwssa_plus_plus.py ../models/enzymatic_futile_cycle/low_s5.json

# Learning-based importance sampling
python run_learning_is.py ../models/enzymatic_futile_cycle/low_s5.json

# Guided wSSA
python run_guidedwssa.py ../models/enzymatic_futile_cycle/low_s5.json
```

## Case Studies

### 1. Enzymatic Futile Cycle
A biochemical network with 6 species and 6 reactions modeling substrate cycling between phosphorylated and dephosphorylated states.
- Model: `models/enzymatic_futile_cycle/`
- Example: Low S5 concentration event

### 2. Motility Regulation
Gene regulatory network with 9 species and 12 reactions controlling bacterial motility.
- Model: `models/motility_regulation/`
- Example: High CodY production

### 3. Yeast Polarization
Model of yeast cell polarization with 7 species and 8 reactions.
- Model: `models/yeast_polarization/`
- Example: High G-protein activation

### 4. Genetic Circuit 0x8E
Complex genetic circuit with 18 species and 15 reactions using Hill functions.
- Model: `models/circuit0x8E/`
- Example: Circuit state transitions

### 5. Michaelis-Menten
Classic enzyme kinetics model with 4 species and 3 reactions.
- Model: `models/michaelis_menten/`
- Example: High product formation

### 6. Pure Decay
Simple decay process for testing and validation.
- Model: `models/pure_decay/`
- Example: Survival probability

## Model Requirements

Models must be written in PRISM language with the following restrictions:
1. All commands in the model must be labeled
2. No module renaming is allowed
3. Each command must have exactly one update

Example PRISM model:
```prism
ctmc

const double k1 = 1.0;
const double k2 = 0.1;

module reaction_network
    X : int init 100;
    Y : int init 0;
    
    [r1] X > 0 -> k1*X : (X'=X-1) & (Y'=Y+1);
    [r2] Y > 0 -> k2*Y : (Y'=Y-1);
endmodule
```

## Output

Each method provides:
- **Probability estimate**: Estimated probability of the rare event
- **Standard error**: Statistical uncertainty of the estimate
- **Confidence intervals**: 95% confidence bounds
- **Computational time**: Runtime statistics
- **Variance reduction**: Improvement over standard simulation
- **Convergence information**: Iteration-wise progress (for iterative methods)

Results are printed to console and can be saved to JSON files for further analysis.

## Advanced Configuration

### Parallel Processing
All methods support parallel execution. Adjust the number of processors in the run scripts:
```python
num_procs = 10  # Number of parallel processes
```

### Hyperparameters
Each method has specific hyperparameters that can be tuned:

- **dwSSA**: `N_train`, `rho`, `K`, `N`
- **dwSSA++**: `N_train`, `rho`, `K`, `N`, `leap_size`
- **Learning IS**: `N_train`, `L`, `K`, `N`, `learning_rate`
- **Guided wSSA**: `N_train`, `K`, `N`, `negative_method`
