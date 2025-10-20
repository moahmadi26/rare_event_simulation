# Repressilator Model

## Description

The Repressilator is a landmark synthetic genetic oscillator first constructed by Elowitz and Leibler in 2000. This three-gene regulatory network creates sustained oscillations through cyclic repression, representing one of the first successful implementations of a biological clock. The circuit demonstrates how simple regulatory motifs can generate complex temporal dynamics in living cells.

## Model Details

### System Architecture

The Repressilator consists of three repressor proteins arranged in a ring topology:

#### Core Components
- **LacI**: Lac repressor protein (represses TetR)
- **TetR**: Tet repressor protein (represses CI)
- **CI**: Lambda CI repressor protein (represses LacI)

#### Promoter Modules
- **P0**: LacI promoter (repressed by CI)
- **P1**: TetR promoter (repressed by LacI)
- **P2**: CI promoter (repressed by TetR)

Each promoter has multiple states representing different binding configurations with RNAP and repressors.

### Key Parameters

```
MAX_AMOUNT = 100     // Maximum protein molecules per species

// Kinetic Parameters
kd = 0.0075         // Protein degradation rate
kr_f = 0.5          // Repressor binding (forward)
kr_r = 1.0          // Repressor binding (reverse)
ka_f = 0.0033       // Activator binding (forward)
ka_r = 1.0          // Activator binding (reverse)
ko_f = 0.033        // RNAP binding (forward)
ko_r = 1.0          // RNAP binding (reverse)
kao_f = 1.0         // Activated RNAP binding (forward)
kao_r = 1.0         // Activated RNAP binding (reverse)

// Transcription Parameters
nc = 2.0            // Cooperativity coefficient
nr = 30.0           // RNAP concentration
ko = 0.05           // Open complex formation rate
kb = 0.0001         // Basal transcription rate
ng = 2.0            // Gene copy number
np = 10.0           // Proteins produced per transcript
ka = 0.25           // Activated transcription rate
```

### Initial Conditions

All protein species start at zero:
```
[LacI, TetR, CI] = [0, 0, 0]
```

Promoters start with 2 copies each (diploid):
```
[P0, P1, P2] = [2, 2, 2]
```

## Biological Significance

### Synthetic Biology Milestone

The Repressilator demonstrated:
- **First synthetic oscillator**: Proof-of-concept for biological clocks
- **Design principles**: Negative feedback with delay
- **Engineering biology**: Rational circuit construction
- **Temporal control**: Programmable biological timing

### Circuit Properties

1. **Cyclic repression**: A→B→C→A inhibition chain
2. **Negative feedback**: Each protein represses the next
3. **Time delays**: mRNA and protein synthesis create phase shifts
4. **Nonlinearity**: Cooperative binding (Hill coefficient = 2)

### Applications

- **Biological timers**: Controlled timing of cellular processes
- **Pulsatile production**: Rhythmic protein expression
- **Synchronization**: Population-level coordination
- **Drug delivery**: Timed therapeutic release
- **Research tool**: Understanding biological rhythms

## System Dynamics

### Oscillatory Behavior

The system exhibits:
- **Limit cycle oscillations**: Stable periodic orbits
- **Phase relationships**: 120° phase shift between proteins
- **Period control**: Tunable through parameter adjustment
- **Amplitude modulation**: Variable expression levels

### Temporal Characteristics

- **Period**: ~150-300 minutes (strain dependent)
- **Amplitude**: Varies with protein degradation rates
- **Phase coherence**: Maintained over multiple cycles
- **Noise sensitivity**: Affected by molecular fluctuations

### Dynamic Regimes

1. **Oscillatory**: Sustained rhythmic expression
2. **Damped**: Decaying oscillations to fixed point
3. **Bistable**: Switching between stable states
4. **Chaotic**: Irregular temporal patterns

## Mathematical Framework

### Deterministic Model

The core dynamics follow:
```
dLacI/dt = α₀/(1 + (CI/K₀)ⁿ) - γLacI
dTetR/dt = α₁/(1 + (LacI/K₁)ⁿ) - γTetR
dCI/dt = α₂/(1 + (TetR/K₂)ⁿ) - γCI
```

Where:
- α: Maximum transcription rates
- K: Repression constants
- n: Hill coefficients (cooperativity)
- γ: Degradation rates

### Oscillation Conditions

For sustained oscillations:
1. **Negative feedback**: Each gene represses the next
2. **Nonlinearity**: Hill coefficient > 1
3. **Time delays**: Finite transcription/translation times
4. **Strong repression**: High repression strength

### Stochastic Effects

- **Intrinsic noise**: From low molecule numbers
- **Extrinsic noise**: Cell-to-cell variability
- **Phase drift**: Random walk in oscillator phase
- **Coherence time**: Duration of synchronized oscillations

## Configuration Files

Several analysis scenarios are provided:

- `high_laci.json`: Peak LacI expression (LacI ≥ 80)
- `high_tetr.json`: Peak TetR expression (TetR ≥ 80)
- `high_ci.json`: Peak CI expression (CI ≥ 80)
- `oscillation_amplitude.json`: Extreme amplitude (LacI ≥ 90)
- `synchronization.json`: Early CI synchronization

## Usage Example

From the `src/` directory:

```bash
# Study oscillation peaks
python run_dwssa.py ../models/repressilator/high_laci.json

# Analyze amplitude variations
python run_learning_is.py ../models/repressilator/oscillation_amplitude.json

# Investigate phase dynamics
python run_guidedwssa.py ../models/repressilator/synchronization.json

# Examine protein expression bursts
python run_dwssa_plus_plus.py ../models/repressilator/high_tetr.json
```

## Rare Events of Interest

1. **Extreme amplitudes**: Unusually high protein peaks
2. **Phase jumps**: Sudden changes in oscillation phase
3. **Period doubling**: Bifurcation to longer periods
4. **Synchronization loss**: Breakdown of coherent oscillations
5. **Stochastic resonance**: Noise-enhanced oscillations

## Model Files

- **PRISM model**: `repressilator.sm`
- **Original source**: [GitHub - fluentverification/CaseStudies_StochasticModelChecking](https://github.com/fluentverification/CaseStudies_StochasticModelChecking/tree/main/GeneticCircuits/Repressilator/Default)
- **Configuration**: Default parameters

## Design Principles

### Requirements for Oscillations

1. **Ring topology**: Circular repression chain
2. **Odd number of nodes**: Prevents stable fixed points
3. **Strong repression**: Cooperative binding (n > 1)
4. **Balanced kinetics**: Comparable timescales
5. **Sufficient nonlinearity**: Sharp activation/repression

### Parameter Sensitivity

- **Degradation rates**: Control oscillation period
- **Repression strength**: Affects amplitude and robustness
- **Hill coefficients**: Determine nonlinearity
- **Copy number**: Influences noise levels

## Experimental Validation

### Key Experiments
- **Time-lapse microscopy**: Single-cell oscillation tracking
- **Flow cytometry**: Population-level dynamics
- **Fluorescent reporters**: Real-time protein levels
- **Perturbation studies**: Response to external stimuli

### Observed Phenomena
- **Heterogeneous periods**: Cell-to-cell variation
- **Phase relationships**: Predicted 120° shifts
- **Noise effects**: Stochastic phase diffusion
- **Growth coupling**: Cell cycle interactions

## Network Motifs

### Oscillator Designs
- **Repressilator**: 3-node negative feedback
- **Activator-repressor**: 2-node mixed feedback
- **Delayed negative feedback**: Single-node with delay
- **Coupled oscillators**: Multiple interacting clocks

### Comparison with Natural Clocks
- **Circadian rhythms**: 24-hour environmental cycles
- **Cell cycle**: DNA replication timing
- **Metabolic oscillations**: Glycolysis dynamics
- **Calcium oscillations**: Signaling rhythms

## Engineering Considerations

### Design Challenges
1. **Parameter tuning**: Achieving sustained oscillations
2. **Noise management**: Maintaining coherence
3. **Resource loading**: Cellular burden effects
4. **Context dependence**: Host strain variations

### Optimization Strategies
- **Modular design**: Orthogonal regulatory components
- **Protein engineering**: Improved degradation tags
- **Promoter tuning**: Balanced expression levels
- **Insulation**: Reducing crosstalk

## Applications and Extensions

### Synthetic Biology Applications
- **Biological pacemakers**: Cardiac rhythm regulation
- **Metabolic cycling**: Periodic pathway activation
- **Population dynamics**: Synchronized cellular behavior
- **Therapeutic timing**: Chronotherapy delivery

### Circuit Extensions
- **Coupled repressilators**: Multi-oscillator systems
- **Tunable periods**: Parameter-controlled timing
- **Gated oscillators**: Conditional activation
- **Memory circuits**: State-dependent oscillations

## Educational Value

The Repressilator teaches:
- **Systems biology**: Network dynamics and emergence
- **Synthetic biology**: Rational circuit design
- **Dynamical systems**: Oscillations and bifurcations
- **Stochastic processes**: Noise in biological systems

## Historical Impact

- **2000**: Original Nature publication by Elowitz & Leibler
- **Paradigm shift**: From parts to systems in biology
- **Field foundation**: Launched synthetic biology
- **Tool development**: Inspired oscillator design principles

## Mathematical Analysis

### Bifurcation Theory
- **Hopf bifurcations**: Onset of oscillations
- **Parameter space**: Oscillatory vs. non-oscillatory regions
- **Stability analysis**: Linear stability of fixed points
- **Period-doubling cascades**: Routes to chaos

### Stochastic Analysis
- **Master equation**: Full probabilistic description
- **Fokker-Planck**: Continuous approximation
- **Chemical Langevin**: Stochastic differential equations
- **Gillespie algorithm**: Exact stochastic simulation

## References

- Elowitz MB, Leibler S (2000). "A synthetic oscillatory network of transcriptional regulators." Nature 403(6767):335-8
- Stricker J et al. (2008). "A fast, robust and tunable synthetic gene oscillator." Nature 456(7221):516-9
- Original model: SBML representation from course materials

## Notes

This model implements the complete Repressilator circuit including:
- All three repressor proteins and their interactions
- Detailed promoter binding states
- Realistic kinetic parameters
- Stochastic molecular dynamics

The implementation captures both deterministic oscillations and stochastic fluctuations, making it ideal for studying rare events in biological timing systems and understanding how molecular noise affects synthetic genetic circuits.