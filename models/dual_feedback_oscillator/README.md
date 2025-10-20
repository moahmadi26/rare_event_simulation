# Dual Feedback Oscillator Model

## Description

The Dual Feedback Oscillator is an advanced synthetic genetic circuit that combines both positive and negative feedback loops to create robust, tunable oscillations. Unlike the simple negative feedback loops of the Repressilator, this design incorporates dual regulatory mechanisms that provide enhanced control over oscillatory dynamics, improved robustness against noise, and more precise temporal programming capabilities.

## Model Details

### System Architecture

The dual feedback oscillator employs a two-node design with mixed feedback:

#### Core Components
- **AraC**: Arabinose regulatory protein (dual regulator)
- **LacI**: Lac repressor protein (negative regulator)

#### Promoter Modules
- **P1**: LacI promoter with dual feedback control
- **P2**: AraC promoter with dual feedback control

Each promoter exhibits complex binding states involving both activators and repressors.

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

// Complex Formation
kc_f = 0.05         // Complex formation (forward)
kc_r = 1.0          // Complex formation (reverse)
kecdiff = 1.0       // Effector complex diffusion
kecd = 0.005        // Effector complex degradation
kmdiff_f = 1.0      // Membrane diffusion (forward)
kmdiff_r = 0.01     // Membrane diffusion (reverse)

// Transcription Parameters
nc = 2.0            // Cooperativity coefficient
nr = 100.0          // RNAP concentration (higher than simple circuits)
ko = 0.05           // Open complex formation rate
kb = 0.0001         // Basal transcription rate
ng = 20.0           // Gene copy number (amplified)
np = 10.0           // Proteins produced per transcript
ka = 0.25           // Activated transcription rate
```

### Initial Conditions

```
[AraC, LacI] = [0, 0]        // Proteins start at zero
[P1, P2] = [20, 20]          // High promoter copy number
```

## Biological Significance

### Advanced Oscillator Design

The dual feedback architecture provides:
- **Positive feedback**: Self-activation for signal amplification
- **Negative feedback**: Cross-repression for oscillatory dynamics
- **Mixed regulation**: Complex temporal control
- **Robustness**: Multiple regulatory mechanisms

### Regulatory Mechanisms

1. **AraC dual function**: Acts as both activator and repressor
2. **Temporal switching**: Context-dependent regulatory roles
3. **Cooperative binding**: Enhanced nonlinearity (nc = 2)
4. **Resource competition**: Shared RNAP pool dynamics

### Applications

- **Precision timing**: Programmable biological clocks
- **Signal processing**: Complex waveform generation
- **Metabolic control**: Coordinated pathway regulation
- **Cellular computation**: Advanced logical operations
- **Therapeutic pulsing**: Controlled drug release patterns

## System Dynamics

### Oscillatory Regimes

The system can exhibit multiple dynamic behaviors:

1. **Stable oscillations**: Robust periodic dynamics
2. **Damped oscillations**: Transient rhythmic behavior
3. **Bistable switching**: Toggle-like state changes
4. **Complex oscillations**: Multi-peak or irregular patterns

### Feedback Interactions

#### Positive Feedback Loop
- **Self-activation**: Enhanced protein production
- **Bistability**: Stable high/low expression states
- **Memory**: Maintenance of expression states
- **Hysteresis**: Path-dependent switching

#### Negative Feedback Loop
- **Cross-repression**: Mutual inhibition
- **Oscillatory dynamics**: Sustained rhythms
- **Phase opposition**: Anti-correlated expression
- **Temporal delays**: Creating phase shifts

### Dynamic Properties

- **Tunable period**: Adjustable through parameter modification
- **Variable amplitude**: Controllable oscillation strength
- **Phase relationships**: Programmable temporal patterns
- **Noise resilience**: Robust against molecular fluctuations

## Mathematical Framework

### Coupled Dynamics

The system follows coupled differential equations:
```
dAraC/dt = α₁ · f₁(AraC, LacI) - γ·AraC
dLacI/dt = α₂ · f₂(AraC, LacI) - γ·LacI
```

Where f₁ and f₂ are complex regulatory functions involving both positive and negative terms.

### Regulatory Functions

For dual feedback with cooperativity:
```
f₁(AraC, LacI) = (k_act·AraC^n)/(K_act + AraC^n) · (K_rep)/(K_rep + LacI^n) + k_basal

f₂(AraC, LacI) = (k_act·LacI^n)/(K_act + LacI^n) · (K_rep)/(K_rep + AraC^n) + k_basal
```

### Stability Analysis

The system exhibits rich dynamics with:
- **Multiple fixed points**: Stable and unstable equilibria
- **Limit cycles**: Stable oscillatory orbits
- **Bifurcations**: Parameter-dependent transitions
- **Basins of attraction**: State-dependent outcomes

## Configuration Files

Several analysis scenarios are provided:

- `high_arac.json`: Peak AraC expression (AraC ≥ 80)
- `high_laci.json`: Peak LacI expression (LacI ≥ 80)
- `extreme_oscillation.json`: Extreme amplitude (AraC ≥ 90)
- `phase_synchronization.json`: Early synchronization events

## Usage Example

From the `src/` directory:

```bash
# Study dual feedback dynamics
python run_dwssa.py ../models/dual_feedback_oscillator/high_arac.json

# Analyze extreme oscillations
python run_learning_is.py ../models/dual_feedback_oscillator/extreme_oscillation.json

# Investigate phase relationships
python run_guidedwssa.py ../models/dual_feedback_oscillator/phase_synchronization.json

# Examine amplitude control
python run_dwssa_plus_plus.py ../models/dual_feedback_oscillator/high_laci.json
```

## Rare Events of Interest

1. **Amplitude spikes**: Unusually high protein expression
2. **Phase locking**: Synchronized oscillatory behavior
3. **Regime transitions**: Switching between dynamic modes
4. **Stochastic resonance**: Noise-enhanced oscillations
5. **Period doubling**: Bifurcation to complex rhythms

## Model Files

- **PRISM model**: `dual_feedback_oscillator.sm`
- **Original source**: [GitHub - fluentverification/CaseStudies_StochasticModelChecking](https://github.com/fluentverification/CaseStudies_StochasticModelChecking/tree/main/GeneticCircuits/Dual_Feedback_Osscillator/Default)
- **Configuration**: Default parameters

## Design Principles

### Requirements for Dual Feedback Oscillations

1. **Balanced feedback**: Appropriate positive/negative strengths
2. **Temporal separation**: Different timescales for feedback loops
3. **Nonlinear regulation**: Cooperative binding (n > 1)
4. **Resource dynamics**: Finite RNAP pool effects
5. **Degradation balance**: Matching protein lifetimes

### Parameter Sensitivity

- **Feedback strengths**: Control oscillation amplitude
- **Cooperativity**: Affects nonlinearity and robustness
- **Copy numbers**: Influence noise and dynamics
- **RNAP levels**: Determine resource competition effects

## Comparison with Other Oscillators

| Property | Repressilator | Dual Feedback | Activator-Repressor |
|----------|---------------|---------------|---------------------|
| Feedback type | Negative only | Mixed | Mixed |
| Node count | 3 | 2 | 2 |
| Complexity | Simple | High | Moderate |
| Robustness | Moderate | High | Moderate |
| Tunability | Limited | High | Moderate |

## Advanced Features

### Resource Competition
- **Shared RNAP pool**: nr = 100 molecules
- **Transcriptional interference**: Competition for polymerase
- **Burden effects**: Resource depletion dynamics
- **Context dependence**: Host cell interactions

### Complex Formation
- **Protein complexes**: Multi-component assemblies
- **Membrane interactions**: Spatial localization effects
- **Diffusion dynamics**: Transport processes
- **Cooperative assembly**: Synergistic interactions

## Experimental Considerations

### Implementation Challenges
1. **Dual regulation**: Engineering bifunctional proteins
2. **Parameter tuning**: Balancing feedback strengths
3. **Noise management**: Maintaining oscillation coherence
4. **Resource burden**: Minimizing cellular stress

### Characterization Methods
- **Fluorescent reporters**: Real-time protein dynamics
- **Flow cytometry**: Population heterogeneity analysis
- **Microfluidics**: Single-cell time-lapse studies
- **Optogenetics**: External control and perturbation

## Applications and Extensions

### Synthetic Biology Applications
- **Biological pacemakers**: Advanced timing circuits
- **Metabolic oscillators**: Rhythmic pathway control
- **Cell cycle regulation**: Artificial division timing
- **Circadian engineering**: Synthetic biological clocks

### Circuit Extensions
- **Multi-node networks**: Complex oscillatory networks
- **Coupled oscillators**: Synchronization studies
- **Hierarchical control**: Multi-level timing systems
- **Adaptive oscillators**: Environment-responsive rhythms

## Mathematical Analysis

### Bifurcation Theory
- **Hopf bifurcations**: Oscillation onset/termination
- **Saddle-node bifurcations**: Bistability regions
- **Period-doubling cascades**: Routes to complex dynamics
- **Heteroclinic cycles**: Slow-fast dynamics

### Stochastic Analysis
- **Noise-induced oscillations**: Stochastic switching
- **Coherence resonance**: Optimal noise levels
- **Phase diffusion**: Stochastic phase dynamics
- **Rare event statistics**: Extreme amplitude events

## Educational Value

The Dual Feedback Oscillator teaches:
- **Advanced circuit design**: Multi-feedback architectures
- **Nonlinear dynamics**: Complex temporal behaviors
- **Systems engineering**: Robust biological timing
- **Synthetic biology**: Next-generation genetic circuits

## Historical Context

- **Evolution from simple oscillators**: Building on Repressilator success
- **Mixed feedback theory**: Integration of positive/negative regulation
- **Robust design principles**: Engineering reliable biological timing
- **Advanced synthetic biology**: Second-generation genetic circuits

## References

- Stricker J et al. (2008). "A fast, robust and tunable synthetic gene oscillator." Nature 456(7221):516-9
- Atkinson MR et al. (2003). "Development of genetic circuitry exhibiting toggle switch or oscillatory behavior in Escherichia coli." Cell 113(5):597-607
- Original model: SBML implementation from course materials

## Notes

This model represents an advanced oscillator design that incorporates:
- Sophisticated regulatory mechanisms
- Resource competition effects
- Complex formation dynamics
- Enhanced parameter control

The dual feedback architecture provides superior performance compared to simple negative feedback oscillators, offering improved robustness, tunability, and dynamic range for synthetic biology applications.