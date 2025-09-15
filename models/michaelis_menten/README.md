# Michaelis-Menten Model

## Description

The Michaelis-Menten model is a fundamental enzyme kinetics model that describes how enzymes catalyze the conversion of substrates into products. This classic model forms the basis for understanding enzyme kinetics and is one of the most important models in biochemistry. It consists of 4 chemical species interacting through 3 reaction channels.

## Model Details

### Chemical Species
- **E**: Free enzyme
- **S**: Substrate
- **C**: Enzyme-substrate complex
- **P**: Product

### Reactions

The model follows the classic Michaelis-Menten mechanism:

```
R₁: E + S ⇌ C     (binding/unbinding)
R₂: C → E + P     (catalysis)
```

Implemented as three elementary reactions:
```
R₁: E + S → C     (rate: k₁)
R₂: C → E + S     (rate: k₋₁)
R₃: C → E + P     (rate: k₂)
```

### Initial Conditions

Default initial populations:
```
[E, S, C, P] = [100, 100, 0, 0]
```

## Kinetic Parameters

The model behavior is characterized by:
- **k₁**: Association rate constant (enzyme-substrate binding)
- **k₋₁**: Dissociation rate constant (complex breakdown)
- **k₂**: Catalytic rate constant (product formation)
- **Kₘ = (k₋₁ + k₂)/k₁**: Michaelis constant
- **Vₘₐₓ = k₂[E]₀**: Maximum reaction velocity

## Biological Significance

The Michaelis-Menten model is foundational for:
- **Drug design**: Understanding drug-enzyme interactions
- **Metabolic engineering**: Optimizing enzymatic pathways
- **Clinical diagnostics**: Enzyme assays and biomarkers
- **Biotechnology**: Industrial enzyme applications

### Key Assumptions
1. Enzyme concentration is much lower than substrate
2. Product does not significantly reverse to substrate
3. Steady-state approximation for complex formation
4. Single substrate, single product reaction

## Configuration Files

Available configuration:
- `high_product.json`: Rare event of high product formation (P ≥ 80)

## Usage Example

From the `src/` directory:

```bash
# Run with dwSSA
python run_dwssa.py ../models/michaelis_menten/high_product.json

# Run with dwSSA++
python run_dwssa_plus_plus.py ../models/michaelis_menten/high_product.json

# Run with learning-based IS
python run_learning_is.py ../models/michaelis_menten/high_product.json

# Run with guided wSSA
python run_guidedwssa.py ../models/michaelis_menten/high_product.json
```

## Rare Events of Interest

1. **High product formation**: P ≥ 80 molecules at early times
2. **Complete substrate depletion**: S = 0
3. **Maximum complex formation**: High C accumulation
4. **Enzyme saturation**: All enzyme bound in complex

## Model File

The PRISM model is defined in `michaelis_menten.sm`

## Stochastic Behavior

At the molecular level, the system exhibits:
- **Stochastic bursts**: Random fluctuations in product formation
- **Substrate depletion effects**: Slowdown as substrate decreases
- **Complex lifetime variations**: Stochastic binding/unbinding events
- **Single-molecule events**: Important at low concentrations

## Mathematical Properties

### Deterministic Limit
In the large-number limit, follows the Michaelis-Menten equation:
```
v = Vₘₐₓ[S]/(Kₘ + [S])
```

### Stochastic Features
- **Coefficient of variation**: Increases at low molecule numbers
- **First passage times**: Time to reach product thresholds
- **Burst statistics**: Distribution of product formation events
- **Quasi-steady state**: Complex concentration equilibration

## Applications

The model is used to study:
- **Enzyme efficiency**: Catalytic effectiveness
- **Substrate specificity**: Selective catalysis
- **Competitive inhibition**: Drug mechanism of action
- **Allosteric regulation**: Cooperative enzyme behavior
- **Metabolic flux**: Pathway bottlenecks

## Extensions

Common extensions include:
- **Competitive inhibition**: Additional inhibitor species
- **Cooperative binding**: Hill-type kinetics
- **Product inhibition**: Reverse reactions
- **Multi-substrate reactions**: Complex mechanisms
- **Compartmentalization**: Spatial effects

## Historical Note

Proposed by Leonor Michaelis and Maud Menten in 1913, this model revolutionized enzyme biochemistry and remains the standard framework for enzyme kinetics analysis over a century later.

## References

This implementation follows the standard Michaelis-Menten formulation widely used in systems biology and biochemistry textbooks.