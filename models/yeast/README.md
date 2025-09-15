# Yeast Polarization Model

## Description

The Yeast Polarization model represents the G-protein signaling cascade involved in yeast cell polarization during mating. This system models how yeast cells establish and maintain polarity in response to pheromone gradients, a fundamental process in cell biology. The network consists of 7 chemical species interacting through 8 reaction channels.

## Model Details

### Chemical Species
- **R**: Receptor (unbound)
- **L**: Ligand (pheromone, e.g., α-factor)
- **RL**: Receptor-ligand complex
- **G**: Heterotrimeric G-protein (inactive)
- **Gₐ**: G-protein α subunit (active)
- **G_bg**: G-protein βγ subunit complex
- **Gd**: G-protein α subunit (deactivated)

### Reactions

The model is defined by the following set of reactions:

```
R₁: ∅ → R                     (rate: 0.0038)
R₂: R → ∅                     (rate: 4.00×10⁻⁴)
R₃: L + R → RL + L            (rate: 0.042)
R₄: RL → R                    (rate: 0.0100)
R₅: RL + G → Gₐ + G_bg        (rate: 0.011)
R₆: Gₐ → Gd                   (rate: 0.100)
R₇: Gd + G_bg → G             (rate: 1.05×10³)
R₈: ∅ → RL                    (rate: 3.21)
```

### Initial Conditions

Default initial populations:
```
[R, L, RL, G, Gₐ, G_bg, Gd] = [50, 2, 0, 50, 0, 0, 0]
```

## Biological Significance

The yeast polarization system is a paradigm for understanding:

- **Cell polarity**: How cells establish spatial organization
- **Signal transduction**: G-protein coupled receptor (GPCR) signaling
- **Gradient sensing**: Detection and response to chemical gradients
- **Mating response**: Cellular preparation for fusion

Key features of the system:
- **Signal amplification**: Small pheromone concentrations trigger robust responses
- **Adaptation**: Cells adjust sensitivity based on signal strength
- **Noise filtering**: Distinguishing true gradients from random fluctuations
- **Polarization commitment**: Establishment of a stable mating projection

## Configuration Files

Available configurations:
- `yeast_polarization.json`: Standard configuration
- `yeast_polarization_new.json`: Alternative parameter set
- `high_gbg.json`: High G-protein βγ subunit scenario

## Usage Example

From the `src/` directory:

```bash
# Run with dwSSA
python run_dwssa.py ../models/yeast/yeast_polarization.json

# Run with dwSSA++
python run_dwssa_plus_plus.py ../models/yeast/yeast_polarization_new.json

# Run with learning-based IS
python run_learning_is.py ../models/yeast/high_gbg.json

# Run with guided wSSA
python run_guidedwssa.py ../models/yeast/yeast_polarization.json
```

## Rare Events of Interest

1. **High G-protein activation**: Excessive Gₐ accumulation
2. **Receptor saturation**: All receptors bound to ligand
3. **Complete signal shutdown**: No active G-proteins
4. **Spontaneous polarization**: Activation without ligand
5. **G_bg accumulation**: High levels of free βγ subunits

## Model File

The PRISM model is defined in `yeast_unb.sm`

## Biological Context

Yeast polarization is essential for:
- **Mating**: Formation of mating projections (shmoos)
- **Budding**: Asymmetric cell division
- **Gradient tracking**: Following pheromone sources
- **Cell-cell communication**: Coordinated mating responses

The model captures the core GPCR signaling module conserved across eukaryotes, making it relevant for understanding:
- Human GPCR signaling (>30% of drugs target GPCRs)
- Cell migration in development and cancer
- Immune cell chemotaxis
- Neuronal guidance

## Mathematical Features

The system exhibits:
- **Ultrasensitivity**: Sharp response to ligand concentration
- **Stochastic focusing**: Noise-induced symmetry breaking
- **Temporal dynamics**: Transient vs. sustained responses
- **Spatial coupling**: When extended to include spatial components

## References

This model is based on the well-studied *Saccharomyces cerevisiae* mating response pathway, particularly the work on G-protein dynamics and polarization establishment.