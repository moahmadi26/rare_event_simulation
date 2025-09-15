# Motility Regulation Model

## Description

The Motility Regulation model represents a gene regulatory network controlling bacterial motility. This system models the complex interactions between transcription factors, genes, and proteins involved in bacterial flagellar synthesis and movement control. The network consists of 9 chemical species interacting through 12 reaction channels.

## Model Details

### Chemical Species
- **codY**: Gene encoding the CodY transcription factor
- **CodY**: CodY protein (transcriptional repressor)
- **flache**: Gene encoding sigma factor SigD
- **SigD**: Sigma factor D protein (activator of flagellar genes)
- **SigD_hag**: SigD bound to hag promoter
- **hag**: Gene encoding flagellin
- **Hag**: Flagellin protein (flagellar filament component)
- **CodY_flache**: CodY bound to flache promoter
- **CodY_hag**: CodY bound to hag promoter

### Reactions

The model is defined by the following set of reactions:

```
R₁:  codY → codY + CodY           (rate: 0.1)
R₂:  CodY → ∅                     (rate: 0.0002)
R₃:  flache → flache + SigD       (rate: 1.0)
R₄:  SigD → ∅                     (rate: 0.0002)
R₅:  SigD_hag → SigD + hag + Hag  (rate: 1.0)
R₆:  Hag → ∅                      (rate: 0.0002)
R₇:  SigD + hag → SigD_hag        (rate: 0.01)
R₈:  SigD_hag → SigD + hag        (rate: 0.1)
R₉:  CodY + flache → CodY_flache  (rate: 0.02)
R₁₀: CodY_flache → CodY + flache  (rate: 0.1)
R₁₁: CodY + hag → CodY_hag        (rate: 0.01)
R₁₂: CodY_hag → CodY + hag        (rate: 0.1)
```

### Initial Conditions

Default initial populations:
```
[codY, CodY, flache, SigD, SigD_hag, hag, Hag, CodY_flache, CodY_hag]
= [1, 10, 1, 10, 1, 1, 10, 1, 1]
```

## Biological Significance

This model captures the regulatory dynamics of bacterial motility:

- **CodY regulation**: Acts as a global regulator responding to nutrient availability
- **SigD cascade**: Controls the expression of late flagellar genes
- **Flagellin synthesis**: Essential for bacterial movement and chemotaxis
- **Feedback loops**: Complex regulatory interactions ensuring proper motility control

The system exhibits:
- **Bistability**: Switching between motile and non-motile states
- **Noise-driven transitions**: Stochastic switching in gene expression
- **Hierarchical regulation**: Coordinated control of flagellar assembly

## Configuration Files

Available configuration:
- `motility_regulation.json`: Standard configuration for rare event analysis

## Usage Example

From the `src/` directory:

```bash
# Run with dwSSA
python run_dwssa.py ../models/motility/motility_regulation.json

# Run with dwSSA++
python run_dwssa_plus_plus.py ../models/motility/motility_regulation.json

# Run with learning-based IS
python run_learning_is.py ../models/motility/motility_regulation.json

# Run with guided wSSA
python run_guidedwssa.py ../models/motility/motility_regulation.json
```

## Rare Events of Interest

1. **High CodY levels**: Excessive repressor accumulation
2. **Complete motility shutdown**: All flagellar genes repressed
3. **Maximum Hag production**: Peak flagellin synthesis
4. **Synchronized gene activation**: Coordinated expression bursts

## Model File

The PRISM model is defined in `motility_unb.sm`

## Biological Context

Bacterial motility regulation is crucial for:
- **Pathogenesis**: Movement toward host cells
- **Biofilm formation**: Transition between planktonic and sessile states
- **Chemotaxis**: Response to chemical gradients
- **Environmental adaptation**: Survival in changing conditions

## References

This model is based on studies of bacterial flagellar regulation, particularly in *Bacillus subtilis* and related species.