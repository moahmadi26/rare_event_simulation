# Enzymatic Futile Cycle Model

## Description

The Enzymatic Futile Cycle is a biochemical reaction network that models substrate cycling between phosphorylated and dephosphorylated states. This system consists of 6 chemical species interacting through 6 reaction channels, representing a fundamental regulatory mechanism in cellular metabolism.

## Model Details

### Chemical Species
- **S₁**: Enzyme 1
- **S₂**: Substrate (unphosphorylated form)
- **S₃**: Enzyme-substrate complex 1
- **S₄**: Enzyme 2  
- **S₅**: Substrate (phosphorylated form)
- **S₆**: Enzyme-substrate complex 2

### Reactions

The model is defined by the following set of reactions:

```
R₁: S₁ + S₂ → S₃       (rate: 1.0)
R₂: S₃ → S₁ + S₂       (rate: 1.0)
R₃: S₃ → S₁ + S₅       (rate: 0.1)
R₄: S₄ + S₅ → S₆       (rate: 1.0)
R₅: S₆ → S₄ + S₅       (rate: 1.0)
R₆: S₆ → S₄ + S₂       (rate: 0.1)
```

### Initial Conditions

Default initial populations:
```
[S₁, S₂, S₃, S₄, S₅, S₆] = [1, 50, 0, 1, 50, 0]
```

## Biological Significance

Futile cycles are important in:
- **Metabolic regulation**: Providing sensitive switching between metabolic states
- **Signal amplification**: Small changes in enzyme activity can produce large metabolic flux changes
- **Heat generation**: Energy dissipation in thermogenesis
- **Ultrasensitivity**: Creating sharp, switch-like responses to signals

## Configuration Files

Several pre-configured scenarios are available:

- `enzymatic_futile_cycle.json`: Standard configuration
- `enzymatic_high_s5.json`: High phosphorylated substrate scenario
- `enzymatic_futile_cycle_30.json`: S₅ ≥ 30 rare event
- `enzymatic_futile_cycle_35.json`: S₅ ≥ 35 rare event
- `enzymatic_futile_cycle_40.json`: S₅ ≥ 40 rare event
- `enzymatic_s5_60_time2.json`: S₅ ≥ 60 at time 2.0

## Usage Example

From the `src/` directory:

```bash
# Run with dwSSA
python run_dwssa.py ../models/enzym/enzymatic_futile_cycle.json

# Run with learning-based IS
python run_learning_is.py ../models/enzym/enzymatic_high_s5.json

# Run with guided wSSA
python run_guidedwssa.py ../models/enzym/enzymatic_futile_cycle_40.json
```

## Rare Events of Interest

1. **High phosphorylation**: Events where S₅ reaches unusually high levels
2. **Low substrate**: Depletion of unphosphorylated substrate S₂
3. **Complex accumulation**: Build-up of enzyme-substrate complexes

## Model File

The PRISM model is defined in `enzym_unb.sm`

## References

This model represents a classic enzymatic futile cycle studied in metabolic regulation and systems biology.