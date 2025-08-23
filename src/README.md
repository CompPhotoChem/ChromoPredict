# ChromoPredict: Rule-Based UV–Vis Prediction

**ChromoPredict** implements digitized Woodward–Fieser, Fieser, and Fieser–Kuhn rules for predicting ππ* absorption maxima (λ<sub>max</sub>) from SMILES.

## Workflow

1. **Rule Set Selection**  
   The input SMILES is matched against chromophore patterns to identify the correct rule set. Users can select original, extended, or refined rules via `chromlib`.

2. **Base Value Assignment**  
   Each chromophore class carries a characteristic base absorption, automatically assigned by the corresponding module (`woodwardfieser`, `fieser`, `fieserkuhn`).

3. **Structural Features**  
   Conjugation extensions, exocyclic double bonds, and ring closures are identified and tagged to adjust λ_max predictions.

4. **Substituent Increments**  
   Type- and position-dependent increments account for electronic effects of substituents, applied to the chromophore and its extensions.

5. **Optional Solvent Corrections**  
   Empirical offsets can correct for solvent effects (e.g., hypsochromic shifts in polar solvents).

Predictions are performed via `cp.predict(smiles="C=CC(=O)C", chromlib="woodward_extended")`, returning λ<sub>max</sub> and optionally detailed contributions for 
interpretable analysis (`verbose=True`) and images of the structure with highlighted structural features (`draw=True`).
