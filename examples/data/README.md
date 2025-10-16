# Datasets

## Woodward Compounds

| Compound Class       | Core Motif Type                                | Substitution Positions               | Substitution Types                                                                                       | No. of Molecules | Stereoisomers Considered                     | Notes                                                                                                          |
|----------------------|------------------------------------------------|---------------------------------------|----------------------------------------------------------------------------------------------------------|------------------|-----------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| Acyclic ketones      | Linear α,β-unsaturated methyl ketones          | α, β, β′                              | Methyl (C), Methoxy (OC), Chloro (Cl), Bromo (Br), Hydroxy (O), Hydrogen (H)                             | 216              | cis/trans for β-, α,β-, β,β′ substitutions    | All mono-, di-, and tri-substitution patterns systematically enumerated                                       |
| Acyclic aldehydes    | Linear α,β-unsaturated aldehydes               | α, β, β′                              | Same as above                                                                                            | 216              | cis/trans for β-, α,β-, β,β′ substitutions    | Same combinatorial approach as ketones                                                                         |
| Acyclic acids        | Linear α,β-unsaturated carboxylic acids        | α, β, β′                              | Same as above                                                                                            | 216              | cis/trans for β-, α,β-, β,β′ substitutions    | Same combinatorial approach as ketones                                                                         |
| Cyclopentenones      | α,β-unsaturated cyclopentenone derivatives     | α, β                                  | Same as above                                                                                            | 36               | —                                             | Mono- and di-substitution patterns; only one β-position available                                              |
| Cyclohexenones       | α,β-unsaturated cyclohexenone derivatives      | α, β                                  | Same as above                                                                                            | 36               | —                                             | Same as cyclopentenones                                                                                         |
| **Total**            | —                                              | —                                     | —                                                                                                        | **720**          | —                                             | All structures encoded as isomeric SMILES, 3D geometries generated with ETKDG, optimized with DFT, TD-DFT data |

## List of Datasets

- [Inference Set (E01-E26): Woodward and *ab initio* (B3LYP/TD-B3LYP, xTB/TD-B3LYP, xTB/HF-TDA, xTB/B3LYP-TDA) predicted maxima](https://github.com/CompPhotoChem/ChromoPredict/blob/main/examples/data/ds_enones_exp_vs_calc_1state.csv)
- [Training Set (Woodward compounds): Woodward and *ab initio* (B3LYP/TD-B3LYP) predicted maxima](https://github.com/CompPhotoChem/ChromoPredict/blob/main/examples/data/ds_enones_b3lyp_woodward_1state.csv)
- [Training Set (Woodward compounds): Woodward and *ab initio* (B3LYP/TD-B3LYP) predicted maxima (seperated by stereochemistry)](https://github.com/CompPhotoChem/ChromoPredict/blob/main/examples/data/data_720_wf_wfr_wfrstereo.csv)
- [Training Set (Woodward compounds): *Ab initio* (B3LYP/TD-B3LYP) predicted maxima (all 10 predicted states)](https://github.com/CompPhotoChem/ChromoPredict/blob/main/examples/data/ds_enones_b3lyp_10states.csv)

## Computational Details

All calculations were performed with [VeloxChem](https://veloxchem.org/docs/intro.html).

- **Geometry Optimization:** DFT (B3LYP/def2-TZVP/D4) or xTB
- **Excited States:** linear response TD-DFT (RPA) or Tamm–Dancoff Approximation (TDA); 10 lowest singlet states computed  
- **State Assignment:** S₁ = nπ*, S₂ = ππ*
- **Analysis Focus:** Vertical excitation energy for S₀ → S₂ (ππ*) transitions
