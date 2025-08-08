# Datasets

## Enones

| Compound Class       | Core Motif Type                                | Substitution Positions               | Substitution Types                                                                                       | No. of Molecules | Stereoisomers Considered                     | Notes                                                                                                          |
|----------------------|------------------------------------------------|---------------------------------------|----------------------------------------------------------------------------------------------------------|------------------|-----------------------------------------------|----------------------------------------------------------------------------------------------------------------|
| Acyclic ketones      | Linear α,β-unsaturated methyl ketones          | α, β, β′                              | Methyl (C), Methoxy (OC), Chloro (Cl), Bromo (Br), Hydroxy (O), Hydrogen (H)                             | 216              | cis/trans for β-, α,β-, β,β′ substitutions    | All mono-, di-, and tri-substitution patterns systematically enumerated                                       |
| Acyclic aldehydes    | Linear α,β-unsaturated aldehydes               | α, β, β′                              | Same as above                                                                                            | 216              | cis/trans for β-, α,β-, β,β′ substitutions    | Same combinatorial approach as ketones                                                                         |
| Acyclic acids        | Linear α,β-unsaturated carboxylic acids        | α, β, β′                              | Same as above                                                                                            | 216              | cis/trans for β-, α,β-, β,β′ substitutions    | Same combinatorial approach as ketones                                                                         |
| Cyclopentenones      | α,β-unsaturated cyclopentenone derivatives     | α, β                                  | Same as above                                                                                            | 36               | —                                             | Mono- and di-substitution patterns; only one β-position available                                              |
| Cyclohexenones       | α,β-unsaturated cyclohexenone derivatives      | α, β                                  | Same as above                                                                                            | 36               | —                                             | Same as cyclopentenones                                                                                         |
| **Total**            | —                                              | —                                     | —                                                                                                        | **720**          | —                                             | All structures encoded as isomeric SMILES, 3D geometries generated with ETKDG, optimized with DFT, TD-DFT data |

### Computational Details

All calculations were performed with [VeloxChem]().

- **Geometry Optimization:** DFT (B3LYP/def2-TZVP/D4)  
- **Excited States:** TD-DFT (VeloxChem), 10 lowest singlet states computed  
- **State Assignment:** S₁ = nπ*, S₂ = ππ*  
- **Analysis Focus:** Vertical excitation energy for S₀ → S₂ (ππ*) transitions  
- **Additional Data:** First five excitations + oscillator strengths available in dataset  

