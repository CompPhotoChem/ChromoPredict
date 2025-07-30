# Properties for enones

# base chromophore units 
woodward_base = {
    "alpha_beta_aldehyd": "[#6]=[#6]-[#6;H1](=[#8])",
    "alpha_beta_ketone": "[#6]=[#6]-[#6](=[#8])-[#6]",
    "alpha_beta_ester": "[#6]=[#6]-[#6](=[#8])-[#8]-[#6]",
    "alpha_beta_acid": "[#6]=[#6][#6](=[#8])[O;H1,-1]",
    "cyclohexanone": "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1(=[#8])",
    "cyclopentenone": "[#6]1:[#6]:[#6]:[#6]:[#6]:1(=[#8])",
}

woodward_base_values = { # doi: textbook
    "alpha_beta_aldehyd": 210,
    "alpha_beta_ketone": 215,
    "alpha_beta_ester": 195,
    "cyclopentenone": 202,
    "cyclohexanone": 215, # doi: 10.1021/jo01164a003
    "alpha_beta_acid" : 195, # doi: 
}

woodward_refine_base_values = { # doi: textbook
    "alpha_beta_aldehyd": 215,
    "alpha_beta_ketone": 206,
    "alpha_beta_ester": 195, #not refined
    "cyclopentenone": 202, #not refined
    "cyclohexanone": 215, #not refined
    "alpha_beta_acid" : 192, 
}

# increments from additional double bonds
factors = {
    "homoanular_cyclodiene": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1",
    "conjugated_double_bond": "[#6]=[#6]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

woodward_factor_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "conjugated_double_bond": 30,
    "exocyclic_double_bond": 5,
    "homoanular_cyclodiene": 39
}

# increments from substituents
woodward_subs = {
    "alkoxy": "[#8]-[#6]",
    "alkyl": "[#6;X4]",
    "bromo": "[Br]",
    "carboxy": "[#6](=[#8])-[#8H]",
    "chloro": "[Cl]",
    "hydroxy": "[#8H]"
}

woodward_sub_values = { # doi: 10.1021/jo01164a003, doi: 10.1021/jo01085a617
    "alpha":  {"alkoxy": 35, "hydroxy": 35, "alkyl": 10, "bromo": 25, "chloro": 15, "carboxy":  0},
    "beta":   {"alkoxy": 30, "hydroxy": 30, "alkyl": 12, "bromo": 30, "chloro": 12, "carboxy":  0},
    "gamma":  {"alkoxy": 17, "hydroxy": 50, "alkyl": 18, "bromo":  0, "chloro": 12, "carboxy": 12},
    "higher": {"alkoxy":  0, "hydroxy":  0, "alkyl": 18, "bromo":  0, "chloro": 12, "carboxy": 12}
}

woodward_refine_sub_values = {
    "alpha":  {"alkoxy": 35, "hydroxy": 43, "alkyl": 14, "bromo": 40, "chloro": 29, "carboxy": 0},
    "beta":   {"alkoxy": 23, "hydroxy": 17, "alkyl": 17, "bromo": 34, "chloro": 24, "carboxy": 0},
    "gamma":  {"alkoxy": 17, "hydroxy": 50, "alkyl": 18, "bromo":  0, "chloro": 12, "carboxy": 12},
    "higher": {"alkoxy":  0, "hydroxy":  0, "alkyl": 18, "bromo":  0, "chloro": 12, "carboxy": 12}
}

# anchors for plotting and structure assignments
woodward_markers = {
    1: "alpha",
    2: "beta",
    3: "gamma",
    4: "higher"
}


