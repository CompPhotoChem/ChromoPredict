general_rules = {
    "woodward_extended": "[#6]=[#6]-[#6](=[#8])-[#6]",
    "woodward": "[#6]=[#6]-[#6]=[#8]",
    "fieser_kuhn": "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]",
    "fieser": "[#6]=[#6]-[#6]=[#6]"
}

# Woodward
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

# Fieser
fieser_base = {
    "homoanular_cyclic": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1",
    "acyclic_diene": "[#6;!R]=[#6;!R]-[#6;!R]=[#6;!R]",
    "heteroanular_cyclic": "[#6]=[#6;R]-[#6;R]=[#6]"
}

fieser_base_values = { # doi: 10.1021/jo01164a003
    "acyclic_diene": 217,
    "homoanular_cyclic": 253,
    "heteroanular_cyclic": 214
}

# Fieser Kuhn
fieser_kuhn_base = {
    "polyene_4_conj": "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"
}

fieser_kuhn_base_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "polyene_4_conj": 114
}
