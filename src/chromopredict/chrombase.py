
general_rules = {
    "woodward_extended": "[#6]=[#6]-[#6](=[#8])-[#6]",
    "woodward": "[#6]=[#6]-[#6]=[#8]",
    "fieser_kuhn": "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]",
    "fieser": "[#6]=[#6]-[#6]=[#6]"
}

woodward_base = {
    "alpha_beta_aldehyd": "[#6]=[#6]-[#6;H1](=[#8])",
    "alpha_beta_ketone": "[#6]=[#6]-[#6](=[#8])-[#6]",
    "alpha_beta_ester": "[#6]=[#6]-[#6](=[#8])-[#8]-[#6]",
    "alpha_beta_acid": "[#6]=[#6][#6](=[#8])[O;H1,-1]",
    "cyclohexanone": "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1(=[#8])",
    "cyclopentenone": "[#6]1:[#6]:[#6]:[#6]:[#6]:1(=[#8])",
}

fieser_base = {
    "homoanular_cyclic": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1",
    "acyclic_diene": "[#6;!R]=[#6;!R]-[#6;!R]=[#6;!R]",
    "heteroanular_cyclic": "[#6]=[#6;R]-[#6;R]=[#6]"
}

fieser_kuhn_base = {
    "polyene_4_conj": "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"
}
