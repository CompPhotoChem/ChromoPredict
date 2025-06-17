# 

# fieser base
fieser_base = {
    "homoanular_cyclic": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1", #[#6;r]=[#6;r]-[#6;r]=[#6;r]
    "acyclic_diene": "[#6;!R]=[#6;!R]-[#6;!R]=[#6;!R]",
    "heteroanular_cyclic": "[#6]=[#6;R]-[#6;R]=[#6]"
}

fieser_base_values = { # doi: 10.1021/jo01164a003
    "acyclic_diene": 217,
    "homoanular_cyclic": 253,
    "heteroanular_cyclic": 214
}

# increments for additional double bonds
factors = {
    "homoanular_cyclodiene": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1",
    "conjugated_double_bond": "[#6]=[#6]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

fieser_factor_values = { # doi: 10.1021/jo01164a003
    "endocyclic_double_bond" : 30,
    "conjugated_double_bond": 30,
    "exocyclic_double_bond": 5,
    "homoanular_cyclodiene": 39
}

# increments from further substituents
fieser_subs = {
    "dialkylamine": "[#7](-[#6])(-[#6])",
    "thioether": "[#16]-[#6]",
    "phenoxy": "[#8]-c1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
    "halogen": "[F,Cl,Br,I]",
    "alkoxy": "[#8]-[#6]",
    "alkyl": "[#6;X4]",
    "carboxy": "[#6](=[#8])-[#8H]",
}

fieser_sub_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "alkoxy": 6,
    "alkyl": 5,
    "carboxy": 0,
    "dialkylamine": 60,
    "halogen": 10,
    "phenoxy": 18,
    "thioether": 30
}
