
factors = {
    "homoanular_cyclodiene": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1",
    "conjugated_double_bond": "[#6]=[#6]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

# Woodward values
woodward_factor_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "conjugated_double_bond": 30,
    "exocyclic_double_bond": 5,
    "homoanular_cyclodiene": 39
}

# Fieser-Kuhn
fieser_kuhn_factors = {
    "alkyl" : "[#6;X4]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

fieser_kuhn_factor_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "endocyclic_double_bond": 1,
    "alkyl" : 1
}


fieser_factor_values = { # doi: 10.1021/jo01164a003
    "endocyclic_double_bond" : 30,
    "conjugated_double_bond": 30,
    "exocyclic_double_bond": 5,
    "homoanular_cyclodiene": 39
}


fieser_kuhn_factors = {
    "alkyl" : "[#6;X4]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}
