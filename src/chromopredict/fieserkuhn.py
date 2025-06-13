# Rules for polynes with more than four conj. double bonds

# base chromophore
fieser_kuhn_base = {
    "polyene_4_conj": "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"
}

fieser_kuhn_base_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "polyene_4_conj": 114
}

# increments further conjugated double bonds
fieser_kuhn_factors = {
    "alkyl" : "[#6;X4]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

fieser_kuhn_factor_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "endocyclic_double_bond": 1,
    "alkyl" : 1
}
