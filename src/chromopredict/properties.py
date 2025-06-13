
woodward_base_values = { # doi: textbook
    "alpha_beta_aldehyd": 210,
    "alpha_beta_ketone": 215,
    "alpha_beta_ester": 195,
    "cyclopentenone": 202,
    "cyclohexanone": 215, # doi: 10.1021/jo01164a003
    "alpha_beta_acid" : 195, # doi: 
}

fieser_base_values = { # doi: 10.1021/jo01164a003
    "acyclic_diene": 217,
    "homoanular_cyclic": 253,
    "heteroanular_cyclic": 214
}

fieser_kuhn_base_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "polyene_4_conj": 114
}

woodward_factor_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "conjugated_double_bond": 30,
    "exocyclic_double_bond": 5,
    "homoanular_cyclodiene": 39
}

fieser_factor_values = { # doi: 10.1021/jo01164a003
    "endocyclic_double_bond" : 30,
    "conjugated_double_bond": 30,
    "exocyclic_double_bond": 5,
    "homoanular_cyclodiene": 39
}

fieser_kuhn_factor_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "endocyclic_double_bond": 1,
    "alkyl" : 1
}

woodward_sub_values = { # doi: 10.1021/jo01164a003, doi: 10.1021/jo01085a617
    "alpha": {"alkoxy": 35, "hydroxy": 35, "alkyl": 10, "bromo": 5, "chloro": 5, "carboxy": 0},
    "beta": {"alkoxy": 30, "hydroxy": 30, "alkyl": 12, "bromo": 10, "chloro": 5, "carboxy": 0},
    "gamma": {"hydroxy": 50, "bromo": 30, "alkyl": 18, "alkoxy": 17, "carboxy": 12, "chloro": 12},
    "higher": {"bromo": 30, "alkyl": 18, "carboxy": 12, "chloro": 12, "alkoxy": 0, "hydroxy": 0}
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
