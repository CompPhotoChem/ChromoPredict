import fieser, fieserkuhn, woodwardfieser

# library of base chromophores
base_library = {
    "fieser": fieser_base,
    "fieser_kuhn": fieser_kuhn_base,
    "woodward_extended": woodward_base,
    "woodward": woodward_base,
}

sub_library = {
    "fieser": fieser_subs,
    "fieser_kuhn": "None",
    "woodward_extended": woodward_subs,
    "woodward": woodward_subs,
}

factor_library = {
    "fieser": factors,
    "fieser_kuhn": fieser_kuhn_factors,
    "woodward_extended": factors,
    "woodward": factors,
}

base_value_library = {
    "fieser": fieser_base_values,
    "fieser_kuhn": fieser_kuhn_base_values,
    "woodward_extended": woodward_base_values,
    "woodward": woodward_base_values,
}

sub_value_library = {
    "fieser": fieser_sub_values,
    "fieser_kuhn": "None",
    "woodward_extended": woodward_sub_values,
    "woodward": woodward_sub_values,
}

factor_value_library = {
    "fieser": fieser_factor_values,
    "fieser_kuhn": fieser_kuhn_factor_values,
    "woodward_extended": woodward_factor_values,
    "woodward": woodward_factor_values,
}
