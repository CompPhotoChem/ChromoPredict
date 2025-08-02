import warnings
import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image

from chromopredict.strucfeatures import *
from chromopredict import *

def __woodward_calc(
        base,
        factor,
        subs,
        mol,
        solvent,
        verbose,
        mol_type='woodward'
        ):

    total = 0
    d_contrib = {}

    # base values
    total += check_values(base, mol_type, base_value_library)
    total += check_values(factor, mol_type, factor_value_library)

    # increments for substituents
    total += max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0)
    total += max([sub["value"] for sub in subs if sub["sub_type"] == "beta"], default=0)
    total += max([sub["value"] for sub in subs if sub["sub_type"] == "gamma"], default=0)
    total += max([sub["value"] for sub in subs if sub["sub_type"] == "higher"], default=0)
    
    # exocyclic bonds
    total += 5 * count_exo_bonds(mol)
    
    # solvent
    total += solvent_values.get(solvent, 0)

    if verbose:
        # base values
        d_contrib['base'] = check_values(base, mol_type, base_value_library)
        d_contrib["factor"] = check_values(factor, mol_type, factor_value_library)

        # increments for substituents
        d_contrib["alpha"], d_contrib['alpha_all']  = __get_max_sub(pos="alpha", subs=subs)
        d_contrib["beta"], d_contrib["beta_all"] = __get_max_sub(pos="beta", subs=subs)
        d_contrib["gamma"], d_contrib["gamma_all"] = __get_max_sub(pos="gamma", subs=subs)
        d_contrib["higher"], d_contrib["higher_all"] = __get_max_sub(pos="higher", subs=subs) 
        
        # exocyclic bonds
        d_contrib["exo"] =  5 * count_exo_bonds(mol)
        
        # solvent
        d_contrib["solvent"] = solvent_values.get(solvent, 0)

    return round(total), d_contrib


def __woodward_extended_calc(
        base,
        factor,
        subs,
        mol,
        solvent,
        verbose,
        mol_type='woodward_extended'
        ):

    total = 0
    d_contrib = {}

    # base value
    total += 212.82

    # substituent contributions
    total += 11.37 * len(subs)
    total += 4.79 * count_exo_bonds(mol)
    
    # solvent
    total += solvent_values.get(solvent, 0)

    if verbose:
        d_contrib['subs'] = len(subs)
        d_contrib['exo'] = count_exo_bonds(mol)
        d_contrib["solvent"] = solvent_values.get(solvent, 0)

    return round(total), d_contrib


def __fieser_calc(
        base,
        factor,
        subs,
        mol,
        solvent,
        verbose,
        mol_type='fieser'
        ):

    total = 0
    d_contrib = {}

    # base values
    total += check_values(base, mol_type, base_value_library)
    total += check_values(factor, mol_type, factor_value_library)
    total += check_values(subs, mol_type, sub_value_library)
    total += count_exo_bonds(mol) * 5

    # solvent
    total += solvent_values.get(solvent, 0)
    
    if verbose:
        # base values
        d_contrib["base"] = check_values(base, mol_type, base_value_library)
        d_contrib["factor"] = check_values(factor, mol_type, factor_value_library)

        # substituents
        d_contrib["subs"] = check_values(subs, mol_type, sub_value_library, single=False)
        d_contrib["exo"] = count_exo_bonds(mol) * 5

    return round(total), d_contrib


def __fieser_kuhn_calc(
        base,
        factor,
        subs,
        mol,
        solvent,
        verbose,
        mol_type='fieser_kuhn'
        ):

    d_contrib = {}

    factors = check_values(factor, mol_type, factor_value_library, False)

    # get ingredients
    m = factors.get("alkyl", 0)
    n = count_conjugated_double_bonds(mol)
    endo = factors.get("endocyclic_double_bond", 0)
    exo = count_exo_bonds(mol)

    total = 114 + 5*m + n*(48 - 1.7*n) - 16.5*endo - 10*exo
    epsilon = 17400 * n

    # solvent
    total += solvent_values.get(solvent, 0)

    if verbose:
        d_contrib["M"] = m
        d_contrib["N"] = n
        d_contrib["endo"] = endo
        d_contrib["exo"] = exo
        d_contrib["epsilon"] = epsilon

    return round(total), d_contrib


def __coumarin_calc(
        base,
        factor,
        subs,
        mol,
        solvent,
        verbose,
        mol_type='woodward'
        ):

    total = 0
    d_contrib = {}

    # base values
    total += check_values(base, mol_type, base_value_library)

    # increments for substituents
    total += max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0)
    total += max([sub["value"] for sub in subs if sub["sub_type"] == "beta"], default=0)
    
    # solvent
    total += solvent_values.get(solvent, 0)

    if verbose:
        # base values
        d_contrib['base'] = check_values(base, mol_type, base_value_library)

        # increments for substituents
        d_contrib["alpha"], d_contrib['alpha_all']  = __get_max_sub(pos="alpha", subs=subs)
        d_contrib["beta"], d_contrib["beta_all"] = __get_max_sub(pos="beta", subs=subs)
        
        # solvent
        d_contrib["solvent"] = solvent_values.get(solvent, 0)

    return round(total), d_contrib


def __get_max_sub(subs, pos='alpha'):

    pos_subs = [sub for sub in subs if sub["sub_type"] == pos]
    pos_max = max(pos_subs, key=lambda x: x["value"], default=None)

    if pos_max:
        max_val, max_sub = pos_max['value'], pos_max['pattern']
    else:
        max_val, max_sub = 0, 'H'

    return max_sub, pos_subs


def draw_images(mol):
    highlight_atom_colors = {}
    highlight_atoms = []

    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.HasProp("used"):
            val = atom.GetProp("used")
            highlight_atoms.append(idx)
            if val == "0":
                highlight_atom_colors[idx] = (0.0, 1.0, 1.0)
            elif val == "1":
                highlight_atom_colors[idx] = (0.0, 1.0, 0.0)
            elif val == "2":
                highlight_atom_colors[idx] = (0.0, 0.0, 1.0)
            elif val == "3":
                highlight_atom_colors[idx] = (1.0, 1.0, 0.0)
            elif val == "4":
                highlight_atom_colors[idx] = (1.0, 0.0, 1.0)

    Chem.rdDepictor.Compute2DCoords(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(500, 500)
    drawer.DrawMolecule(mol, highlightAtoms = highlight_atoms, highlightAtomColors = highlight_atom_colors)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()

    return Image.open(io.BytesIO(png))


def __select_mol_type(mol_type_auto, chromlib=None):

    # Define categories as sets
    categories = {
        'woodward': {'woodward', 'woodward_extended', 'woodward_refine'},
        'fieser': {'fieser'},
        'fieser_kuhn': {'fieser_kuhn'},
    }

    # Reverse lookup: map each rule to its category name
    rule_to_category = {
        rule: cat
        for cat, rules in categories.items()
        for rule in rules
    }

    if not chromlib:
        return mol_type_auto

    detected_cat = rule_to_category.get(mol_type_auto)
    requested_cat = rule_to_category.get(chromlib)

    if detected_cat == requested_cat:
        return chromlib
    else:
        warnings.warn(
            f"Requested rule set '{chromlib}' is incompatible with detected rule set '{mol_type_auto}'. "
            f"Proceeding with detected type '{mol_type_auto}'."
        )
        return mol_type_auto


def predict(
        smiles, 
        solvent, 
        draw=False, 
        verbose=False, 
        chromlib=None
        ):

    """
    Parameters
    ----------
    smiles ... isomeric SMILES string of a molecule
    solvent ... SMILES string of a solvent
    draw ... Boolean either return image or not
    verbose ... Boolean, if True returns dictionary of individual contributions
              i.e. structural descriptors and their increments

    Returns
    -------
    Prediction of absorption maximum in nm and image object with
    highlighted base chromophore and other structural fragments 
    contributing as increments to the final computation

    Examples
    ---------

    smi = 'CC(=O)C=CC1=C(C)CCCC1(C)C'
    solv = 'CCCCCCC'

    nm, img = predict(smi, solv)

    """

    mol_type_auto, mol = classify_type(smiles, general_rules)
    
    # handle rule set selection with user constraints
    mol_type = __select_mol_type(mol_type_auto, chromlib)

    # get baselib and factor lib
    base = get_libData(mol, mol_type, base_library, 0)
    factor = get_libData(mol, mol_type, factor_library, 1)

    if mol_type == "woodward_extended" and factor is not None:
        mol_type = "woodward"

    # get substituent libraries
    if mol_type == 'fieser':
        subs = get_libData(mol, mol_type, sub_library, 2)
        pred, contrib = __fieser_calc(base, factor, subs, mol, solvent, verbose)

    if mol_type == 'fieser_kuhn':
        subs = get_libData(mol, mol_type, sub_library, 2)
        pred, contrib = __fieser_kuhn_calc(base, factor, subs, mol, solvent, verbose)

    elif mol_type == 'woodward':
        subs = get_woodward_sub_values(mol)
        pred, contrib = __woodward_calc(base, factor, subs, mol, solvent, verbose)

    elif mol_type == 'woodward_refine':
        subs = get_woodward_sub_values(mol, sub_values_lib=woodward_refine_sub_values)
        pred, contrib = __woodward_calc(base, factor, subs, mol, solvent, verbose, mol_type='woodward_refine')
        
    elif mol_type == 'woodward_extended':
        subs = get_woodward_sub_values(mol)
        pred, contrib = __woodward_extended_calc(base, factor, subs, mol, solvent, verbose)

    elif mol_type == 'woodward_coumarin':
        subs = get_woodward_sub_values(mol, sub_values_lib=woodward_refine_sub_values)
        pred, contrib = __coumarin_calc(base, factor, subs, mol, solvent, verbose, mol_type='woodward_refine')

    im = draw_images(mol)
    if draw:
        im.show()

    return pred, contrib, im


