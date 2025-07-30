import rdkit
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image

from chromopredict.strucfeatures import *
from chromopredict import *

def combine_data(
        base, 
        factor, 
        subs, 
        mol_type, 
        mol, 
        solvent, 
        debug, 
        extended
        ):

    """
    Parameters
    ----------

    Returns
    -------

    Prediction of absorption maximum in nm and image object with
    highlighted base chromophore and other structural fragments 
    contributing as increments to the final computation


    Examples
    --------

    ...

    """

    total = 0
    d_contrib = {}
    
    if mol_type == "fieser":
        total += check_values(base, mol_type, base_value_library)
        total += check_values(factor, mol_type, factor_value_library)
        total += check_values(subs, mol_type, sub_value_library)
        total += count_exo_bonds(mol) * 5

        if debug:
            d_contrib["base"] = check_values(base, mol_type, base_value_library)
            d_contrib["factor"] = check_values(factor, mol_type, factor_value_library)
            d_contrib["subs"] = check_values(subs, mol_type, sub_value_library)
            d_contrib["exo"] = count_exo_bonds(mol) * 5
    
    elif mol_type == "fieser_kuhn":
        factors = check_values(factor, mol_type, factor_value_library, False)
        m = factors.get("alkyl", 0)
        n = count_single_double(mol)
        endo = factors.get("endocyclic_double_bond", 0)
        exo = count_exo_bonds(mol)
        total = 114 + 5 * m + 48 * n - 1.7 * n ** 2 - 16.5 * endo - 10 * exo
        epsilon = 17400 * n
        print("e_max: ", epsilon)

        if debug:
            d_contrib["M"] = m
            d_contrib["N"] = n
            d_contrib["endo"] = endo
            d_contrib["exo"] = exo
    
#    elif multi and mol_type in ("woodward", "woodward_extended"):
#        total += check_values(base, mol_type, base_value_library)
#        total += check_values(factor, mol_type, factor_value_library)
#        total += max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0)
#        total += sum([sub["value"] for sub in subs if sub["sub_type"] == "beta"])
#        total += sum([sub["value"] for sub in subs if sub["sub_type"] == "gamma"])
#        total += sum([sub["value"] for sub in subs if sub["sub_type"] == "higher"])
#        total += 5 * count_exo_bonds(mol)
#
#        if debug:
#            d_contrib['base'] = check_values(base, mol_type, base_value_library)
#            d_contrib["factor"] = check_values(factor, mol_type, factor_value_library)
#            d_contrib["alpha"] =  max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0)
#            d_contrib["beta"] = sum([sub["value"] for sub in subs if sub["sub_type"] == "beta"])
#            d_contrib["gamma"] = sum([sub["value"] for sub in subs if sub["sub_type"] == "gamma"])
#            d_contrib["higher"] = sum([sub["value"] for sub in subs if sub["sub_type"] == "higher"])
#            d_contrib["exo"] =  5 * count_exo_bonds(mol)


    elif not extended and mol_type in ("woodward", "woodward_extended"):
        total += check_values(base, mol_type, base_value_library)
        total += check_values(factor, mol_type, factor_value_library)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "beta"], default=0)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "gamma"], default=0)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "higher"], default=0)
        total += 5 * count_exo_bonds(mol)

        if debug:
            d_contrib['base'] = check_values(base, mol_type, base_value_library)
            d_contrib["factor"] = check_values(factor, mol_type, factor_value_library)
            d_contrib["exo"] =  5 * count_exo_bonds(mol)

            def get_max_sub(pos='alpha'):

                pos_subs = [sub for sub in subs if sub["sub_type"] == pos]
                pos_max = max(pos_subs, key=lambda x: x["value"], default=None)

                if pos_max:
                    max_val, max_sub = pos_max['value'], pos_max['pattern']
                else:
                    max_val, max_sub = 0, 'H'

                return max_sub

            d_contrib["alpha"] = get_max_sub("alpha")
            d_contrib["beta"] = get_max_sub("beta")
            d_contrib["gamma"] = get_max_sub("gamma")
            d_contrib["higher"] = get_max_sub("higher")

            #d_contrib["beta"] = max([sub["value"] for sub in subs if sub["sub_type"] == "beta"], default=0)
            #d_contrib["gamma"] = max([sub["value"] for sub in subs if sub["sub_type"] == "gamma"], default=0)
            #d_contrib["higher"] = max([sub["value"] for sub in subs if sub["sub_type"] == "higher"], default=0)
    
    elif extended and mol_type in ("woodward", "woodward_extended"):
        total += 212.82 
        total += 11.37 * len(subs)
        total += 4.79 * count_exo_bonds(mol)

        if debug:
            d_contrib['subs'] = len(subs)
            d_contrib['exo'] = count_exo_bonds(mol)
    
    total += solvent_values.get(solvent, 0)

    if debug:
        d_contrib["solvent"] = solvent_values.get(solvent, 0)

    return round(total), d_contrib

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

def predict(
        smiles, 
        solvent, 
        draw=False, 
        debug = False, 
        extended = True,
        refine = None
        ):

    """
    Parameters
    ----------
    smiles ... isomeric SMILES string of a molecule
    solvent ... SMILES string of a solvent
    draw ... Boolean either return image or not
    debug ... Boolean, if True returns dictionary of individual contributions
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

    mol_type, mol = classify_type(smiles, general_rules)
    
    if mol_type is None:
        return None

    base = get_libData(mol, mol_type, base_library, 0)
    factor = get_libData(mol, mol_type, factor_library, 1)

    if mol_type == "woodward_extended" and factor is not None:
        mol_type = "woodward"

    if mol_type not in ["woodward", "woodward_extended"]:
        subs = get_libData(mol, mol_type, sub_library, 2)
    else:
        subs = get_woodward_sub_values(mol)
        
    im = draw_images(mol)
    
    if draw:
        im.show()

    pred, contrib = combine_data(base, factor, subs, mol_type, mol, solvent, debug, extended)

    return pred, contrib, im


def woodward_refine_predict(
        smiles, 
        solvent, 
        draw=False, 
        debug = False, 
        ):

    """
    Parameters
    ----------

    smiles ... isomeric SMILES string of a molecule
    solvent ... SMILES string of a solvent
    draw ... Boolean either return image or not
    debug ... Boolean, if True returns dictionary of individual contributions
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
    
    if 'woodward' in mol_type_auto:
        base = get_libData(mol, 'woodward_refine', base_library, 0)
        factor = get_libData(mol, 'woodward_refine', factor_library, 1)
        subs = get_woodward_sub_values(mol, sub_values_lib=woodward_refine_sub_values)
        
        im = draw_images(mol)
        if draw:
            im.show()

        mol_type = 'woodward'
        pred, contrib = combine_data(base, factor, subs, mol_type, mol, solvent, debug, extended=False)
    
        return pred, contrib, im

    else:
        return None



