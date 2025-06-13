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
    
    if mol_type == "fieser":
        total += check_values(base, mol_type, base_value_library)
        total += check_values(factor, mol_type, factor_value_library)
        total += check_values(subs, mol_type, sub_value_library)
        total += count_exo_bonds(mol) * 5

        if debug:
            print("Fieser-Calculation:")
            print("Base:", check_values(base, mol_type, base_value_library))
            print("Factor", check_values(factor, mol_type, factor_value_library))
            print("Subs:", check_values(subs, mol_type, sub_value_library))
            print("Exo-Value:", count_exo_bonds(mol) * 5)
    
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
            print("Fieser-Kuhn-Calculation:")
            print("M:", m)
            print("N:", n)
            print("Endo:", endo)
            print("Exo:", exo)
    
    elif not extended and mol_type in ("woodward", "woodward_extended"):
        total += check_values(base, mol_type, base_value_library)
        total += check_values(factor, mol_type, factor_value_library)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "beta"], default=0)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "gamma"], default=0)
        total += max([sub["value"] for sub in subs if sub["sub_type"] == "higher"], default=0)
        total += 5 * count_exo_bonds(mol)

        if debug:
            print("Woodward-Calculation (mol-Type =", mol_type, "):")
            print("Base:", check_values(base, mol_type, base_value_library))
            print("Factor:", check_values(factor, mol_type, factor_value_library))
            print("Alpha:", max([sub["value"] for sub in subs if sub["sub_type"] == "alpha"], default=0))
            print("Beta:", max([sub["value"] for sub in subs if sub["sub_type"] == "beta"], default=0))
            print("Gamma:", max([sub["value"] for sub in subs if sub["sub_type"] == "gamma"], default=0))
            print("Higher:", max([sub["value"] for sub in subs if sub["sub_type"] == "higher"], default=0))
            print("Exo:", 5 * count_exo_bonds(mol))

    
    elif extended and mol_type in ("woodward", "woodward_extended"):
        total += 212.82 
        total += 11.37 * len(subs)
        total += 4.79 * count_exo_bonds(mol)

        if debug:
            print("Woodward-Extended-Calculation:")
            print("Subs:", len(subs))
            print("Exo:", count_exo_bonds(mol))
    
    total += solvent_values.get(solvent, 0)

    if debug:
        print("Solvent:", solvent_values.get(solvent, 0))

    return round(total)

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
        extended = True
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
        mark_atoms(mol)
        subs = get_woodward_sub_values(mol)
        
    im = draw_images(mol)
    
    if draw:
        im.show()

    return combine_data(base, factor, subs, mol_type, mol, solvent, debug, extended), im


