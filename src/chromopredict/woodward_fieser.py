from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import io
from PIL import Image
# General informations before use: 
# the elements of each dictionary are to be sorted from the highest down, to avoid false calculation

# SMARTS

general_rules = {
    # "cyanide": "[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1-,:[#7](-[#6]-[#6])-,:[#6]",
    "acridine": "c1cc2[n,s]c3ccccc3[c,n,s]c2cc1",
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
    "homoanular_cyclic": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1", #[#6;r]=[#6;r]-[#6;r]=[#6;r]
    "acyclic_diene": "[#6;!R]=[#6;!R]-[#6;!R]=[#6;!R]",
    "heteroanular_cyclic": "[#6]=[#6;R]-[#6;R]=[#6]"
}

fieser_kuhn_base = {
    "polyene_4_conj": "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"
}

cyanide_base = {
    "temp": "[#6]=[#6]"
}

factors = {
    "homoanular_cyclodiene": "[#6]1=[#6]-[#6]=[#6]-[#6]-[#6]1",
    "conjugated_double_bond": "[#6]=[#6]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

fieser_kuhn_factors = {
    "alkyl" : "[#6;X4]",
    "endocyclic_double_bond": "[#6;R]=[#6;R]"
}

cyanide_factors = {
    "temp": "[#6]=[#6]"
}

fieser_subs = {
    "dialkylamine": "[#7](-[#6])(-[#6])",
    "thioether": "[#16]-[#6]",
    "phenoxy": "[#8]-c1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
    "halogen": "[F,Cl,Br,I]",
    "alkoxy": "[#8]-[#6]",
    "alkyl": "[#6;X4]",
    "carboxy": "[#6](=[#8])-[#8H]",
}

woodward_subs = {
    "alkoxy": "[#8]-[#6]",
    "alkyl": "[#6;X4]",
    "bromo": "[Br]",
    "carboxy": "[#6](=[#8])-[#8H]",
    "chloro": "[Cl]",
    "hydroxy": "[#8H]"
}

cyanide_subs = {
    "temp": "[#6]=[#6]"
}

# Werte:

woodward_base_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "alpha_beta_aldehyd": 210,
    "alpha_beta_ketone": 215,
    "alpha_beta_ester": 195,
    "cyclopentenone": 202,
    "cyclohexanone": 215, # doi: 10.1021/jo01164a003
    "alpha_beta_acid" : 195, # nur eine halbwegs verlässliche quelle mit 195 (https://www.researchgate.net/publication/372535169_The_Woodward_Fisher_Regulation_for_Calculating_Absorption_Maxima)
}

fieser_base_values = { # doi: 10.1021/jo01164a003
    "acyclic_diene": 217,
    "homoanular_cyclic": 253,
    "heteroanular_cyclic": 214
}

fieser_kuhn_base_values = { #https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "polyene_4_conj": 114
}

cyanide_base_values = {
    "temp": 0
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

cyanide_factor_values = {
    "temp": 0
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

cyanide_sub_values = {
    "temp": 0
}

solvent_values = { #doi: 10.1021/ja01849a066, https://de.wikipedia.org/wiki/Woodward-Fieser-Regeln
    "O": -8,         
    "CO": -5,
    "CCO": -1,
    "ClC(Cl)Cl": 0,
    "O1CCOCC1": 5,
    "CCOCC": 7,
    "CCCCCC": 11,
    "C1CCCCC1": 11   
}

# Dictionarys

base_library = {
    "fieser": fieser_base,
    "fieser_kuhn": fieser_kuhn_base,
    "woodward_extended": woodward_base,
    "woodward": woodward_base,
    "cyanide": cyanide_base
}

sub_library = {
    "fieser": fieser_subs,
    "fieser_kuhn": "None",
    "woodward_extended": woodward_subs,
    "woodward": woodward_subs,
    "cyanide": cyanide_subs
}

factor_library = {
    "fieser": factors,
    "fieser_kuhn": fieser_kuhn_factors,
    "woodward_extended": factors,
    "woodward": factors,
    "cyanide": cyanide_factors
}

base_value_library = {
    "fieser": fieser_base_values,
    "fieser_kuhn": fieser_kuhn_base_values,
    "woodward_extended": woodward_base_values,
    "woodward": woodward_base_values,
    "cyanide": cyanide_base_values
}

sub_value_library = {
    "fieser": fieser_sub_values,
    "fieser_kuhn": "None",
    "woodward_extended": woodward_sub_values,
    "woodward": woodward_sub_values,
    "cyanide": cyanide_sub_values
}

factor_value_library = {
    "fieser": fieser_factor_values,
    "fieser_kuhn": fieser_kuhn_factor_values,
    "woodward_extended": woodward_factor_values,
    "woodward": woodward_factor_values,
    "cyanide": cyanide_factor_values
}

# Extras

woodward_markers = {
    1: "alpha",
    2: "beta",
    3: "gamma",
    4: "higher"
}

hetero_atoms = [6, 7, 8, 9, 15, 16, 17, 35, 53]

# Functions:

def check_pattern(mol, patterns, mark, version, mol_type):
    result = []
    addFactor = 0
    somethingFound = True

    if version > 0:
        if mol_type == "woodward":
            addFactor = 1
        curProps = [str(i) for i in range(version + addFactor)]
        newProp = str(version)

    while somethingFound:
        somethingFound = False
        big_break = False

        for patternName, smarts in patterns.items():
            pattern = Chem.MolFromSmarts(smarts)

            if pattern is None:
                raise ValueError("Smarts to Mol conversion error")
            
            matches = mol.GetSubstructMatches(pattern)

            # print("Matches: ", matches, " Version: ", version, " PatternName: ", patternName)

            for match in matches:
                if version == 0: 
                    result.append(patternName)

                    if mark:
                        for idx in match:
                            if patternName != "homoanular_cyclic":  # "homoanular_cyclodiene"
                                mol.GetAtomWithIdx(idx).SetProp("used", "0")
                            else:
                                atom = mol.GetAtomWithIdx(idx)
                                in_double_bond = False
                                for bond in atom.GetBonds():
                                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        in_double_bond = True
                                        break
                                if in_double_bond:
                                    mol.GetAtomWithIdx(idx).SetProp("used", "0")
                    big_break = True
                    break
                elif version in [1, 2]:
                    if all(not mol.GetAtomWithIdx(idx).HasProp("used") for idx in match):
                        accepted_neighbor = False
                        for idx in match:
                            atom = mol.GetAtomWithIdx(idx)
                            for neighbor in atom.GetNeighbors():
                                if neighbor.HasProp("used") and neighbor.GetProp("used") in curProps:
                                    accepted_neighbor = True
                                    break
                            if accepted_neighbor:
                                break
                        if accepted_neighbor:
                            result.append(patternName)

                            matches = tuple(x for x in matches if x != match)
                            somethingFound = True
                            big_break = True

                            if mark:
                                for idx in match:
                                    mol.GetAtomWithIdx(idx).SetProp("used", newProp)
                    else:
                        matches = tuple(x for x in matches if x != match)

                if big_break:
                    break
            if big_break:
                break

    if not result:
        return None
    
    if len(result) == 1:
        result = result[0]

    return result

def use_library(bookTitle, library):
    book = library.get(bookTitle)

    if book:
        return book
    else:
        print(bookTitle)
        print("\n")
        print(library)
        raise ValueError("No book found")

def classify_type(smiles, rules):
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Smiles to Mol conversion error")
    
    type = check_pattern(mol, rules, False, 0, None)

    return type, mol

def check_values(var, mol_type, library, single = True):
    if var is None:
        if single:
            return 0
        else:
            return{"Empty": 0}
    
    data = use_library(mol_type, library)

    if isinstance(var, list):
        temp_dict = {}
        var2 = list(set(var))
        for temp in var2:
            temp_dict[temp] = data[temp] * var.count(temp)
        if single:
            return sum(temp_dict.values())
        return temp_dict
    else:
        if single:
            return data[var]
        else:
            return {var: data[var]}

def get_libData(mol, mol_type, library, version):
    data = use_library(mol_type, library)

    if data == "None":
        return None
    return check_pattern(mol, data, True, version, mol_type)

def count_exo_bonds(mol):
    count = 0

    for bond in mol.GetBonds():

        if bond.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue

        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()

        if a1.GetAtomicNum() != 6 or a2.GetAtomicNum() != 6:
            continue

        if a1.IsInRing() != a2.IsInRing():
            count += 1
            continue

        if a1.IsInRing() and a2.IsInRing():
            ring_info = mol.GetRingInfo()
            atom_rings = ring_info.AtomRings()
            for ring in atom_rings:
                if ((a1.GetIdx() in ring) != (a2.GetIdx() in ring)):
                    count += 1
    return count

def count_single_double(mol):
    smarts = "[#6]=[#6]-[#6]=[#6]-[#6]=[#6]-[#6]=[#6]"
    count = 10
    exists = True

    while exists:
        pattern = Chem.MolFromSmarts(smarts)

        matches = mol.GetSubstructMatches(pattern)

        if matches:
            smarts += "-[#6]=[#6]"
            count += 2
        else:
            exists = False
    return count

def identify_atom(mol, pos, pos_type, double_indexes):
    atom = mol.GetAtomWithIdx(pos)

    if pos in double_indexes:
        if atom.HasProp("used") and atom.GetProp("used") in ["0", "1"]:
            if any((atom.GetNeighbors()[i].GetAtomicNum() in hetero_atoms) and not (atom.GetNeighbors()[i].HasProp("sub_type")) for i in range(len(atom.GetNeighbors()))):
                atom.SetProp("sub_type", woodward_markers.get(pos_type))
                
                if pos_type < 4:
                    pos_type += 1

                for neighbor in atom.GetNeighbors():
                    if neighbor.HasProp("sub_type"):
                        continue
                    if neighbor.GetAtomicNum() == 6:
                        identify_atom(mol, neighbor.GetIdx(), pos_type, double_indexes)

def get_double_idx(mol):
    idx = set()

    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            idx.add(bond.GetBeginAtomIdx())
            idx.add(bond.GetEndAtomIdx())
    
    return idx

def mark_atoms(mol):
    start = Chem.MolFromSmarts("[#6]=[#8]")

    matches = mol.GetSubstructMatches(start)

    doubIdx = get_double_idx(mol)

    for match in matches:
        if all(mol.GetAtomWithIdx(idx).HasProp("used") and (mol.GetAtomWithIdx(idx).GetProp("used") == "0") for idx in match):
            break
    
    if mol.GetAtomWithIdx(match[0]).GetAtomicNum() == 6:
        startIdx = match[0]
    else:
        startIdx = match[1]
    
    mol.GetAtomWithIdx(startIdx).SetProp("sub_type", "0")

    for neighbor in mol.GetAtomWithIdx(startIdx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            identify_atom(mol, neighbor.GetIdx(), 1, doubIdx)
                
def get_woodward_sub_values(mol):
    result = []

    for patternName, smarts in woodward_subs.items():
        pattern = Chem.MolFromSmarts(smarts)

        if pattern is None:
            raise ValueError("Smarts to Mol conversion error")
        
        matches = mol.GetSubstructMatches(pattern)

        for match in matches:
            useful = False

            if any(mol.GetAtomWithIdx(idx).HasProp("used") for idx in match):
                continue
            for idx in match:
                for neighbor in mol.GetAtomWithIdx(idx).GetNeighbors():
                    if neighbor.HasProp("sub_type"):
                        if neighbor.GetProp("sub_type") == "0":
                            continue
                        # print("Name: ", patternName, ", Position: ", neighbor.GetProp("sub_type"), ", Wert: ", woodward_sub_values[neighbor.GetProp("sub_type")][patternName])
                        result.append({"value": woodward_sub_values[neighbor.GetProp("sub_type")][patternName], "sub_type": neighbor.GetProp("sub_type")})
                        useful = True
                        break
            if useful:
                for idx in match:
                    mol.GetAtomWithIdx(idx).SetProp("used", "2")
    return result      


def combine_data(base, factor, subs, mol_type, mol, solvent, debug, extended):
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
    
    elif mol_type == "woodward" or  (not extended and mol_type == "woodward_extended"):
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

    
    elif mol_type == "woodward_extended":
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

    drawer = rdMolDraw2D.MolDraw2DCairo(1500, 1500)
    drawer.DrawMolecule(mol, highlightAtoms = highlight_atoms, highlightAtomColors = highlight_atom_colors)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()

    return Image.open(io.BytesIO(png))

def main(smiles, solvent, draw=False, debug = False, extended = True):
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

    if draw:
        im = draw_images(mol)
        im.show()

    return combine_data(base, factor, subs, mol_type, mol, solvent, debug, extended)  #, im

smiles = "CC=CC=CC=CC(=O)O"
solvent = "CCO"

print(main(smiles, solvent, debug = True, extended = False))

# neue test smiles:
#I CC(=O)C1=C(C)CCCC1	CCO	247	
#II	CC(C)=C1CCC(C)CC1=O	CCO	252	
#III	CC1=CC(=O)CC(C)(C)C1	CCO	235	
#IV	CC1=C2CCCCC2CCC1=O	CCCCCCC	239	
#V	CC=CC=CC=O	CCCCCC	261	
#VI	CC(=O)C=CC1=C(C)CCCC1(C)C	CCCCCCC	283	
#VII	CC=CC=CC=CC(=O)O	CCO	297	
#VIII	CC(C)=CC(=O)C=C(C)C	CCCCCC	260	
#IX	O=CC(Cl)=C(Cl)C(=O)O	O	263	
#X	COC(=O)C1=C(C)C(C)=C(C(=O)OC)CC1	CCCCCC	300	
#XI	O=C(O)C=CC=CC(=O)O	CO	259	
#XII	CC(C=CC=O)=CC(=O)O	CO	271	
#XIII	CC(C=O)=CC=CC(C)=CC(=O)O
   
   
   
   
   
   
   
   
   
   
    
#print(main(smiles, solvent))

# smiles = [
#     "CC=CC=CC", 227
#     "C=C1CCCC=C1", 229
#     "CC13CCCCC1=CC=C2CCCCC23", 283 
#     "CC13CCC=CC1=CCC2CCCCC23", 234
#     "CC(=O)OC1=CCC3(C)C(=C1)C=CC2CCCCC23", #309
#     "C=C1CCCCC1=C", 234
#     "[N+](C)=C1C=CC=CC1=CC=CC=C2C=CC=CC2=[N+](C)C", 407
#
#     "CC(C=CC=O)Cl", woodward 228, has 222
#     "OCC(C=CC=O)C" woodward 222, has 222
#     "C=C1C(Cl)=C(C(C(O)=CCl)=C(\O)C=O)C(Cl)CC1Br" woodward 403, has 405
#     "C=C(C)C(Cl)=C(Br)C(Br)=C(C=O)OC" woodward 363, has 375
#     "C=CCC=CC=CC=O" 240, has 258
#     ]

# used_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp("used")]
# print(f"Match found for pattern '{patternName}', marked atoms: {used_atoms}")
# img = Draw.MolToImage(mol, highlightAtoms=used_atoms, size=(300, 300))
# img.show()

# used_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if (atom.HasProp("used") and atom.GetProp("used") == "0")]
# used_atom_colors = {idx: (0.0, 1.0, 0.0) for idx in used_atoms}
# img = Draw.MolToImage(mol, highlightAtoms=used_atoms, highlightAtomColors=used_atom_colors, size=(300, 300))
# img.show()
