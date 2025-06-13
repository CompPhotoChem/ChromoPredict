from typing import Callable, Dict, Optional, Sequence, Union
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from chromopredict import *

def check_pattern(
        mol, 
        patterns: Dict, 
        mark, 
        version: int, 
        mol_type: str,
        ):

    """
    Identify substructure patterns and get increments for patterns

    Parameters
    ----------

    mol: rdkit Mol object
    patterns: dictionary of SMARTS and labels to be matched
    mark:
    version:
    mol_type: 

    Returns
    -------

    Examples
    --------

    Find ...

    """
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
                    break
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


