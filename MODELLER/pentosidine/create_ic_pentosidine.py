from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms
import numpy as np

def is_in_ring(atom):
    return atom.IsInRing()

def get_ring_size(atom):
    mol = atom.GetOwningMol()
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    atom_idx = atom.GetIdx()
    
    min_size = float('inf')
    for ring in atom_rings:
        if atom_idx in ring and len(ring) < min_size:
            min_size = len(ring)
    
    return min_size if min_size != float('inf') else 0

def is_ring_junction(atom):
    ring_info = atom.GetOwningMol().GetRingInfo()
    atom_rings = ring_info.AtomRings()
    atom_idx = atom.GetIdx()
    ring_count = sum(1 for ring in atom_rings if atom_idx in ring)
    return ring_count > 1

def assign_residue_atom_names(mol):
    """
    Enhanced residue atom naming with special handling for rings and branches
    """
    names = {}
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    ring_atoms = set()
    for ring in rings:
        ring_atoms.update(ring)
    
    chain_atoms = []
    ring_systems = []
    other_atoms = []
    
    # First pass: categorize atoms
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if idx in ring_atoms:
            # Create new ring system or add to existing one
            assigned = False
            for ring_sys in ring_systems:
                if any(idx in ring for ring in rings if any(x in ring_sys for x in ring)):
                    ring_sys.add(idx)
                    assigned = True
                    break
            if not assigned:
                ring_systems.append({idx})
        elif atom.GetSymbol() == 'C':
            chain_atoms.append(idx)
        else:
            other_atoms.append(idx)
    
    # Name ring system atoms
    ring_count = 1
    for ring_sys in ring_systems:
        for i, idx in enumerate(sorted(ring_sys)):
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            if symbol == 'C':
                names[idx] = f'CR{ring_count}{chr(65+i)}'  # CR1A, CR1B, etc.
            elif symbol == 'N':
                names[idx] = f'NR{ring_count}{chr(65+i)}'  # NR1A, NR1B, etc.
            elif symbol == 'O':
                names[idx] = f'OR{ring_count}{chr(65+i)}'  # OR1A, OR1B, etc.
        ring_count += 1
    
    # Name chain atoms
    for i, idx in enumerate(chain_atoms):
        if idx not in names:  
            if i == 0:
                names[idx] = 'CA'
            elif i == 1:
                names[idx] = 'CB'
            elif i == 2:
                names[idx] = 'CG'
            elif i == 3:
                names[idx] = 'CD'
            elif i == 4:
                names[idx] = 'CE'
            elif i == 5:
                names[idx] = 'CZ'
            else:
                names[idx] = f'C{i+1}'
    
    # Name other atoms
    hetero_counts = {'O': 1, 'N': 1, 'H': 1}
    for idx in other_atoms:
        if idx not in names:  
            atom = mol.GetAtomWithIdx(idx)
            symbol = atom.GetSymbol()
            count = hetero_counts.get(symbol, 1)
            names[idx] = f'{symbol}{count}'
            hetero_counts[symbol] = count + 1
    
    # Name hydrogens based on their parent atoms
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H' and atom.GetIdx() not in names:
            parent = atom.GetNeighbors()[0]
            parent_name = names[parent.GetIdx()]
            h_count = sum(1 for n in parent.GetNeighbors() if n.GetSymbol() == 'H')
            if h_count == 1:
                names[atom.GetIdx()] = f'H{parent_name}'
            else:
                h_num = sum(1 for n in parent.GetNeighbors() 
                          if n.GetSymbol() == 'H' and n.GetIdx() < atom.GetIdx()) + 1
                names[atom.GetIdx()] = f'H{h_num}{parent_name}'
    
    return names

def get_charmm_atom_type(atom, mol):
    """
    Enhanced CHARMM atom typing with specific types for different chemical environments
    """
    atomic_num = atom.GetSymbol()
    idx = atom.GetIdx()
    
    is_aromatic = atom.GetIsAromatic()
    hybridization = str(atom.GetHybridization())
    num_neighbors = len(atom.GetNeighbors())
    in_ring = is_in_ring(atom)
    ring_size = get_ring_size(atom) if in_ring else 0
    is_junction = is_ring_junction(atom)
    
    # Carbon typing
    if atomic_num == 'C':
        if is_aromatic:
            if is_junction:
                return 'CJ'    # Junction aromatic carbon
            return 'CA'        # Regular aromatic carbon
        elif in_ring:
            if is_junction:
                return 'CT2'   # Ring junction carbon
            elif hybridization == 'SP3':
                return 'CT1'   # Ring SP3 carbon
            else:
                return 'CR'    # Generic ring carbon
        elif hybridization == 'SP3':
            oxygen_neighbors = sum(1 for n in atom.GetNeighbors() if n.GetSymbol() == 'O')
            if oxygen_neighbors > 1:
                return 'CT1H'  # Carbon with multiple OH groups (sugar-like)
            elif oxygen_neighbors == 1:
                return 'CT1'   # Carbon adjacent to oxygen
            else:
                return 'CT3'   # SP3 carbon (methyl-like)
        elif hybridization == 'SP2':
            if any(n.GetSymbol() == 'O' for n in atom.GetNeighbors()):
                return 'C'     # Carbonyl carbon
            elif any(n.GetSymbol() == 'N' for n in atom.GetNeighbors()):
                return 'CN'    # Carbon adjacent to nitrogen
            else:
                return 'CE1'   # SP2 carbon
    
    # Oxygen typing
    elif atomic_num == 'O':
        if num_neighbors == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetHybridization() == Chem.HybridizationType.SP2:
                return 'O'     # Carbonyl oxygen
            else:
                if in_ring:
                    return 'OH2'  # Ring hydroxyl oxygen
                return 'OH1'   # Regular hydroxyl oxygen
        elif num_neighbors == 2:
            if in_ring:
                return 'OR'    # Ring ether oxygen
            return 'OS'        # Regular ether oxygen
    
    # Nitrogen typing
    elif atomic_num == 'N':
        if is_aromatic:
            if is_junction:
                return 'NJ'    # Junction aromatic nitrogen
            return 'NR'        # Regular aromatic nitrogen
        elif hybridization == 'SP3':
            if num_neighbors == 4:
                return 'NH3'   # Ammonium nitrogen
            elif in_ring:
                return 'NT'    # Ring tertiary amine
            else:
                return 'NH1'   # Regular SP3 amine
        elif hybridization == 'SP2':
            if in_ring:
                if is_junction:
                    return 'NN'  # Ring junction SP2 nitrogen
                return 'NR2'     # Ring SP2 nitrogen
            return 'NH2'         # Regular SP2 nitrogen
    
    # Hydrogen typing
    elif atomic_num == 'H':
        parent = mol.GetAtomWithIdx(atom.GetNeighbors()[0].GetIdx())
        parent_symbol = parent.GetSymbol()
        if parent_symbol == 'C':
            if parent.GetIsAromatic():
                return 'HP'    # Aromatic H
            elif is_in_ring(parent):
                return 'HR'    # Ring H
            elif any(n.GetSymbol() == 'O' for n in parent.GetNeighbors()):
                return 'HT'    # H on carbon next to OH
            else:
                return 'HA'    # Aliphatic H
        elif parent_symbol == 'O':
            return 'HO'        # Hydroxyl H
        elif parent_symbol == 'N':
            if parent.GetIsAromatic():
                return 'HNR'   # Aromatic amine H
            elif is_in_ring(parent):
                return 'HNT'   # Ring amine H
            return 'HN'        # Regular amine H
    
    # Default/generic types for other atoms
    return atomic_num

def smiles_to_3d(smiles_string):
    """Convert SMILES to 3D structure using RDKit"""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    mol = Chem.AddHs(mol)
    
    # Generate 3D conformation
    success = AllChem.EmbedMolecule(mol, randomSeed=42)
    if success == -1:
        raise ValueError("Could not generate 3D coordinates")
    
    AllChem.MMFFOptimizeMolecule(mol)
    
    return mol

def get_geometric_parameters(mol):
    """Extract bond distances, angles, and dihedrals from molecule, excluding hydrogens"""
    conf = mol.GetConformer()
    positions = conf.GetPositions()
    
    # Get bonds (excluding H)
    bonds = []
    for bond in mol.GetBonds():
        a1 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        a2 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if a1.GetSymbol() != 'H' and a2.GetSymbol() != 'H':
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            dist = np.linalg.norm(positions[idx1] - positions[idx2])
            bonds.append((idx1, idx2, dist))
    
    # Get angles (excluding H)
    angles = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        idx2 = atom.GetIdx()
        neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if len(neighbors) >= 2:
            for i in range(len(neighbors)):
                for j in range(i+1, len(neighbors)):
                    idx1 = neighbors[i].GetIdx()
                    idx3 = neighbors[j].GetIdx()
                    v1 = positions[idx1] - positions[idx2]
                    v2 = positions[idx3] - positions[idx2]
                    angle = np.degrees(np.arccos(np.dot(v1, v2) / 
                                    (np.linalg.norm(v1) * np.linalg.norm(v2))))
                    angles.append((idx1, idx2, idx3, angle))
    
    # Get dihedrals (excluding H)
    dihedrals = []
    for bond in mol.GetBonds():
        a2 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        a3 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if a2.GetSymbol() == 'H' or a3.GetSymbol() == 'H':
            continue
            
        idx2 = bond.GetBeginAtomIdx()
        idx3 = bond.GetEndAtomIdx()
        
        for n1 in a2.GetNeighbors():
            if n1.GetSymbol() == 'H' or n1.GetIdx() == idx3:
                continue
            idx1 = n1.GetIdx()
            
            for n4 in a3.GetNeighbors():
                if n4.GetSymbol() == 'H' or n4.GetIdx() == idx2 or n4.GetIdx() == idx1:
                    continue
                idx4 = n4.GetIdx()
                
                p1, p2, p3, p4 = [positions[i] for i in (idx1, idx2, idx3, idx4)]
                v1 = p2 - p1
                v2 = p3 - p2
                v3 = p4 - p3
                n1 = np.cross(v1, v2)
                n2 = np.cross(v2, v3)
                angle = np.degrees(np.arctan2(
                    np.dot(np.cross(n1, n2), v2/np.linalg.norm(v2)),
                    np.dot(n1, n2)))
                dihedrals.append((idx1, idx2, idx3, idx4, angle))
    
    return bonds, angles, dihedrals

def get_proper_dihedrals(mol):
    """Identify proper dihedral angles in the molecule, excluding hydrogens"""
    proper_dihedrals = []
    for bond in mol.GetBonds():
        a2 = mol.GetAtomWithIdx(bond.GetBeginAtomIdx())
        a3 = mol.GetAtomWithIdx(bond.GetEndAtomIdx())
        if a2.GetSymbol() == 'H' or a3.GetSymbol() == 'H':
            continue
            
        atom2_idx = bond.GetBeginAtomIdx()
        atom3_idx = bond.GetEndAtomIdx()
        
        atom2_neighbors = [n.GetIdx() for n in a2.GetNeighbors() 
                         if n.GetSymbol() != 'H' and n.GetIdx() != atom3_idx]
        atom3_neighbors = [n.GetIdx() for n in a3.GetNeighbors() 
                         if n.GetSymbol() != 'H' and n.GetIdx() != atom2_idx]
        
        for atom1_idx in atom2_neighbors:
            for atom4_idx in atom3_neighbors:
                proper_dihedrals.append((atom1_idx, atom2_idx, atom3_idx, atom4_idx))
    
    return proper_dihedrals

def get_improper_dihedrals(mol):
    """Identify improper dihedral angles in the molecule, excluding hydrogens"""
    improper_dihedrals = []
    
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
        
        neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        if len(neighbors) >= 3:
            central_idx = atom.GetIdx()
            neighbor_indices = [n.GetIdx() for n in neighbors]
            
            for i in range(len(neighbor_indices)):
                for j in range(i + 1, len(neighbor_indices)):
                    for k in range(j + 1, len(neighbor_indices)):
                        improper_dihedrals.append(
                            (neighbor_indices[i], neighbor_indices[j], central_idx, neighbor_indices[k]))
    
    return improper_dihedrals

def get_all_angles(mol):
    """Get all bond angles in the molecule, excluding hydrogens"""
    angles = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            continue
            
        central_idx = atom.GetIdx()
        neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() != 'H']
        
        if len(neighbors) >= 2:
            for i in range(len(neighbors)):
                for j in range(i + 1, len(neighbors)):
                    angles.append((neighbors[i].GetIdx(), central_idx, neighbors[j].GetIdx()))
    
    return angles

def process_smiles(smiles_string, output_file="structure.str"):
    """Main function to process SMILES and generate CHARMM parameters"""
    try:
        # Generate 3D structure
        mol = smiles_to_3d(smiles_string)
        
        bonds, angles, dihedrals = get_geometric_parameters(mol)
        
        write_charmm_parameters(mol, bonds, angles, dihedrals, output_file)
        
        print(f"Successfully generated CHARMM structure file: {output_file}")
        
        return {
            "n_atoms": mol.GetNumAtoms(),
            "n_bonds": len(bonds),
            "n_angles": len(angles),
            "n_dihedrals": len(dihedrals)
        }
        
    except Exception as e:
        print(f"Error processing SMILES: {str(e)}")
        return None

def calculate_complete_ic_table(mol, atom_names):
    """
    Calculate complete internal coordinates including all angles
    Returns formatted IC table entries
    """
    conf = mol.GetConformer()
    proper_dihedrals = get_proper_dihedrals(mol)
    improper_dihedrals = get_improper_dihedrals(mol)
    all_angles = get_all_angles(mol)
    ic_entries = []
    
    def get_bond_length(idx1, idx2):
        pos1 = conf.GetAtomPosition(idx1)
        pos2 = conf.GetAtomPosition(idx2)
        return np.linalg.norm(pos1 - pos2)
    
    def get_angle(idx1, idx2, idx3):
        pos1 = conf.GetAtomPosition(idx1)
        pos2 = conf.GetAtomPosition(idx2)
        pos3 = conf.GetAtomPosition(idx3)
        v1 = pos1 - pos2
        v2 = pos3 - pos2
        angle = np.degrees(np.arccos(np.dot(v1, v2) / 
                        (np.linalg.norm(v1) * np.linalg.norm(v2))))
        return angle
    
    def get_dihedral(idx1, idx2, idx3, idx4):
        pos1 = conf.GetAtomPosition(idx1)
        pos2 = conf.GetAtomPosition(idx2)
        pos3 = conf.GetAtomPosition(idx3)
        pos4 = conf.GetAtomPosition(idx4)
        
        v1 = pos2 - pos1
        v2 = pos3 - pos2
        v3 = pos4 - pos3
        
        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)
        
        angle = np.degrees(np.arctan2(
            np.dot(np.cross(n1, n2), v2/np.linalg.norm(v2)),
            np.dot(n1, n2)))
        return angle

    for a1, a2, a3 in all_angles:
        r12 = get_bond_length(a1, a2)
        angle = get_angle(a1, a2, a3)
        r23 = get_bond_length(a2, a3)
        
        entry = (f"IC {atom_names[a1]:<4} {atom_names[a2]:<4} {atom_names[a3]:<4} "
                f"{atom_names[a3]:<4} {r12:6.3f} {angle:7.2f} {0.00:8.2f} "
                f"{0.00:7.2f} {0.00:6.3f}")
        ic_entries.append(entry)
    
    for a1, a2, a3, a4 in proper_dihedrals:
        r12 = get_bond_length(a1, a2)
        angle123 = get_angle(a1, a2, a3)
        dihedral = get_dihedral(a1, a2, a3, a4)
        angle234 = get_angle(a2, a3, a4)
        r34 = get_bond_length(a3, a4)
        
        entry = (f"IC {atom_names[a1]:<4} {atom_names[a2]:<4} {atom_names[a3]:<4} "
                f"{atom_names[a4]:<4} {r12:6.3f} {angle123:7.2f} {dihedral:8.2f} "
                f"{angle234:7.2f} {r34:6.3f}")
        ic_entries.append(entry)
    
    for a1, a2, c, a4 in improper_dihedrals:
        r12 = get_bond_length(a1, a2)
        angle123 = get_angle(a1, a2, c)
        dihedral = get_dihedral(a1, a2, c, a4)
        angle234 = get_angle(a2, c, a4)
        r34 = get_bond_length(c, a4)
        
        entry = (f"IC {atom_names[a1]:<4} {atom_names[a2]:<4} *{atom_names[c]:<4} "
                f"{atom_names[a4]:<4} {r12:6.3f} {angle123:7.2f} {dihedral:8.2f} "
                f"{angle234:7.2f} {r34:6.3f}")
        ic_entries.append(entry)
    
    return ic_entries

def write_charmm_parameters(mol, bonds, angles, dihedrals, filename="structure.str"):
    """Write geometric parameters in CHARMM stream file format, excluding hydrogens"""
    # Get CHARMM atom types and residue names (excluding H)
    charmm_types = {atom.GetIdx(): get_charmm_atom_type(atom, mol) 
                   for atom in mol.GetAtoms() if atom.GetSymbol() != 'H'}
    residue_names = assign_residue_atom_names(mol)
    
    with open(filename, 'w') as f:
        f.write("* Structure stream file generated from SMILES\n*\n\n")
        
        f.write("MASS -1 BEGIN 1.0\n\n")
        
        f.write("read sequence card\n* Sequence\n*\n")
        f.write("1\n")  # number of residues
        f.write("MOL  ")  # residue name
        f.write("\n\n")
        
        f.write("generate MOL setup\n\n")
        
        f.write("! Coordinates\n")
        conf = mol.GetConformer()
        for atom in mol.GetAtoms():
            if atom.GetSymbol() != 'H':
                i = atom.GetIdx()
                pos = conf.GetAtomPosition(i)
                atom_type = charmm_types[i]
                atom_name = residue_names[i]
                f.write(f"coord atom {atom_name} MOL {pos.x:.3f} {pos.y:.3f} {pos.z:.3f}  ! {atom_type}\n")
        f.write("\n")
        
        f.write("! Bonds\n")
        for idx1, idx2, dist in bonds:
            name1 = residue_names[idx1]
            name2 = residue_names[idx2]
            f.write(f"BOND {name1} {name2} ! {dist:.3f}\n")
        f.write("\n")
        
        f.write("! Angles\n")
        for idx1, idx2, idx3, angle in angles:
            name1 = residue_names[idx1]
            name2 = residue_names[idx2]
            name3 = residue_names[idx3]
            f.write(f"angle MOL {name1} MOL {name2} MOL {name3} ! {angle:.2f}\n")
        f.write("\n")
        
        f.write("! Dihedrals\n")
        for idx1, idx2, idx3, idx4, angle in dihedrals:
            name1 = residue_names[idx1]
            name2 = residue_names[idx2]
            name3 = residue_names[idx3]
            name4 = residue_names[idx4]
            f.write(f"IMPR {name1} {name2} {name3} {name4} ! {angle:.2f}\n")
        f.write("\n")
        
        f.write("! Internal Coordinates Table\n")
        ic_entries = calculate_complete_ic_table(mol, residue_names)
        for entry in ic_entries:
            if not any(name.startswith('H') for name in entry.split()[1:5]):  
                f.write(entry + "\n")
        f.write("\n")

if __name__ == "__main__":
    smiles = "CCCCNc2nc1ncccc1[nH]2"
    result = process_smiles(smiles, "pentosidine_ic.str")
    if result:
        print("\nStructure statistics:")
        for key, value in result.items():
            print(f"{key}: {value}")
