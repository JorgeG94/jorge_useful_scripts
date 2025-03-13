import json
import numpy as np
from numpy import linalg
import sys

COULOMB_CONSTANT_KJMOL = 138.935458


def read_protein_exess_topology(protein_exess_topology):
    with open(protein_exess_topology, "r") as f:
        data = json.load(f)
    
    natoms = max(data["fragments"][-1])+1

    frag_id = [0 for _ in range(0, natoms)]

    frag_count = 0
    for fragment in data["fragments"]:
        for atom in fragment:
            frag_id[atom] = frag_count
        frag_count += 1
    
    new_fragments = data["fragments"].copy()

    # for connection in data["connectivity"]: # broken bonds, append hydrogen caps
    #     atom1, atom2, _ = connection
    #     frag1 = frag_id[atom1]
    #     frag2 = frag_id[atom2]

    #     new_fragments[frag1].append(atom2 + natoms)
    #     new_fragments[frag2].append(atom1 + natoms)

    return new_fragments, natoms

def determine_protein_charges(protein_gmx_topology, natoms):
    with open(protein_gmx_topology, "r") as f:
        contents = f.readlines()

    start = False
    charges = [0 for _ in range(0, natoms)]
    for line in contents:
        # print(line)
        if not start and "[ atoms ]" in line:
            start = True
            continue
        elif start and not line.startswith(";") and len(line.split()) > 0: 
            # print(line)
            line_split = line.split()
            iatom = int(line_split[0]) - 1
            atom_charge = float(line_split[6])
            charges[iatom] = atom_charge
        elif start and len(line.split()) == 0:
            break
    
    return charges

def read_complex_pdb(pdb_file):
    with open(pdb_file, "r") as f:
        contents = f.readlines()

    coordinates = []
    for line in contents:
        if line.startswith("ATOM "):
            line_split = line.split()
            # iatom = int(line_split[1]) - 1

            if " UNL " in line:
                x = float(line_split[-5])
                y = float(line_split[-4])
                z = float(line_split[-3])
            else:
                x = float(line_split[-6])
                y = float(line_split[-5])
                z = float(line_split[-4])
            coordinates.append([x,y,z])
    return coordinates

def read_ligand_exess_topology(ligand_exess_topolgy, protein_natoms):
    with open(ligand_exess_topolgy, "r") as f:
        data = json.load(f)
    
    natoms = max(data["fragments"][-1])+1

    frag_id = [0 for _ in range(0, natoms)]

    frag_count = 0
    for fragment in data["fragments"]:
        for atom in fragment:
            frag_id[atom] = frag_count
        frag_count += 1
    
    new_fragments = data["fragments"][0].copy()
    new_fragments = [x + protein_natoms for x in new_fragments]
    return new_fragments, natoms

def determine_ligand_charges(ligand_gmx_topology, natoms):
    with open(ligand_gmx_topology, "r") as f:
        contents = f.readlines()

    start = False
    charges = [0 for _ in range(0, natoms)]
    for line in contents:
        if not start and "[ atoms ]" in line:
            start = True
            continue
        elif start and not line.startswith(";") and len(line.split()) > 0: 
            line_split = line.split()
            iatom = int(line_split[0]) - 1
            atom_charge = float(line_split[6])
            charges[iatom] = atom_charge
        elif start and len(line.split()) == 0:
            break
    
    return charges

def distance_pair(atom1, atom2, coordinates):
    v1 = np.array(coordinates[atom1])
    v2 = np.array(coordinates[atom2])
    return linalg.norm(v1 - v2)

def coulomb_interaction_energy(protein_fragments, protein_atom_charges, protein_natoms, ligand_atom_charges, ligand_natoms, coordinates):

    ref_frag = len(protein_fragments)
    dimer_energy = {}
    for ifrag, fragment in enumerate(protein_fragments):
        coulomb_energy = 0.0
        for patom in fragment:
            for latom in range(0, ligand_natoms):
                r = distance_pair(patom, latom + protein_natoms, coordinates) * 0.1
                patom_charge = protein_atom_charges[patom]
                latom_charge = ligand_atom_charges[latom]

                contribution = (COULOMB_CONSTANT_KJMOL * patom_charge * latom_charge) / r
                # print(contribution)
                coulomb_energy += contribution
        dimer_energy[(min([ifrag, ref_frag]), max([ifrag, ref_frag]))] = coulomb_energy
    return dimer_energy

def protein_fragment_dipoles(protein_fragments, protein_atom_charges, coordinates):
    protein_dipole_info = {}


    for ifrag, fragment in enumerate(protein_fragments):
        r_vector = np.array([0.0, 0.0, 0.0])
        dipole = np.array([0.0, 0.0, 0.0])
        fragment_charge = 0.0
        for iatom in fragment:
            iatom_charge = protein_atom_charges[iatom]
            coordinate = np.array(coordinates[iatom])

            dipole += (iatom_charge * coordinate)
            r_vector += (coordinate)
            fragment_charge += iatom_charge
        if fragment_charge == 0.0:
            r_vector = np.array([0.0, 0.0, 0.0])
        else:
            r_vector = (1.0 / len(fragment)) * r_vector
        protein_dipole_info[ifrag] = [dipole, r_vector]
    return protein_dipole_info

def ligand_dipoles(ligand_natoms, protein_natoms, ligand_atom_charges, coordinates):

    r_vector = np.array([0.0, 0.0, 0.0])
    dipole = np.array([0.0, 0.0, 0.0])
    ligand_charge = 0.0
    for iatom in range(0, ligand_natoms):
        iatom_charge = ligand_atom_charges[iatom]
        coordinate = np.array(coordinates[iatom + protein_natoms])

        dipole += (iatom_charge * coordinate)
        r_vector += (1 * coordinate)
        ligand_charge += iatom_charge
    r_vector = (1.0 / ligand_natoms) * r_vector

    return [dipole, r_vector]


def calculate_dipole_dipole_interactions(protein_dipole_info, ligand_dipole_info):

    ligand_dipole = ligand_dipole_info[0]
    ligand_rvector = ligand_dipole_info[1]
    dimer_dipole_dipole_interactions = {}
    for ifrag, dipole_info in protein_dipole_info.items():
        frag_dipole = dipole_info[0]
        frag_rvector = dipole_info[1]

        r_vector = frag_rvector - ligand_rvector
        r = np.linalg.norm(r_vector)

        r_hat = r_vector / r

        p1_dot_p2 = np.dot(frag_dipole, ligand_dipole)
        p1_dot_r_hat = np.dot(frag_dipole, r_hat)
        p2_dot_r_hat = np.dot(ligand_dipole, r_hat)

        interaction_energy = COULOMB_CONSTANT_KJMOL * (1 / r**3) * (p1_dot_p2 - 3 * (p1_dot_r_hat * p2_dot_r_hat) / r**2)
        dimer_dipole_dipole_interactions[ifrag] = interaction_energy
    return dimer_dipole_dipole_interactions



if __name__ == "__main__":
    protein_exess_topology = sys.argv[1]
    protein_gmx_topology = sys.argv[2]
    pdb_file = sys.argv[3]
    ligand_exess_topolgy = sys.argv[4]
    ligand_gmx_topology = sys.argv[5]

    protein_fragments, protein_natoms = read_protein_exess_topology(protein_exess_topology)
    # print(f"protein_natoms: {protein_natoms}")
    protein_charges = determine_protein_charges(protein_gmx_topology, protein_natoms)
    coordinates = read_complex_pdb(pdb_file)
    # print(f"len coordinates {len(coordinates)}")
    ligand_fragment, ligand_natoms = read_ligand_exess_topology(ligand_exess_topolgy, protein_natoms)
    # print(f"ligand_natoms {ligand_natoms}")
    ligand_charges = determine_ligand_charges(ligand_gmx_topology, ligand_natoms)
    # print(protein_charges)
    # print(ligand_charges)
    dimer_coulomb = coulomb_interaction_energy(protein_fragments, protein_charges, protein_natoms, ligand_charges, ligand_natoms, coordinates)

    protein_dipole_info = protein_fragment_dipoles(protein_fragments, protein_charges, coordinates)
    ligand_dipole_info = ligand_dipoles(ligand_natoms, protein_natoms, ligand_charges, coordinates)

    dimer_dipole_dipole_interactions = calculate_dipole_dipole_interactions(protein_dipole_info, ligand_dipole_info)


    fragment_energies = []

    for dimer, energy in dimer_coulomb.items():
        coulomb_energy = energy
        dipole_dipole_interaction_energy = dimer_dipole_dipole_interactions[min(list(dimer))]
        dimer_info = {
            "fragments": list(dimer),
            "energy": coulomb_energy + dipole_dipole_interaction_energy
        }
        fragment_energies.append(dimer_info)

    with open("dimer_coulomb_dipole_dipole_energy.json", "w") as f:
        json.dump(fragment_energies, f, indent=4)







