import json
import sys
import os
from math import sqrt
from itertools import combinations

covalent_radii = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66,
    'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
    'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Sc': 1.70, 'Ti': 1.60, 'V': 1.53, 'Cr': 1.39,
    'Mn': 1.39, 'Fe': 1.32, 'Co': 1.26, 'Ni': 1.24, 'Cu': 1.32, 'Zn': 1.22, 'Ga': 1.22, 'Ge': 1.20,
    'As': 1.19, 'Se': 1.20, 'Br': 1.20,
    # Add more elements as needed
}

def bond_order_estimation(atom1, atom2, distance):
    # Define bond reduction factors (percentages) for double and triple bonds
    bond_reduction_factors = {
        'single': 1.0,
        'double': 0.91,  # e.g., double bond length is 90% of the sum of covalent radii
        'triple': 0.7,  # e.g., triple bond length is 80% of the sum of covalent radii
    }
    # Check for hydrogen involvement, enforce single bond
    tolerance_factor = 0.4  # Adjust based on your dataset for better accuracy

    # Immediate return for hydrogen involvement, but first check if they're bonded
    if 'H' in [atom1, atom2]:
        sum_radii = covalent_radii.get(atom1, 0) + covalent_radii.get(atom2, 0)
        if distance <= sum_radii + tolerance_factor:
            return 1  # If bonded, it's a single bond
        else:
            return 0  # No bond
 
    # Sum of covalent radii for a single bond
    single_bond_length = covalent_radii.get(atom1, 0) + covalent_radii.get(atom2, 0)
    
    # Calculate expected lengths for double and triple bonds
    expected_double_bond_length = single_bond_length * bond_reduction_factors['double']
    expected_triple_bond_length = single_bond_length * bond_reduction_factors['triple']

    # Estimate bond order based on actual distance
    if distance <= expected_triple_bond_length:
        return 3  # Triple bond
    elif distance <= expected_double_bond_length:
        return 2  # Double bond
    elif distance <= single_bond_length:
        return 1  # Single bond or assume single if unsure
    else:
        return 0 # no bond


def calculate_distance(coord1, coord2):
    return sqrt(sum((c2 - c1) ** 2 for c1, c2 in zip(coord1, coord2)))

def calculate_unique_distances_and_connectivity(molecules):
    results = []
    for symbols, coordinates in molecules:
        connectivity = []
        for (i1, coord1), (i2, coord2) in combinations(enumerate(coordinates), 2):
            distance = calculate_distance(coord1, coord2)
            # Ensure both symbols are in the covalent radii dictionary
            if symbols[i1] in covalent_radii and symbols[i2] in covalent_radii:
                atom1, atom2 = symbols[i1], symbols[i2]
                bond_order = bond_order_estimation(atom1, atom2, distance)
                if bond_order > 0:  # A bond is detected
                    connectivity.append([i1, i2, bond_order])
        results.append({
            "symbols": symbols,
            "geometry": [coord for triplet in coordinates for coord in triplet],
            "connectivity": connectivity
        })
    return results


def print_molecules_as_json(molecules_data):
    formatted_data = {"topologies": molecules_data}
    print(json.dumps(formatted_data, indent=2))

def read_json(filepath):
    with open(filepath, 'r') as file:
        data = json.load(file)
    molecules = []
    for topology in data.get("topologies", []):
        symbols = topology.get("symbols", [])
        geometry = topology.get("geometry", [])
        coordinates = [geometry[i:i+3] for i in range(0, len(geometry), 3)]
        molecules.append((symbols, coordinates))
    return molecules

def read_xyz(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()
    symbols = []
    coordinates = []
    for line in lines[2:]:
        parts = line.strip().split()
        symbols.append(parts[0])
        coordinates.append(list(map(float, parts[1:4])))
    return [(symbols, coordinates)]


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <file.json or file.xyz>")
        sys.exit(1)
    
    filepath = sys.argv[1]
    extension = os.path.splitext(filepath)[1]
    
    if extension == '.json':
        molecules = read_json(filepath)
    elif extension == '.xyz':
        molecules = read_xyz(filepath)
    else:
        print("Unsupported file format.")
        sys.exit(1)
    
    molecules_data = calculate_unique_distances_and_connectivity(molecules)
    print_molecules_as_json(molecules_data)

