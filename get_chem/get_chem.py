import json
import sys
from collections import Counter
import argparse

def read_json(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)
    return data.get("topologies", [])

def read_xyz(xyz_file):
    topologies = []
    with open(xyz_file, 'r') as file:
        lines = file.readlines()
    # Assuming XYZ file format doesn't include 'topologies' but is a single molecule
    symbols = []
    geometry = []
    for line in lines[2:]:  # Skip the first two lines
        parts = line.split()
        symbols.append(parts[0])
        geometry.extend(map(float, parts[1:4]))
    topologies.append({"symbols": symbols, "geometry": geometry})
    return topologies

def calculate_and_print_electrons_and_formula(symbols, atomic_numbers):
    formula_counter = Counter(symbols)
    formula_with_spaces = ' '.join(f"{atom}{count if count > 1 else ''}" for atom, count in formula_counter.items())
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'MG': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
        'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'MN': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        # Add more elements as needed
    }

    print("Chemical formula:", formula_with_spaces)
    print("Total number of atoms: ", len(symbols))

    total_electrons = sum(atomic_numbers[atom] * count for atom, count in formula_counter.items())
    print("Total number of electrons:", total_electrons)

def print_xyz(symbols, geometry):
    print(len(symbols))
    print("Comment line")
    for symbol, (x, y, z) in zip(symbols, zip(*[iter(geometry)]*3)):
        print(f"{symbol} {x} {y} {z}")

def main(file_path, print_xyz_flag):
    atomic_numbers = {
        # Your mapping of element symbols to their atomic numbers
    }
    
    if file_path.endswith('.json'):
        topologies = read_json(file_path)
        if print_xyz_flag:
            for topology in topologies:
                symbols = topology.get("symbols", [])
                geometry = topology.get("geometry", [])
                print_xyz(symbols, geometry)
                print("-" * 20)  # Separator for multiple topologies or molecules
    elif file_path.endswith('.xyz'):
        topologies = read_xyz(file_path)
    else:
        print("Unsupported file format.")
        sys.exit(1)

    for topology in topologies:
        symbols = topology.get("symbols", [])
        geometry = topology.get("geometry", [])
        calculate_and_print_electrons_and_formula(symbols, atomic_numbers)
        print("-" * 20)  # Separator for multiple topologies or molecules

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process a .json or .xyz file.")
    parser.add_argument("file_path", type=str, help="Path to the file to be processed.")
    parser.add_argument("--print-xyz", action="store_true", help="Print the XYZ format if the file is .json")
    
    args = parser.parse_args()

    main(args.file_path, args.print_xyz)

