import json
import sys
from collections import Counter

def read_json_and_print_xyz_and_electrons(json_file):
    # Mapping of element symbols to their atomic numbers
    atomic_numbers = {
        'H': 1, 'He': 2, 'Li': 3, 'Be': 4, 'B': 5, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Ne': 10,
        'Na': 11, 'Mg': 12, 'Al': 13, 'Si': 14, 'P': 15, 'S': 16, 'Cl': 17, 'Ar': 18, 'K': 19, 'Ca': 20,
        'Sc': 21, 'Ti': 22, 'V': 23, 'Cr': 24, 'Mn': 25, 'Fe': 26, 'Co': 27, 'Ni': 28, 'Cu': 29, 'Zn': 30,
        # Add more elements as needed
    }

    with open(json_file, 'r') as file:
        data = json.load(file)

    for topology in data.get("topologies", []):
        symbols = topology.get("symbols", [])
        geometry = topology.get("geometry", [])
        
        print(len(symbols))
        print("Comment line")
        for symbol, (x, y, z) in zip(symbols, zip(*[iter(geometry)]*3)):
            print(f"{symbol} {x} {y} {z}")
        
        # Calculate and print the chemical formula with spaces
        formula_counter = Counter(symbols)
        formula_with_spaces = ' '.join(f"{atom}{count if count > 1 else ''}" for atom, count in formula_counter.items())
        print("Chemical formula:", formula_with_spaces)

        # Calculate the total number of electrons
        total_electrons = sum(atomic_numbers[atom] * count for atom, count in formula_counter.items())
        print("Total number of electrons:", total_electrons)

        print("-" * 20)  # Separator for multiple topologies
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py file.json")
        sys.exit(1)
    
    json_file = sys.argv[1]
    read_json_and_print_xyz_and_electrons(json_file)
