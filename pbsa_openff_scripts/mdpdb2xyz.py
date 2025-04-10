from pathlib import Path
import json
import re
import os
import sys

valid_symbols = {
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr"
}

ligand_carbon_atoms = {'CAF', 'CAP', 'CAL', 'CAU', 'CAR', 'CAK', 'CAE', 'CAS', 'CAD', 'CAQ', 'CAG', 'CAM', 'CAC', 'CAI', 'CAH', 'CAY', 'CAB', 'CAX', 'CBA', 'CBC', 'CAW', 'C23', 'CAZ', 'CBE', 'CAO', 'CBB', 'CAT', 'CAV', 'CAA', 'CAJ'}
ligand_nitrogen_atoms = {'NBD', 'NAL', 'NAP', 'NAK', 'NAV', 'NAO', 'NAN', 'NAM'}
ligand_oxygen_atoms = {'OAG', 'OAQ'}


def read_topology(topo_file):
    # returns charges

    start = False
    charges = []

    with open(topo_file, "r") as f:
        contents = f.readlines()

    for line in contents:
        if "[ atoms ]" in line:
            start = True
            continue
        elif start and len(line.split()) >= 8 and not line.startswith(";"):
            charge = float(line.split()[6])
            charges.append(charge)
        elif "[ bonds ]" in line:
            start = False
            break
    return charges    

def read_pdb(pdb_file):
    start = False

    with open(pdb_file, "r") as f:
        contents = f.readlines()
    
    for line in contents:
        if line.startswith("MODEL"):
            start = True
            
            complex_geometry = []
            complex_symbols = []
            continue
        elif line.startswith("TER") and start:
            start = False
            yield {"symbols": complex_symbols, "geometry": complex_geometry}
            continue
        elif start and line.startswith("ATOM"):
            line_split = line.split()
            residue = line_split[3]
            atom_name = line_split[2]

            x = float(line_split[6])
            y = float(line_split[7])
            z = float(line_split[8])

            if residue == "UNL":
                symbol = re.sub(r'\d+', '', atom_name)
                if len(symbol) == 2 and symbol not in valid_symbols:
                    print(symbol)
                    if symbol == "OG":
                        symbol = "O"
                    elif symbol == "CB" or symbol == "CA":
                        symbol = "C"
                    elif symbol == "HA" or symbol == "HB":
                        symbol = "H"
                    else:
                        symbol = f"{symbol[0]}{symbol[1].lower()}"

                elif symbol == "OXT":
                    symbol = "O"
                elif symbol == "SAR":
                    symbol = "S"
                elif symbol in {"CLAB", "CLA", "CLAH"}:
                    symbol = "Cl"
                elif symbol in ligand_carbon_atoms:
                    symbol = "C"
                elif symbol in ligand_oxygen_atoms:
                    symbol = "O"
                elif symbol in ligand_nitrogen_atoms:
                    symbol = "N"

            else:
                symbol = line_split[-1].rstrip("\n")
            
            if symbol not in valid_symbols:
                # print(f"{symbol} is not a valid symbol.\nProblematic line: {line}", file=open("mdpdb2qdxf_error.dat", "a"))
                raise SystemExit(f"{symbol} is not a valid symbol.\nProblematic line: {line}")

            complex_symbols.append(symbol)
            complex_geometry.append([x, y, z])

def print_json_to_xyz(topo, frame, name):
    symbols = topo["symbols"]
    geometry = topo["geometry"]

    natoms = len(symbols)
    xyz_string = f"{natoms}\n\n"
    for i, symbol in enumerate(symbols):
        coordinate = geometry[i]
        x = coordinate[0]
        y = coordinate[1]
        z = coordinate[2]
        xyz_string += f"{symbol:<4}{x:<10.3f}{y:<10.3f}{z:<10.3f}"
        if i != natoms - 1:
            xyz_string += "\n"

    with open(f"{name}_{frame}.xyz", "w") as f:
        f.write(xyz_string)

if __name__ == "__main__":
    wd = Path(os.getcwd())
    pdb_file = sys.argv[1]
    name = pdb_file.split(".pdb")[0]

    frame = 0
    for complex_topo in read_pdb(pdb_file):
        print_json_to_xyz(complex_topo, frame, name)
        frame += 1
