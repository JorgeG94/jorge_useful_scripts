from pathlib import Path
import json
import re
import os

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
            protein_geometry = []
            protein_symbols = []
            ligand_geometry = []
            ligand_symbols = []
            continue
        elif line.startswith("TER") and start:
            start = False
            yield {"symbols": complex_symbols, "geometry": complex_geometry}, {"symbols": protein_symbols, "geometry": protein_geometry}, {"symbols": ligand_symbols, "geometry": ligand_geometry}
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



                ligand_symbols.append(symbol)
                ligand_geometry.extend([x, y, z])
            else:
                symbol = line_split[-1].rstrip("\n")
                protein_symbols.append(symbol)
                protein_geometry.extend([x, y, z])
            
            if symbol not in valid_symbols:
                # print(f"{symbol} is not a valid symbol.\nProblematic line: {line}", file=open("mdpdb2qdxf_error.dat", "a"))
                raise SystemExit(f"{symbol} is not a valid symbol.\nProblematic line: {line}")

            complex_symbols.append(symbol)
            complex_geometry.extend([x, y, z])

if __name__ == "__main__":
    wd = Path(os.getcwd())
    md_path = wd / "md"
    ligands = [item for item in md_path.iterdir() if item.is_dir()]

    for ligand in ligands:
        ligand_name = str(ligand.stem)
        ligand_path = wd / "md" / ligand_name
        protein_topo_file = ligand_path / f"topol.{ligand_name}.top"
        ligand_topo_file = ligand_path / f"{ligand_name}_GMX.itp"

        protein_charges = read_topology(protein_topo_file)
        ligand_charges = read_topology(ligand_topo_file)

        complex_charges = protein_charges + ligand_charges

        for i in range(2, 4):
            output_complex_topo_models = []
            output_protein_topo_models = []
            output_ligand_topo_models = []

            output_dir = wd / "pbsa" / ligand_name / f"{i}"
            output_dir.mkdir(parents=True, exist_ok=True)

            repeat_path = ligand_path / f"{i}"
            pdb_file = repeat_path / f"md.{ligand_name}.cluster.pdb"

            complex_output_path = output_dir / "complex.qdx.json"
            protein_output_path = output_dir / "protein.qdx.json"
            ligand_output_path = output_dir / "ligand.qdx.json"

            if complex_output_path.is_file() and protein_output_path.is_file() and ligand_output_path.is_file():
                continue

            for complex_topo, protein_topo, ligand_topo in read_pdb(pdb_file):
                complex_topo["partial_charges"] = complex_charges
                protein_topo["partial_charges"] = protein_charges
                ligand_topo["partial_charges"] = ligand_charges

                output_complex_topo_models.append({"topology": complex_topo})
                output_protein_topo_models.append({"topology": protein_topo})
                output_ligand_topo_models.append({"topology": ligand_topo})
                        
            with open(output_dir / "complex.qdx.json", "w") as f:
                json.dump(output_complex_topo_models, f, indent=4)

            with open(output_dir / "protein.qdx.json", "w") as f:
                json.dump(output_protein_topo_models, f, indent=4)

            with open(output_dir / "ligand.qdx.json", "w") as f:
                json.dump(output_ligand_topo_models, f, indent=4)

            
            





