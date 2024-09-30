from openbabel import openbabel as ob
import sys
import networkx as nx
from pathlib import Path
import json

def read_xyz(xyz_file):
    with open(xyz_file, 'r+') as f:
        data = f.read()

    return data

def topo_information(xyz_contents):
    lines = xyz_contents.split("\n")
    lines = lines[2:]

    geometry = []
    symbols = []
    for line in lines:
        if len(line.split()) != 4:
            continue
        atom, x, y, z = line.split()
        symbols.append(atom)
        geometry.extend([float(x), float(y), float(z)])
    
    return symbols, geometry


def main():
# create an OBMol object
    mol = ob.OBMol()

    # read the XYZ string and set the coordinates in the OBMol object
    obConversion = ob.OBConversion()
    obConversion.SetInFormat("xyz")
    
    if len(sys.argv) < 2:
        raise SystemExit("Provide .xyz file")
    xyz_file = Path(sys.argv[1])
    xyz_name = xyz_file.stem

    xyz_string = read_xyz(xyz_file)
    symbols, geometry = topo_information(xyz_string)

    obConversion.ReadString(mol, xyz_string)

    # perceive the bonds in the molecule
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()

    # num_molecules = mol.NumConformers()

    
    graph = nx.Graph()

    # iterate through the bonds and print their order
    connectivity = []
    for bond in ob.OBMolBondIter(mol):
        bond_order = bond.GetBondOrder()
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        graph.add_edge(begin_idx, end_idx, weight=3)
        connectivity.append([begin_idx, end_idx, bond_order])
        
    #     print(f"Atom indices: {begin_idx}, {end_idx}, Bond order: {bond_order}")

    fragment_list = []
    connected_components = [graph.subgraph(x) for x in nx.connected_components(graph)]
    nfrags = len(connected_components)
    fragment_formal_charges = [0 for _ in range(nfrags)]

    print(f"Number of molecules: {nfrags}")

    for connected_component in connected_components:
        fragment_nodes = list(connected_component.nodes)
        fragment_list.append(fragment_nodes)

    topology = {
        "symbols": symbols,
        "geometry": geometry,
        "connectivity": connectivity,
        "fragments": fragment_list,
        "fragment_formal_charges": fragment_formal_charges,
    }
    sys.stderr.write("WARNING: fragment formal charges are set to 0 by default.\n")

    with open(f"{xyz_name}.qdxf", "w") as f:
        json.dump(topology, f, indent=4)


if __name__ == "__main__":
    main()