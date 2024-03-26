import argparse
import json

def generate_fragments(nfrag, n_atom, n_atom_reference):
    fragments = {"fragments": []}
    for i in range(nfrag):
        if i == 0 and n_atom_reference > 0:
            start = 0
            end = start + n_atom_reference
            fragment = list(range(start, end))
            fragments["fragments"].append(fragment)
            
        start = i * n_atom + n_atom_reference
        end = start + n_atom
        fragment = list(range(start, end))
        fragments["fragments"].append(fragment)
    
    return fragments

def generate_fragment_charges(nfrag):
    fragment_charges = {"fragment_formal_charges": []}
    for i in range(nfrag):
        fragment_charges["fragment_formal_charges"].append(0)

    return fragment_charges

def read_json_input(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def check_and_append(data, fragments, fragment_charges, index_of_topology, overwrite):
    if "topologies" in data and index_of_topology < len(data["topologies"]):
        topology = data["topologies"][index_of_topology]
        if not overwrite and ("fragments" in topology or "fragment_formal_charges" in topology):
            raise ValueError("Existing 'fragments' or 'fragment_formal_charges' found. Use --overwrite to replace.")
        
        topology["fragments"] = fragments["fragments"]
        topology["fragment_formal_charges"] = fragment_charges["fragment_formal_charges"]
    else:
        print(f"Index {index_of_topology} out of range in 'topologies'.")
    return data

def write_json_output(file_path, data):
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)

def main():
    parser = argparse.ArgumentParser(description="Generate fragment lists and charges, with options to update an existing JSON.")
    parser.add_argument("--natoms-per-frag", type=int, required=True, help="Number of atoms per fragment")
    parser.add_argument("--nfrags", type=int, required=True, help="Number of fragments")
    parser.add_argument("--natom-reference", type=int, default=0, help="Number of atoms in the reference fragment")
    parser.add_argument("--input-file", type=str, help="Path to input JSON file to append data")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing fragments and charges if present")

    args = parser.parse_args()

    fragments = generate_fragments(args.nfrags, args.natoms_per_frag, args.natom_reference)
    fragment_charges = generate_fragment_charges(args.nfrags)

    if args.input_file:
        data = read_json_input(args.input_file)
        index_of_topology = 0 # Assuming you want to append to the first topology; modify as needed
        try:
            data = check_and_append(data, fragments, fragment_charges, index_of_topology, args.overwrite)
            write_json_output(args.input_file, data)
            print("Updated JSON data written to", args.input_file)
        except ValueError as e:
            print(e)
    else:
        print(fragments)
        print(fragment_charges)

if __name__ == "__main__":
    main()
