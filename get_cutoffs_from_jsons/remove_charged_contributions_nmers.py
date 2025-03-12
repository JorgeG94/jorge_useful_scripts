import json 
import sys

def read_input_json(json_file):
    with open(json_file, "r") as f:
        data = json.load(f)

    protein_topo = data["topologies"][0]
    reference_frag = data["keywords"]["frag"]["reference_fragment"]

    charged_frags = set()
    charged_frags = set()
    for i, x in enumerate(protein_topo["fragment_formal_charges"]):
        if x != 0:
            charged_frags.add(i)
    return charged_frags, reference_frag

def generate_neutral_output_only_nmer(json_file, nmer, reference_frag, charged_frags):
    with open(json_file, "r") as f:
        data = json.load(f)

    new_qmmbe_nmers = []
    for x in data["qmmbe"]["nmers"][nmer-1]:
        fragments = x["fragments"]
        other_frags = [y for y in fragments if y != int(reference_frag)]

        accept = True
        for frag in other_frags:
            if frag in charged_frags:
                accept = False
                break
        
        if accept:
            new_qmmbe_nmers.append(x)
    
    data["qmmbe"]["nmers"][nmer-1] = new_qmmbe_nmers
    name = json_file.split(".json")[0]
    with open(f"{name}_filtered_{nmer}.json", "w") as g:
        json.dump(data, g, indent=4)
    
if __name__ == "__main__":
    input_json_file = sys.argv[1]
    output_json_file = sys.argv[2]
    nmer = int(sys.argv[3])

    charged_frags, reference_frag = read_input_json(input_json_file)
    generate_neutral_output_only_nmer(output_json_file, nmer, reference_frag, charged_frags)