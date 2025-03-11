import json
import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import os

def analyze_nmers(file_path):
    # Read the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Split the filename from its directory and extension
    # Extracting the reference fragment and initializing data structures
    reference_fragment = data["qmmbe"]["reference_fragment"]
    monomers = data["qmmbe"]["nmers"][0]
    dimers = data["qmmbe"]["nmers"][1]
    trimers = data["qmmbe"]["nmers"][2] if len(data["qmmbe"]["nmers"]) > 2 else None
    monomer_energies = {}
    dimer_deltaEs = {}
    dimer_distances = {}
    dimer_deltaEs_print = {}
    dimer_distances_print = {}
    conversion_to_kj=2625.5
    small = 1e-5
    less_small=9e-1
    # Process monomers
    for monomer in monomers:
        id = monomer["fragments"][0]
        total_energy = monomer["hf_energy"] + monomer["mp2_os_correction"] + monomer["mp2_ss_correction"]
        monomer_energies[id] = total_energy

    # Process dimers, only including those containing the reference fragment
    dimer_results = []
    dimer_distances_deltaEs = defaultdict(list)
    
    distances_for_cutoff_consideration = []
    for dimer in dimers:
        fragments = tuple(sorted(dimer["fragments"]))
        fragments = tuple(sorted(dimer["fragments"]))
        E_IJ = dimer["hf_energy"] + dimer["mp2_os_correction"] + dimer["mp2_ss_correction"]
        E_I = monomer_energies.get(fragments[0], 0)
        E_J = monomer_energies.get(fragments[1], 0)
        deltaE_IJ = E_IJ - (E_I + E_J)
        dimer_deltaEs[fragments] = deltaE_IJ
        dimer_distances[fragments] = dimer["fragment_distance"]
        if reference_fragment in fragments:
            E_IJ = dimer["hf_energy"] + dimer["mp2_os_correction"] + dimer["mp2_ss_correction"]
            E_I = monomer_energies.get(fragments[0], 0)
            E_J = monomer_energies.get(fragments[1], 0)
            deltaE_IJ = E_IJ - (E_I + E_J)
            dimer_distances_print[fragments] = dimer["fragment_distance"]

            dimer_deltaEs_print[fragments] = deltaE_IJ
            if conversion_to_kj * deltaE_IJ > 0.1:
                distances_for_cutoff_consideration.append(dimer["fragment_distance"])
            
    print(f"dimer cutoff: {max(distances_for_cutoff_consideration)}")     

    # Check if trimers exist and process them, only including those containing the reference fragment
    if trimers:
        trimer_results = []
        trimer_distances_deltaEs = defaultdict(list)
        trimer_distances_for_cutoff_consideration = []
        for trimer in trimers:
            fragments = tuple(sorted(trimer["fragments"]))
            if reference_fragment in fragments:
                E_IJK = trimer["hf_energy"] + trimer["mp2_os_correction"] + trimer["mp2_ss_correction"]
                deltaE_components = []
                distances = []
                for i in range(3):
                    for j in range(i+1, 3):
                        dimer_key = tuple(sorted([fragments[i], fragments[j]]))
                        if dimer_key in dimer_distances:
                            distances.append(dimer_distances[dimer_key])
                            deltaE_components.append(dimer_deltaEs[dimer_key])
                E_I = monomer_energies.get(fragments[0], 0)
                E_J = monomer_energies.get(fragments[1], 0)
                E_K = monomer_energies.get(fragments[2], 0)
                deltaE_IJK = E_IJK - sum(deltaE_components) - (E_I + E_J + E_K)
                max_distance = round(min(distances) if distances else 0, 16)

                trimer_distances_deltaEs[max_distance].append(deltaE_IJK)

                trimer_results.append({"distance": max_distance, "deltaE": deltaE_IJK})
                if deltaE_IJK * conversion_to_kj > 0.1:
                    trimer_distances_for_cutoff_consideration.append(max_distance)
        print(f"trimer cutoff: {max(trimer_distances_for_cutoff_consideration)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    analyze_nmers(filename)
