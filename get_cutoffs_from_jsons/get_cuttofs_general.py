import json
import matplotlib.pyplot as plt
from collections import defaultdict
import sys
import os

def analyze_and_plot_nmers(file_path):
    # Read the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Split the filename from its directory and extension
    base_name = os.path.basename(file_path)  # Extracts the filename from a path
    file_name_without_ext = os.path.splitext(base_name)[0]  # Removes the file extension



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


        

        
        
    # Plotting setup for dimers
    plt.figure(figsize=(12, 8))
    #plt.scatter(dimer_distances_plot, dimer_avg_deltaEs_plot, color='blue', label='Dimers \u0394E')
    plt.axhline(y=0.1, color='black', linestyle='--', label='y=0.1')
    plt.axhline(y=-0.1, color='black', linestyle='--', label='y=-0.1')
    dimer_distances_plot = [dimer_distances[fragments] for fragments in dimer_deltaEs_print]
    dimer_abs_deltaEs_plot = [(small - dimer_deltaEs_print[fragments])*conversion_to_kj for fragments in dimer_deltaEs_print]
    plt.scatter(dimer_distances_plot, dimer_abs_deltaEs_plot, color='blue', label='Dimers')

    # Check if trimers exist and process them, only including those containing the reference fragment
    if trimers:
        trimer_results = []
        trimer_distances_deltaEs = defaultdict(list)
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



        trimer_distances_plot = [result["distance"] for result in trimer_results]
        trimer_abs_deltaEs_plot = [result["deltaE"]*conversion_to_kj for result in trimer_results]
        plt.scatter(trimer_distances_plot, trimer_abs_deltaEs_plot, color='red', label='Trimers')

        #plt.scatter(trimer_distances_plot, trimer_avg_deltaEs_plot, color='red', label='Trimers \u0394E')

    # Finalize plotting
    #plt.title(f'Log Scale Abs(Avg \u0394E) vs. Distance for Nmers Containing Fragment {reference_fragment}')
    plt.ylim(-0.5, 0.5)
    #plt.xlim(0,180)
    plt.xlabel('Distance')
    plt.ylabel('\u0394E (kJ/mol)')
    plt.legend()
    plt.grid(True)
    # Construct the new plot filename
    if trimers:
        plot_file_name = f"{file_name_without_ext}_nmers_plot.png"
    else:
        plot_file_name = f"{file_name_without_ext}_dimers_plot.png"

    plt.savefig(plot_file_name, dpi=500)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 script.py <filename>")
        sys.exit(1)

    filename = sys.argv[1]
    analyze_and_plot_nmers(filename)
