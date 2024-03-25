import json
import matplotlib.pyplot as plt

def analyze_and_sort_trimers_with_corrected_formula(file_path):
    # Read the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Extracting the nmers array under the specified structure
    monomers = data["qmmbe"]["nmers"][0]
    dimers = data["qmmbe"]["nmers"][1]
    trimers = data["qmmbe"]["nmers"][2]

    # Preprocess monomers to store their id and total energy
    monomer_energies = {}
    for monomer in monomers:
        id = monomer["fragments"][0]
        total_energy = monomer["hf_energy"] #+ monomer["mp2_os_correction"] + monomer["mp2_ss_correction"]
        monomer_energies[id] = total_energy

    # Preprocess dimers to calculate deltaE_IJ and store their distances
    dimer_deltaEs = {}
    dimer_distances = {}
    for dimer in dimers:
        fragments = tuple(sorted(dimer["fragments"]))
        E_IJ = dimer["hf_energy"] #+ dimer["mp2_os_correction"] + dimer["mp2_ss_correction"]
        E_I = monomer_energies.get(fragments[0], 0)
        E_J = monomer_energies.get(fragments[1], 0)
        deltaE_IJ = E_IJ - (E_I + E_J)
        dimer_deltaEs[fragments] = deltaE_IJ
        dimer_distances[fragments] = dimer["fragment_distance"]

    # Analyze trimers for corrected deltaE_IJK calculation with calculated distance
    trimer_results = []
    for trimer in trimers:
        fragments = tuple(sorted(trimer["fragments"]))
        if 0 in fragments:
            E_IJK = trimer["hf_energy"] #+ trimer["mp2_os_correction"] + trimer["mp2_ss_correction"]
            deltaE_components = []
            distances = []
            for i in range(3):
                for j in range(i+1, 3):
                    dimer_key = tuple(sorted([fragments[i], fragments[j]]))
                    distances.append(dimer_distances.get(dimer_key, 0))
                    deltaE_components.append(dimer_deltaEs.get(dimer_key, 0))
            E_I = monomer_energies.get(fragments[0], 0)
            E_J = monomer_energies.get(fragments[1], 0)
            E_K = monomer_energies.get(fragments[2], 0)
            # deltaE_IJK = (trimer["delta_hf_energy"] +
            #           trimer["delta_mp2_os_correction"] +
            #           trimer["delta_mp2_ss_correction"])
            deltaE_IJK = E_IJK - sum(deltaE_components) - (E_I + E_J + E_K)
            
            # if(abs(deltaE_IJK) > 10):
            #      print(trimer["mp2_os_correction"], trimer["mp2_ss_correction"],  sum(deltaE_components), E_I, E_J, E_K, fragments[0], fragments[1], fragments[2])
            max_distance = min(distances) if distances else 0
            #if(abs(deltaE_IJK) < 1E-1 ):
            trimer_results.append({"distance": max_distance, "deltaE": deltaE_IJK})

    # Sorting results by distance
    sorted_trimer_results = sorted(trimer_results, key=lambda x: x["distance"])

    # Plotting
    distances = [result["distance"] for result in sorted_trimer_results]
    abs_deltaEs = [abs(result["deltaE"]) for result in sorted_trimer_results]

    plt.figure(figsize=(10, 6))
    plt.semilogy(distances, abs_deltaEs, 'bo')
    plt.title('Log Scale Abs(DeltaE) vs. Calculated Distance for Trimers Containing Fragment 0')
    plt.xlabel('Calculated Distance')
    plt.ylabel('Abs(DeltaE) (log scale)')
    plt.grid(True)
    plt.savefig('trimer_corrected_formula_plot.png', dpi=500)

    # Filtering and printing results where deltaE < 3E-7
    # for trimer in sorted_trimer_results:
    #     if abs(trimer["deltaE"]) > 3e-1:
    #         print(f'Calculated Distance: {trimer["distance"]}, deltaE: {trimer["deltaE"]}')

    return sorted_trimer_results

# Replace 'your_file_path_here' with the actual path to your 'nmers.json' file
# Uncomment the line below to run the function and execute the final part
sorted_trimer_deltas = analyze_and_sort_trimers_with_corrected_formula('trimers.json')

