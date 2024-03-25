import json
import matplotlib.pyplot as plt

def analyze_and_sort_dimers_with_monomers(file_path):
    # Read the JSON file
    with open(file_path, 'r') as file:
        data = json.load(file)

    # Extracting the dimers and monomers array under the specified structure
    monomers = data["qmmbe"]["nmers"][0]
    dimers = data["qmmbe"]["nmers"][1]

    # Preprocess monomers to store their id and total energy
    monomer_energies = {}
    for monomer in monomers:
        id = monomer["fragments"][0]  # Assuming each monomer has a single id in the 'fragments' list
        total_energy = monomer["hf_energy"] + monomer["mp2_os_correction"] + monomer["mp2_ss_correction"]
        monomer_energies[id] = total_energy

    # Analyze dimers for deltaE calculation
    dimer_results = []
    for dimer in dimers:
        if 0 in dimer["fragments"]:
            fragments = dimer["fragments"]
            dimer_energy = dimer["hf_energy"] + dimer["mp2_os_correction"] + dimer["mp2_ss_correction"]
            if all(frag in monomer_energies for frag in fragments):  # Ensure both fragments have preprocessed energies
                E_IJ = dimer_energy
                E_I = monomer_energies[fragments[0]]
                E_J = monomer_energies[fragments[1]]
                deltaE = E_IJ - (E_I + E_J)
                dimer_results.append({"distance": dimer["fragment_distance"], "deltaE": deltaE})
    print(len(dimers))
    # Sorting results by distance
    sorted_dimer_results = sorted(dimer_results, key=lambda x: x["distance"])
    #print(sorted_dimer_results)

    # Iterate through sorted_dimer_results and print when deltaE < 3E-7
    # Plotting deltaE vs. distance
    distances = [result["distance"] for result in sorted_dimer_results]
    abs_deltaEs = [abs(result["deltaE"]) for result in sorted_dimer_results]

    plt.figure(figsize=(10, 6))
    plt.semilogy(distances, abs_deltaEs, 'bo')  # 'bo' for blue circle markers
    plt.title('Log Scale Abs(DeltaE) vs. Distance for Dimers Containing Fragment 0')
    plt.xlabel('Distance')
    plt.ylabel('Abs(DeltaE) (log scale)')
    plt.grid(True)


    #plt.plot(distances, deltaEs, marker='o', linestyle='-')
    #plt.title('DeltaE vs. Distance for Dimers Containing Fragment 0')
    #plt.xlabel('Distance')
    #plt.ylabel('DeltaE')
    #plt.grid(True)
    plt.savefig('output_plot.png', dpi=500)
    for dimer in sorted_dimer_results:
        if abs(dimer["deltaE"]) < 3e-7:
            print(f'Distance: {dimer["distance"]}, deltaE: {dimer["deltaE"]}')

    return sorted_dimer_results

# Replace 'your_file_path_here' with the actual path to your 'nmers.json' file
# Uncomment the line below to run the function and execute the final part
sorted_dimer_deltas = analyze_and_sort_dimers_with_monomers('dimers.json')

