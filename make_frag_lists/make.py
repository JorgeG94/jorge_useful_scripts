def generate_fragments(nfrag, n_atom, n_atom_reference):
    # Calculate the total numbers needed to fill all fragments
    total_numbers = nfrag * n_atom + n_atom_reference
    
    # Check if the total_numbers exceeds the maximum allowed value
    if total_numbers > nfrag * n_atom + n_atom_reference:
        raise ValueError("The total number of elements exceeds the maximum allowed value of 126*n_atom.")
    
    # Initialize the fragments dictionary
    fragments = {"fragments": []}
    
    # Generate the numbers for each fragment
    for i in range(nfrag):
        if(i == 0):
                start = 0
                end = start + n_atom_reference
                fragment = list(range(start,end))
                fragments["fragments"].append(fragment)
                                
        start = i * n_atom + n_atom_reference
        end = start + n_atom
        fragment = list(range(start, end))
        fragments["fragments"].append(fragment)
    
    return fragments
def generate_fragment_charges(nfrag):
		fragment_charges = {"fragment_formal_charges": []}
		for i in range(nfrag):
				fragment_charges_l = 0
				fragment_charges["fragment_formal_charges"].append(fragment_charges_l)

		return fragment_charges
# Example usage
nfrag = 18  # Number of fragments
n_atom = 20  # Size of each fragment
n_atom_reference = 0
fragments = generate_fragments(nfrag, n_atom, n_atom_reference)
print(fragments)
fragment_charges = generate_fragment_charges(nfrag+1)
print(fragment_charges)

