import json
import argparse

# Import the analysis framework
from fraganalysis import *

def parse_arguments():
    parser = argparse.ArgumentParser(description='Fragment analysis tool')
    parser.add_argument('input', type=str, action='store', help='Path to input json')
    args = parser.parse_args()
    return args

args = parse_arguments()

dinosaur()

input_file = None
with open(args.input) as f:
    input_file = json.load(f)

# Fragmentation level is inferred from number of cutoffs
cutoffs = [1000.0, 1000.0]

# FragmentedSystem takes a topology json
system = FragmentedSystem(input_file, cutoffs)
#system = FragmentedSystem(input_file['topologies'][0], cutoffs)

ref_monomer = system.get_reference_monomer()
print("Reference fragment: ", ref_monomer)

# Form polymers by list of monomer ids
tetramer = system.get_n_mer([0,1,2,3])
tetramer.output_xyz('tetramer.xyz')

# Largest dimer/trimer
basis = 'cc-pvdz'

print("Monomers:", system.n_mer_count(1))
monomer = system.largest_n_mer(1, basis)
print("Largest monomer:", monomer.id, "with", monomer.nbas(basis), "basis functions and", monomer.n_atoms(), "atoms")
Fragment([monomer]).output_xyz('monomer.xyz')

print("Dimers:", system.n_mer_count(2))
dimer = system.largest_n_mer(2, basis)
print("Largest dimer:", dimer.monomer_ids(), "with", dimer.nbas(basis), "basis functions and", dimer.n_atoms(), "atoms")
dimer.output_xyz('dimer.xyz')

print("Trimers:", system.n_mer_count(3))
trimer = system.largest_n_mer(3, basis)
print("Largest trimer:", trimer.monomer_ids(), "with", trimer.nbas(basis), "basis functions and", trimer.n_atoms(), "atoms")
trimer.output_xyz('trimer.xyz')
