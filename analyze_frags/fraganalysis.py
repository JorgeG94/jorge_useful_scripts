#!/bin/python3
from pyscf import gto
from itertools import combinations

########################################################################################################
############################################# INFRASTRUCTURE ###########################################
########################################################################################################

a_no_to_symbol = {1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O", 9:"F", 10:"Ne",
    11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P", 16:"S", 17:"Cl", 18:"Ar", 19:"K", 20:"Ca",
    21:"Sc", 22:"Ti", 23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu", 30:"Zn",
    31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br", 36:"Kr", 37:"Rb", 38:"Sr", 39:"Y", 40:"Zr",
    41:"Nb", 42:"Mo", 43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd", 49:"In", 50:"Sn",
    51:"Sb", 52:"Te", 53:"I", 54:"Xe", 55:"Cs", 56:"Ba", 57:"La", 58:"Ce", 59:"Pr", 60:"Nd",
    61:"Pm", 62:"Sm", 63:"Eu", 64:"Gd", 65:"Tb", 66:"Dy", 67:"Ho", 68:"Er", 69:"Tm", 70:"Yb",
    71:"Lu", 72:"Hf", 73:"Ta", 74:"W", 75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au", 80:"Hg",
    81:"Tl", 82:"Pb", 83:"Bi", 84:"Po", 85:"At", 86:"Rn", 87:"Fr", 88:"Ra", 89:"Ac", 90:"Th",
    91:"Pa", 92:"U", 93:"Np", 94:"Pu", 95:"Am", 96:"Cm", 97:"Bk", 98:"Cf", 99:"Es", 100:"Fm",
    101:"Md", 102:"No", 103:"Lr", 104:"Rf", 105:"Db", 106:"Sg", 107:"Bh", 108:"Hs", 109:"Mt", 110:"Ds",
    111:"Rg"}

symbol_to_a_no = {}
for a_no in a_no_to_symbol:
    symbol_to_a_no[a_no_to_symbol[a_no]] = a_no


def cart(x,y):
    return ((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)**0.5
    
class Atom:
    def __init__(self, coord, a_no, atom_id):
        self.coord = coord
        self.a_no = a_no   
        self.id = atom_id
        
class Monomer:
    def __init__(self, atoms, broken_bonds = [], charge = 0, ix = -1):
        self.id = ix
        self.atoms = []
        self.charge = charge
        self.broken_bonds = broken_bonds
        self.centroid = [0.0,0.0,0.0]
        self.basis_cache = {}
        for atom in atoms:
            self.atoms.append(atom)
            for i in range(3):
                self.centroid[i] += atom.coord[i]

        for i in range(3):
            self.centroid[i] /= len(self.atoms)


    def n_atoms(self, caps=False):
        if caps:
            return len(self.atoms) + len(self.broken_bonds)
        else:
            return len(self.atoms)
            
    def fingerprint(self):
        elements = {}
        for atom in self.atoms:
            symbol = a_no_to_symbol[atom.a_no]
            if symbol not in elements:
                elements[symbol] = 0
            elements[symbol] += 1
            
        el_list = list(elements.items())
        el_list.sort(key=lambda x:x[0])
        
        ret = ""
        for el in el_list:
            ret += el[0]
            ret += str(el[1])
        return ret

    def nbas(self, basis_set):
        if basis_set in self.basis_cache:
            return self.basis_cache[basis_set]

        mol_str = ""
        for atom in self.atoms:
            atom_str = f'{a_no_to_symbol[atom.a_no]} {atom.coord[0]} {atom.coord[1]} {atom.coord[2]}; '
            mol_str += atom_str

        for bond in self.broken_bonds:
            atom_str = f'1 {bond.coord[0]} {bond.coord[1]} {bond.coord[2]}; '
            mol_str += atom_str

        mol = gto.Mole()
        mol.charge = self.charge
        mol.build(atom = mol_str, basis = basis_set)
        self.basis_cache[basis_set] = mol.nao_nr()
        return self.basis_cache[basis_set]

H2 = Monomer([Atom([0.0,0.0,0.0],1,-1), Atom([1.0,0.0,0.0],1,-1)])

class Fragment:
    def __init__(self, monomers):
        self.monomers = monomers

        self.centroid = [0.0,0.0,0.0]
        n_atoms = 0
        for frag in self.monomers:
            n_atoms += frag.n_atoms()
            for i in range(3):
                self.centroid[i] += frag.centroid[i]*frag.n_atoms()

        for i in range(3):
            self.centroid[i] /= n_atoms

        self.restored_bonds = 0
        self.h_caps = []
        atom_ids = [atom.id for atom in self.atoms()]
        for frag in self.monomers:
            for cap in frag.broken_bonds:
                if cap.id in atom_ids:
                    self.restored_bonds += 1
                else:
                    self.h_caps.append(cap)
        self.restored_bonds = self.restored_bonds//2

    def nbas(self, basis_set):
        ret = 0 
        for frag in self.monomers:
            ret += frag.nbas(basis_set)

        ret -= H2.nbas(basis_set)*self.restored_bonds
        return ret

    def monomer_ids(self):
        return sorted([mon.id for mon in self.monomers])

    def n_atoms(self, caps=False):
        n = sum([mon.n_atoms() for mon in self.monomers])
        if caps:
            n += len(self.h_caps)
        return n

    def atoms(self):
        ret = []
        for frag in self.monomers:
            ret.extend(frag.atoms)
        return ret

    def output_xyz(self, filename):
        with open(filename, 'w') as f:
            print(len(self.atoms())+len(self.h_caps), file=f)
            print('calum was here', file=f)
            for atom in self.atoms():
                print(f'{a_no_to_symbol[atom.a_no]} {atom.coord[0]} {atom.coord[1]} {atom.coord[2]}', file=f)
            for cap in self.h_caps:
                print(f'H {cap.coord[0]} {cap.coord[1]} {cap.coord[2]}', file=f)

# Represents topology component of input json
class Topology:
    def __init__(self):
        self.atoms = []
        self.connectivity = []
        self.fragments = []
        self.fragment_formal_charges = []
        self.xyz = None
        
    def assemble_fragments(self):
        connectivity_map = {}
        for atom in self.atoms:
            connectivity_map[atom.id] = set()

        for bond in self.connectivity:
            connectivity_map[bond[0]].add(bond[1])
            connectivity_map[bond[1]].add(bond[0])


        ret = []
        for (i,(atom_ixs,charge)) in enumerate(zip(self.fragments, self.fragment_formal_charges)):
            broken_bonds = set.union(*[connectivity_map[ix] for ix in atom_ixs])
            broken_bonds = broken_bonds.difference(set(atom_ixs)) 
            broken_bonds = [self.atoms[ix] for ix in broken_bonds]
            frag = Monomer([self.atoms[ix] for ix in atom_ixs], broken_bonds, charge, ix=i)
            ret.append(frag)
        return ret
        

    def from_json(self, topology_json):

        symbols = []
        raw_coords = []
        if 'symbols' in topology_json and 'geometry' in topology_json:
            symbols = topology_json['symbols']
            raw_coords = topology_json['geometry']


        elif 'xyz' in topology_json:
            self.xyz = topology_json['xyz']
            with open(self.xyz) as xyz:
                xyz.readline()
                xyz.readline()
                for line in xyz.readlines():
                    dat = line.split()
                    symbols.append(dat[0])
                    raw_coords.append(float(dat[1]))
                    raw_coords.append(float(dat[2]))
                    raw_coords.append(float(dat[3]))
        else:
            raise Exception("Bad input")

        coords = []
        for i in range(0, len(raw_coords), 3):
            coords.append([raw_coords[i],raw_coords[i+1],raw_coords[i+2]])

        for (i, (symbol, coord)) in enumerate(zip(symbols,coords)):
            self.atoms.append(Atom(coord, symbol_to_a_no[symbol], i))

        if 'connectivity' in topology_json:
            self.connectivity = topology_json['connectivity']
        else:
            self.connectivity = []

        self.fragments = topology_json['fragments']
        self.fragment_formal_charges = topology_json["fragment_formal_charges"]

        
    def to_json(self):
        topology_json = {}
        topology_json['fragments'] = self.fragments
        topology_json['fragment_formal_charges'] = self.fragment_formal_charges
        if self.connectivity:
            topology_json['connectivity'] = self.connectivity
        
        if self.xyz == None:
            topology_json['geometry'] = []
            topology_json['symbols'] = []
            for atom in self.atoms:
                topology_json['symbols'].append(a_no_to_symbol[atom.a_no])
                topology_json['geometry'].extend(atom.coord)
        else:
            topology_json['xyz'] = self.xyz
        return topology_json

    def nfrag(self):
        return len(self.fragments)
    
    def natoms(self):
        return len(self.atoms)

class FragmentedSystem:
    def __init__(self, topology_json, cutoffs):
        self.topology = Topology()
        self.topology.from_json(topology_json)

        self.polymers = [None, self.topology.assemble_fragments()]
        monomers = self.polymers[1]

        for (i, cutoff) in enumerate(cutoffs):
            n = i+2
            n_mers = []

            neighbour_map = []
            for (i,A) in enumerate(monomers):
                neighbours = set()
                for (j,B) in enumerate(monomers):
                    if i == j:
                        continue

                    if cart(A.centroid, B.centroid) < cutoff:
                        neighbours.add(j)

                neighbour_map.append(neighbours)

            for (i,mon) in enumerate(monomers):
                for nbrs in combinations(neighbour_map[i], n-1):
                    good = True
                    for (a,b) in combinations(nbrs, 2):
                        if cart(monomers[a].centroid, monomers[b].centroid) >= cutoff:
                            good = False
                            break

                    if not good:
                        continue

                    mons = [monomers[k] for k in nbrs]
                    mons.append(mon)
                    n_mers.append(Fragment(mons))

            self.polymers.append(n_mers)

    def get_n_mer(self, monomers):
        return Fragment([self.polymers[1][i] for i in monomers])

    def n_mer_count(self, n):
        return len(self.polymers[n])

    def largest_n_mer(self, n, basis_set):
        max_ix = 0
        max_nbas = 0
        for ix,frag in enumerate(self.polymers[n]):
            nbas = frag.nbas(basis_set)
            if nbas > max_nbas:
                max_nbas = nbas
                max_ix = ix

        return self.polymers[n][max_ix]

def dinosaur():
    print('''                                                     ___._       ''')
    print('''                                                   .'  <0>'-.._''')
    print('''                                                  /  /.--.____")''')
    print('''                                                 |   \\   __.-'~''')
    print('''                                                 |  :  -'/''')
    print('''                                                /:.  :.-' ''')
    print('''__________                                     | : '. |''')
    print(''''--.____  '--------.______       _.----.-----./      :/''')
    print('''        '--.__            `'----/       '-.      __ :/''')
    print('''              '-.___           :           \\   .'  )/''')
    print('''                    '---._           _.-'   ] /  _/''')
    print('''                         '-._      _/     _/ / _/''')
    print('''                             \\_ .-'____.-'__< |  \\___''')
    print('''                               <_______.\\    \\_\\_---.7''')
    print('''                              |   /'=r_.-'     _\\\\ =/''')
    print('''                          .--'   /            ._/'>''')
    print('''                        .'   _.-' ''')
    print('''                       / .--' ''')
    print('''                      /,/''')
    print('''                      |/`)''')
    print('''                      'c=,''')
    print()
