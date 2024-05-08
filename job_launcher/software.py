from queue import *
import subprocess
import os

a_no_to_symbol = {1:"H", 2:"He", 3:"Li", 4:"Be", 5:"B", 6:"C", 7:"N", 8:"O",
                  9:"F", 10:"Ne", 11:"Na", 12:"Mg", 13:"Al", 14:"Si", 15:"P",
                  16:"S", 17:"Cl", 18:"Ar", 19:"K", 20:"Ca", 21:"Sc", 22:"Ti",
                  23:"V", 24:"Cr", 25:"Mn", 26:"Fe", 27:"Co", 28:"Ni", 29:"Cu",
                  30:"Zn", 31:"Ga", 32:"Ge", 33:"As", 34:"Se", 35:"Br",
                  36:"Kr", 37:"Rb", 38:"Sr", 39:"Y", 40:"Zr", 41:"Nb", 42:"Mo",
                  43:"Tc", 44:"Ru", 45:"Rh", 46:"Pd", 47:"Ag", 48:"Cd",
                  49:"In", 50:"Sn", 51:"Sb", 52:"Te", 53:"I", 54:"Xe", 55:"Cs",
                  56:"Ba", 57:"La", 58:"Ce", 59:"Pr", 60:"Nd", 61:"Pm",
                  62:"Sm", 63:"Eu", 64:"Gd", 65:"Tb", 66:"Dy", 67:"Ho",
                  68:"Er", 69:"Tm", 70:"Yb", 71:"Lu", 72:"Hf", 73:"Ta", 74:"W",
                  75:"Re", 76:"Os", 77:"Ir", 78:"Pt", 79:"Au", 80:"Hg",
                  81:"Tl", 82:"Pb", 83:"Bi", 84:"Po", 85:"At", 86:"Rn",
                  87:"Fr", 88:"Ra", 89:"Ac", 90:"Th", 91:"Pa", 92:"U", 93:"Np",
                  94:"Pu", 95:"Am", 96:"Cm", 97:"Bk", 98:"Cf", 99:"Es",
                  100:"Fm", 101:"Md", 102:"No", 103:"Lr", 104:"Rf", 105:"Db",
                  106:"Sg", 107:"Bh", 108:"Hs", 109:"Mt", 110:"Ds", 111:"Rg"}

symbol_to_a_no = {}
for a_no in a_no_to_symbol:
    symbol_to_a_no[a_no_to_symbol[a_no]] = a_no

class Atom:
    def __init__(self, symbol, coord):
        self.symbol = symbol
        self.coord = coord

class Geometry:
    def __init__(self):
        self.atoms = []

    def from_xyz(self, xyz_file):
        raw_symbols = []
        raw_coords = []

        with open(xyz_file) as xyz:
            xyz.readline()
            xyz.readline()
            for line in xyz.readlines():
                dat = line.split()
                raw_symbols.append(dat[0])
                raw_coords.append(dat[1])
                raw_coords.append(dat[2])
                raw_coords.append(dat[3])

        self.atoms = []
        for i in range(0, len(raw_coords), 3):
            self.atoms.append(Atom(raw_symbols[i//3], [raw_coords[i], raw_coords[i+1], raw_coords[i+2]]))


class Orca:
    def __init__(self):
        self.name = "orca"

    def create_input_file(self, geometry, resource_request):
        file = ""
        file += "! RI-MP2 cc-pVDZ cc-pVDZ/C TIGHTSCF NoFrozenCore\n"
        file += f"%PAL NPROCS {resource_request.ncpus} END\n"
        file += "*xyz 0 1\n"
        for atom in geometry.atoms:
            file += f"{atom.symbol}    {atom.coord[0]}    {atom.coord[1]}    {atom.coord[2]}\n"
        file += "*\n"
        return file

    def create_run_script(self, resource_request, input_file):
        file = ""
        file += 'module load orca\n'
        file += f'input="{input_file}"\n'
        file += 'cp ${input} $PBS_JOBFS/\n'
        file += 'cd $PBS_JOBFS\n'
        file += '\n'
        file += 'CMD=`which orca`\n'
        file += '${CMD} ${input}\n'
        return file


class Qchem:
    def __init__(self):
        self.name = "qchem"

    def create_input_file(self, geometry, resource_request):
        file = ""
        file += "$molecule\n"
        file += "0 1\n"
        for atom in geometry.atoms:
            file += f"  {atom.symbol}    {atom.coord[0]}    {atom.coord[1]}    {atom.coord[2]}\n"
        file += "$end\n\n"
    
        file += "$rem\n"
        file += "  JOBTYPE      sp\n"
        file += "  EXCHANGE     hf\n"
        file += "  METHOD       rimp2\n"
        file += "  BASIS        cc-pvdz\n"
        file += "  AUX_BASIS    rimp2-cc-pvtz\n"
        file += f"  MEM_TOTAL    {resource_request.mem*1000}\n"
        file += "  N_FROZEN_CORE 0\n"
        file += "  PURECART     2222\n"
        file += "  THRESH       8\n"
        file += "$end\n"
        return file

    def create_run_script(self, resource_request, input_file):
        file = ""
        file += 'module load qchem\n'
        file += f'qchem -nt {resource_request.ncpus} {input_file}\n'
        return file

class Gamess:
    def __init__(self):
        self.name = "gamess"

    def create_input_file(self, geometry, resource_request):
        file = ""
        file += " $CONTRL SCFTYP=RHF RUNTYP=ENERGY MPLEVL=2 MAXIT=30 MULT=1 ISPHER=1 $END\n"
        file += " $SYSTEM mwords=16000 memddi=1000 $END\n"
        file += " $BASIS GBASIS=CCD $END\n"
        file += " $SCF DIRSCF=.TRUE. $END\n"
        file += "\n"
        file += " $INTGRL INTOMP=1 $END\n"
        file += " $mp2    NACORE=0 code=rimp2 $end\n"
        file += " $rimp2  othaux=.f. ivmtd=2 gosmp=.f. usedm=.true. $end\n"
        file += " $auxbas cabnam=ccd $end\n"
        file += " $DATA\n"
        file += "Title\n"
        file += "C1\n"
        for atom in geometry.atoms:
            file += f"{atom.symbol} {symbol_to_a_no[atom.symbol]}.0    {atom.coord[0]}    {atom.coord[1]}    {atom.coord[2]}\n"
        file += "$END\n"
        return file

    def create_run_script(self, resource_request, input_file):
        file = ""
        file += 'module load gamess\n'
        file += f'rungms {input_file} 00 {resource_request.ncpus}\n'
        return file


class NWchem:
    def __init__(self):
        self.name = "nwchem"

    def create_input_file(self, geometry, resource_request):
        file = ""
        file += "geometry units au\n"
        for atom in geometry.atoms:
            file += f"  {atom.symbol}    {atom.coord[0]}    {atom.coord[1]}    {atom.coord[2]}\n"
        file += "  symmetry c1\n"
        file += "end\n\n"

        mb_per_cpu = 1000*resource_request.mem/resource_request.ncpus
        file += f"memory stack {round(0.4*mb_per_cpu)} mb heap {round(0.1*mb_per_cpu)} mb global {round(0.5*mb_per_cpu)} mb\n"
    
        file += "basis\n"
        file += " * library cc-pvdz\n"
        file += "end\n\n"
        file += 'basis "ri-mp2 basis"\n'
        file += " * library cc-pvdz-ri\n"
        file += "end\n\n"
        file += "task rimp2\n"
        return file

    def create_run_script(self, resource_request, input_file):
        file = ""
        file += f"cp {input_file} $PBS_JOBFS/\n"
        file += "cd $PBS_JOBFS\n"
        file += "module unload openmpi\n"
        file += "module load nwchem\n"
        file += f"mpirun -np {resource_request.ncpus} nwchem {input_file}\n"
        return file

# JOB SYSTEMS 

class PBS:
    def create_job_script_header(self, resource_request):
        file = ""
        file += "#!/bin/bash\n"
        file += f"#PBS -q {resource_request.queue}\n"
        file += "#PBS -j oe\n"
        file += "#PBS -l wd\n"
        file += f"#PBS -l walltime={resource_request.minutes//60}:{resource_request.minutes%60:02d}:00\n"
        file += f"#PBS -l mem={resource_request.mem}GB\n"
        file += f"#PBS -l jobfs={resource_request.jobfs}GB\n"
        file += f"#PBS -l ncpus={resource_request.ncpus}\n"
        if resource_request.ngpus > 0:
            file += f"#PBS -l ngpus={resource_request.ngpus}\n"
        file += "\n"
        return file

    def launch_job(self, directory, job_filename):
        subprocess.run(['qsub', job_filename], cwd=directory)


def run_input(job_name, directory, job_system, resource_request, software, geometry, repeats=1):
    if os.path.exists(directory) and len(os.listdir(directory)) > 0:
        raise Exception(f"Target directory {directory} exists and is non-empty")

    os.makedirs(directory, exist_ok=True)

    input_filename = f"{job_name}.inp"
    job_filename = f"batch_{job_name}.sh"

    with open(f"{directory}/{input_filename}",'w') as input_file:
        input_file.write(software.create_input_file(geometry, resource_request))

    with open(f"{directory}/{job_filename}",'w') as job_file:
        job_file.write(job_system.create_job_script_header(resource_request))
        job_file.write(software.create_run_script(resource_request, input_filename))

    for _ in range(repeats):
        job_system.launch_job(directory, job_filename)




#xyz = Geometry()
#xyz.from_xyz("w40.xyz")
#resources = queues["normal"].resource_request(24)
#
#run_input("w40_test_normal", "output/w40_test_normal/", PBS(), resources, Orca(), xyz, repeats=3)
