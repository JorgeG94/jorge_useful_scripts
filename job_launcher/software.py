from queue import *
import subprocess
import os

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

    input_filename = f"{job_name}.in"
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
