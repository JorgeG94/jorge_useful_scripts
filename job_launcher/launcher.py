from queue import *
from software import *
import argparse
import os

def launch(experiment_name, xyz_filename, software_list, job_system, resources, repeats=3):
    exp_directory = f"experiments/{experiment_name}_{resources.queue}_{resources.ncpus}"
    postfix = 1
    while os.path.exists(exp_directory) and len(os.listdir(exp_directory)) > 0:
        exp_directory = f"experiments/{experiment_name}_{resources.queue}_{resources.ncpus}_{postfix}"
        postfix += 1
    os.makedirs(exp_directory, exist_ok=True)

    xyz = Geometry()
    xyz.from_xyz(xyz_filename)

    for software in software_list:
        software_dir = f"{exp_directory}/{software.name}"
        os.mkdir(software_dir)

        run_input(experiment_name, software_dir, job_system, resources, software, xyz, repeats=repeats);

# python3 <> -i input.xyz -j pbs -q normal -c 48 -w 120

def parse_arguments():
    parser = argparse.ArgumentParser(description='Systematically run computational chemistry jobs')
    parser.add_argument('-i', required=True, type=str, action='store', help='Path to input .xyz file.')
    parser.add_argument('-j', required=True, type=str, action='store', help='Job system', choices=['pbs','slurm'])
    parser.add_argument('-q', required=True, type=str, action='store', help='Job queue')
    parser.add_argument('-c', required=True, type=int, action='store', help='# cpu cores')
    parser.add_argument('-w', type=int, action='store', help='Walltime (minutes)', default=120)
    parser.add_argument('--name', required=True, type=str, action='store', help='Unique experiment name')
    
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    resources = queues[args.q].resource_request(ncpus=args.c,minutes=args.w)

    if not os.path.isfile(args.i):
        raise Exception(f"Could not find input file {args.i}")

    job_system = None 
    if args.j == 'pbs':
        job_system = PBS()
    else:
        raise Exception("Only PBS currently supported")

    launch(args.name, args.i, [Orca(), Qchem(), Gamess()], job_system, resources)

    


if __name__ == "__main__":
    main()
