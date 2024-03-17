import os

def generate_pbs_script(input_file, output_file):
    script = f"""#!/bin/bash
#PBS -l ncpus=52
#PBS -l mem=512GB
#PBS -l jobfs=200GB
#PBS -q normalsr
#PBS -P kx58
#PBS -l walltime=01:00:00
#PBS -l wd
module load qchem
qchem -nt 52 {input_file}
"""
    with open(output_file, 'w') as f:
        f.write(script)

def main():
    input_directory = "./"
    output_directory = "./"

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    for filename in os.listdir(input_directory):
        if filename.endswith(".inp"):
            input_file = os.path.join(input_directory, filename)
            output_file = os.path.join(output_directory, f"{os.path.splitext(filename)[0]}.pbs")
            generate_pbs_script(input_file, output_file)

if __name__ == "__main__":
    main()

