import os

template_content = """$molecule
0 1

{coordinates}
$end

$rem
   JOBTYPE OPT
   METHOD        RIMP2
   BASIS         cc-pvdz
   AUX_BASIS     rimp2-cc-pvdz
   PURECART 2222
   N_FROZEN_CORE 0
   MEM_TOTAL 10000
$end
"""

# Iterate over each .xyz file in the directory
for filename in os.listdir("./xyzs"):
    if filename.endswith(".xyz"):
        xyz_file_path = os.path.join("./xyzs", filename)

        # Read the coordinates from the .xyz file
        with open(xyz_file_path, "r") as xyz_file:
            lines = xyz_file.readlines()
            coordinates = "".join(lines[2:])

        # Create a unique output file based on the .xyz file name
        output_file = filename.replace(".xyz", "_opt.inp")

        # Write the final content to the output .inp file
        with open(output_file, "w") as final_inp_file:
            final_inp_file.write(template_content.format(coordinates=coordinates))

