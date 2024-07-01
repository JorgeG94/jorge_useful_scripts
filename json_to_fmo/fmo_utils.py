def make_fmoxyz(symbols, element_numbers, geometry):
    # Ensure the lengths of symbols, element_numbers, and geometry are compatible
    if len(symbols) != len(element_numbers) or len(geometry) % 3 != 0:
        raise ValueError("Invalid lengths of input lists")

    # Create the fmoxyz structure
    fmoxyz_structure = ""
    for i in range(0, len(geometry), 3):
        symbol = symbols[i // 3]
        element_number = element_numbers[i // 3]
        x, y, z = geometry[i:i + 3]

        # Append the formatted line to the structure
        fmoxyz_structure += f"{symbol} {element_number} {x} {y} {z}"

        # Add a newline character if it's not the last element
        if i < len(geometry) - 3:
            fmoxyz_structure += "\n"

    return fmoxyz_structure

def make_data(symbols, element_numbers):
    unique_atoms = set(zip(symbols, element_numbers))

    # Create a formatted string for the unique atom information
    unique_atoms_string = ""
    for i, (symbol, element_number) in enumerate(unique_atoms):
        unique_atoms_string += f"{symbol}   {element_number}"

        # Add a newline character if it's not the last element
        if i < len(unique_atoms) - 1:
            unique_atoms_string += "\n"

    return unique_atoms_string

def make_icharg(charge_values):
    # Create a formatted string for the icharg information
    charge_information_string = "icharg(1) ="
    for charge in charge_values:
        charge_information_string += f" {charge},"
    charge_information_string = charge_information_string.rstrip(",")

    return charge_information_string

def make_indat(fragments_list, atoms_per_line=10):
    # Create a mapping of atom indices to fragment indices
    atom_to_fragment_map = sorted([(atom + 1, fragment_index) for fragment_index, fragment in enumerate(fragments_list, start=1) for atom in fragment], key=lambda x: x[0])

    # Create a formatted string for indat information
    indat_information_string = f"nfrag = {len(fragments_list)}\nindat(1) ="

    for index, (atom_index, fragment_index) in enumerate(atom_to_fragment_map, start=1):
        # Print the fragment index for each atom
        indat_information_string += f" {fragment_index},"

        # Break the line every atoms_per_line atoms
        if index % atoms_per_line == 0:
            indat_information_string += "\n           "

    # Remove the trailing comma and extra line
    indat_information_string = indat_information_string.rstrip(",\n") + "\n"

    return indat_information_string




def make_fmobnd(connectivity_list):
    # Create a formatted string for connectivity information
    connectivity_information_string = ""
    
    for i, connectivity_entry in enumerate(connectivity_list):
        # Extract the first two elements and the third for error checking
        first, second = connectivity_entry[:2]
        first += 1
        second += 1
        # Print the first and second elements
        connectivity_information_string += f"-{first} {second}"

        # Add a newline character if it's not the last element
        if i < len(connectivity_list) - 1:
            connectivity_information_string += "\n"

        # Check for errors
        #if error_check != 1:
        #    raise ValueError(f"Error in connectivity: Third element is not 1 (found {error_check})")

    return connectivity_information_string


def write_fmo_input_file(input_file_path, fmoxyz_content, data_content, icharg_content, indat_content, fmobnd_content, level, efmo):
    if(efmo):
        if(level > 2):
            print("EFMO can't do trimers, defaulting to dimers! \n")
        level = 2
        efmo_var = 1
    else:
        efmo_var = 0

    with open(input_file_path, "w") as output_file:
        output_file.write(" $contrl nprint=-5 runtyp=energy ispher=-1 local=ruednbrg $end\n")
        output_file.write(" $contrl maxit=200 $end\n")
        output_file.write(" $system mwords=2000 memddi=2000 $end\n")
        output_file.write(" $GDDI NGROUP=10 $END\n")
        output_file.write(" $SCF DIRSCF=.TRUE. DIIS=.t. $END\n")
        output_file.write(" $basis gbasis=n31 ngauss=6 ndfunc=1 $end\n")
        output_file.write(f" $FMO nbody={level} iefmo={efmo_var} rafo(1)=1,1,1\n")
        if(efmo):
            output_file.write(" MODEFM(1)=0, 16, 1, 1, 1\n")
        output_file.write(icharg_content + "  \n")
        output_file.write(indat_content + "  \n")
        output_file.write(" $end\n")
        output_file.write(" $fmoxyz\n")
        output_file.write(fmoxyz_content + "\n")
        output_file.write(" $end\n")
        output_file.write(" $FMOBND\n")
        output_file.write(fmobnd_content + "\n")
        output_file.write(" $end\n")
        output_file.write(" $DATA\n")
        output_file.write("Title\nC1\n")
        output_file.write(data_content + "\n")
        output_file.write(" $end\n")
