function xyz_to_geometry(xyzfile)
  input_file = []
  #println(joinpath(@__DIR__, xyzfile))
  open(joinpath(@__DIR__, xyzfile)) do file
    input_file = readlines(file)
  end


  natoms = parse(Int, input_file[1])
  geometry = input_file[3:natoms+2]

  coordinates = []
  symbols = []
  for geometry_line in geometry
    parts = split(strip(geometry_line))
    push!(symbols, parts[1])
    x = parse(Float64, parts[2])
    y = parse(Float64, parts[3])
    z = parse(Float64, parts[4])
    push!(coordinates, x)
    push!(coordinates, y)
    push!(coordinates, z)
  end

  charge = 0

  return coordinates, symbols, charge
end


function create_input_rhf(input, basis)
  #== write input json file ==#
  out = split(input, ".")
  jsonout = out[1]
  println("Transforming $jsonout to json using $basis")
  output = open("$jsonout.json", "w") do file
    write(file, "{\n")
    #== write in molecule information ==#
    write(file, "  \"topologies\": [{\n")
    #== write each fragment's molecule information ==#
    #write(file, "    \"", input, "\": {\n")

    coords, symbols, charge = xyz_to_geometry(input)
    #== write fragment geometry ==#
    write(file, "      \"geometry\" : [\n")
    for icoord = 1:length(coords)
      write(file, "        ")

      value = coords[icoord]
      if icoord == length(coords)
        write(file, "$value\n")
      elseif icoord % 3 == 0
        write(file, "$value,\n")
      else
        write(file, "$value,")
      end
    end
    write(file, "      ],\n")

    #== write fragment symbols ==#
    write(file, "      \"symbols\" : [\n")
    for iatom = 1:length(symbols)
      write(file, "        ")
      symbol = symbols[iatom]
      if iatom == length(symbols)
        write(file, "\"$symbol\"\n")
      elseif iatom % 5 == 0
        write(file, "\"$symbol\",\n")
      else
        write(file, "\"$symbol\",")
      end
    end
    write(file, "      ],\n")
    #== write fragment charge ==#
    write(file, "      \"fragment_formal_charges\": [$charge]\n")
    #write(file, "    },\n")
    write(file, "  }],\n")

    #== write calculation driver and model information ==#
    write(file, "  \"driver\": \"Energy\",\n")
    write(file, "  \"model\": {\n")
    write(file, "    \"method\": \"RestrictedHF\",\n")
    write(file, "    \"basis\": \"$basis\",\n")
    write(file, "    \"aux_basis\": \"$basis\"\n")
    write(file, "  },\n")

    #== write keywords ==#
    write(file, "  \"keywords\": {\n")
    write(file, "    \"scf\": {\n")
    write(file, "      \"max_iters\":50,\n")
    write(file, "      \"max_diis_history_length\":10,\n")
    write(file, "      \"convergence_threshold\":1E-6,\n")
    write(file, "      \"density_threshold\":1E-10,\n")
    write(file, "      \"convergence_metric\":\"Diis\"\n")
    write(file, "    }\n")
    write(file, "  },\n")
    write(file, "  \"system\": {\n")
    write(file, "  \"max_gpu_memory_mb\": 1000 \n")
    write(file, "  }\n")
    write(file, "}")
  end
end

input_directory = "/home/jorgegv/work/scripts_jorge/get_xyz_from_quick/xyz_files"
output_directory = "/home/jorgegv/work/scripts_jorge/make_jsons_from_xyz/json_inputs"

# Loop through XYZ files
for file in readdir(input_directory)
    if endswith(file, ".xyz")
        input_path = joinpath(input_directory, file)
        output_path = joinpath(output_directory, replace(file, ".xyz" => ".output"))

        # Call your function for each XYZ file
        create_input_rhf(input_path, "cc-pVDZ")

        # Optionally, you can save the output or perform additional actions
        # For example, save the output to a file
        # write(output_path, "output_data_here")
    end
end

