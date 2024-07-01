import argparse
import json

from atom_maps import * 
from fmo_utils import *



def read_json_to_gamess(json_data, output_file, level, efmo):

    #symbols = json_data["topology"]["symbols"]
    #geometry = json_data["topology"]["geometry"]
    #connectivity = json_data["topology"]["connectivity"]
    #fragments = json_data["topology"]["fragments"]
    #fragment_charges = json_data["topology"]["fragment_charges"]
    symbols = json_data["symbols"]
    geometry = json_data["geometry"]
    connectivity = json_data["connectivity"]
    fragments = json_data["fragments"]
    fragment_charges = json_data["fragment_formal_charges"]


    symbol_to_element_map = create_symbol_to_element_map()
    element_numbers = map_symbols_to_elements(symbols, symbol_to_element_map)

    fmoxyz_content = make_fmoxyz(symbols, element_numbers, geometry)
    data_content = make_data(symbols, element_numbers)
    icharg_content = make_icharg(fragment_charges)
    indat_content = make_indat(fragments)
    fmobnd_content = make_fmobnd(connectivity)



    write_fmo_input_file(output_file, fmoxyz_content, data_content, icharg_content, indat_content, fmobnd_content, level, efmo )


def read_qdxf_to_gamess(json_data, output_file, level, efmo):

    symbols = json_data["symbols"]
    geometry = json_data["geometry"]
    connectivity = json_data["connectivity"]
    fragments = json_data["fragments"]
    fragment_charges = json_data["fragment_formal_charges"]


    symbol_to_element_map = create_symbol_to_element_map()
    element_numbers = map_symbols_to_elements(symbols, symbol_to_element_map)

    fmoxyz_content = make_fmoxyz(symbols, element_numbers, geometry)
    data_content = make_data(symbols, element_numbers)
    icharg_content = make_icharg(fragment_charges)
    indat_content = make_indat(fragments)
    fmobnd_content = make_fmobnd(connectivity)

    write_fmo_input_file(output_file, fmoxyz_content, data_content, icharg_content, indat_content, fmobnd_content, level, efmo )


def convert_to_gamess(json_data, output_file, level, efmo):
    # Your conversion logic from JSON to GAMESS format goes here
    # Replace the following print statement with your actual conversion code
    read_json_to_gamess(json_data, output_file, level, efmo)
    print("Converted from JSON to GAMESS format")

def convert_to_inp(json_data):
    # Your conversion logic from JSON to INP format goes here
    # Replace the following print statement with your actual conversion code
    
    print("Converted from JSON to INP format")

def convert_to_gamess_from_qdxf(qdxf_data, output_file, level, efmo):
    # Your conversion logic from qdxf to GAMESS format goes here
    # Replace the following print statement with your actual conversion code
    read_qdxf_to_gamess(qdxf_data, output_file, level, efmo)
    print("Converted from qdxf to GAMESS format")

def main():
    parser = argparse.ArgumentParser(description="Convert fragmented system files.")
    parser.add_argument("--from", dest="from_format", choices=["json", "qdxf"], required=True,
                        help="Input format (json or qdxf)")
    parser.add_argument("--to", dest="to_format", choices=["gamess", "inp"], required=True,
                        help="Output format (gamess or inp)")
    parser.add_argument("--input", dest="input_file", required=True,help="Input file" )
    parser.add_argument("--output", dest="output_file", required=True, help="Output input file")
    parser.add_argument("--frag_level", dest="level", type=int, default=2, help="Fragmentation level")
    parser.add_argument("--use_efmo", dest="efmo", action="store_true", help="Use efmo")


    args = parser.parse_args()

    # Read the input file
    with open(args.input_file, "r") as input_file:
        if args.from_format == "json":
            input_data = json.load(input_file)
        elif args.from_format == "qdxf":
            input_data = json.load(input_file)
        else:
            print("Invalid input format specified.")
            return

    # Perform the conversion based on the specified formats
    if args.to_format == "gamess":
        if args.from_format == "json":
            convert_to_gamess(input_data, args.output_file, args.level, args.efmo)
        elif args.from_format == "qdxf":
            convert_to_gamess_from_qdxf(input_data, args.output_file, args.level, args.efmo)
    elif args.to_format == "inp":
        convert_to_inp(input_data)
    else:
        print("Invalid output format specified.")

if __name__ == "__main__":
    main()
