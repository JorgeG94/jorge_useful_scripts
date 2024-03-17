#include <iostream>
#include <fstream>
#include <filesystem> // Requires C++17
#include <molecule.hpp>
#include "utilities.hpp"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_json_file> [output_json_file]" << std::endl;
        return 1;
    }

    try
    {
        std::string inputFilePath = argv[1];
        std::string outputFilePath;

        // Generate the output file path as previously described
        if (argc == 3)
        {
            outputFilePath = argv[2];
        }
        else
        {
            std::filesystem::path inputPath(inputFilePath);
            std::string baseName = inputPath.stem().string();
            std::filesystem::path outputPath = "json_logs/" + baseName + ".json";
            outputFilePath = outputPath.string();
            std::filesystem::create_directory("json_logs");
        }

        std::vector<Molecule::Atom> atoms;
        std::filesystem::path inputPath(inputFilePath);
        bool is_input_xyz = false;
        if (inputPath.extension() == ".json")
        {
            atoms = createAtomsFromJson(inputFilePath);
        }
        else if (inputPath.extension() == ".xyz")
        {
            is_input_xyz = true;
            atoms = createAtomsFromXyz(inputFilePath);
        }
        else
        {
            std::cerr << "Unsupported file format. Please use .json or .xyz files." << std::endl;
            return 1;
        }

        auto bonds = Molecule::calculateBonds(atoms);
        // Adjust saveMoleculesToJson if necessary to handle creating JSON from XYZ data
        if (is_input_xyz)
        {
            saveMoleculeToJsonFromXyz(atoms, bonds, outputFilePath);
        }
        else
        {
            saveMoleculesToJson(inputFilePath, bonds, outputFilePath);
        }

        std::cout << "Processed and saved to " << outputFilePath << std::endl;
        return 0;

        std::cout << "Processed and saved to " << outputFilePath << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
