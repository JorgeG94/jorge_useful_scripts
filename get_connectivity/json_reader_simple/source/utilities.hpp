#pragma once 
#include <iostream>
#include <fstream>
#include <nlohmann/json.hpp>
#include <filesystem> // Requires C++17
#include <molecule.hpp>
using json = nlohmann::json;
std::vector<Molecule::Atom> createAtomsFromJson(const std::string &filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "Could not open the file " << filePath << std::endl;
        throw std::runtime_error("File opening failed");
    }

    json j;
    file >> j;
    file.close();

    std::vector<Molecule::Atom> atoms;

    // Assuming each topology is a separate molecule
    for (const auto &topology : j["topologies"])
    {
        auto symbols = topology["symbols"];
        auto geometry = topology["geometry"];

        for (size_t i = 0; i < symbols.size(); ++i)
        {
            Molecule::Atom atom;
            atom.element = symbols[i];
            atom.position = {geometry[i * 3], geometry[i * 3 + 1], geometry[i * 3 + 2]};
            atoms.push_back(atom);
        }
    }

    return atoms;
}
std::vector<Molecule::Atom> createAtomsFromXyz(const std::string &filePath)
{
    std::ifstream file(filePath);
    if (!file.is_open())
    {
        std::cerr << "Could not open the file " << filePath << std::endl;
        throw std::runtime_error("File opening failed");
    }

    std::vector<Molecule::Atom> atoms;
    std::string line;
    std::getline(file, line); // Skip the first line containing the atom count
    std::getline(file, line); // Skip the comment line

    std::string element;
    double x, y, z;
    while (file >> element >> x >> y >> z)
    {
        atoms.push_back(Molecule::Atom{element, {x, y, z}});
    }

    file.close();
    return atoms;
}

void saveMoleculesToJson(const std::string &inputFilePath, const std::vector<Molecule::Bond> &bonds, const std::string &outputPath)
{
    // Read the original JSON file to preserve the input structure
    std::ifstream inputFile(inputFilePath);
    if (!inputFile.is_open())
    {
        std::cerr << "Could not open the input file " << inputFilePath << std::endl;
        throw std::runtime_error("File opening failed");
    }

    json originalJson;
    inputFile >> originalJson;
    inputFile.close();

    // Assuming each topology is a separate molecule, update it with connectivity
    for (auto &topology : originalJson["topologies"])
    {
        std::vector<json> connectivityJson;
        for (const auto &bond : bonds)
        {
            connectivityJson.push_back({bond.atomIndex1, bond.atomIndex2, bond.bondOrder});
        }
        topology["connectivity"] = connectivityJson;
    }

    // Write the updated JSON with connectivity back to a new file
    std::ofstream outputFile(outputPath);
    if (!outputFile.is_open())
    {
        std::cerr << "Failed to open " << outputPath << " for writing." << std::endl;
        return;
    }
    outputFile << originalJson.dump(4); // Pretty printing with an indent of 4 spaces
    outputFile.close();
}

void saveMoleculeToJsonFromXyz(const std::vector<Molecule::Atom> &atoms, const std::vector<Molecule::Bond> &bonds, const std::string &outputPath)
{
    json topology;
    topology["symbols"] = json::array();
    topology["geometry"] = json::array();

    topology["connectivity"] = json::array();

    for (const auto &atom : atoms)
    {
        topology["symbols"].push_back(atom.element);
        // Append each position component to the geometry array
        for (const auto &coord : atom.position)
        {
            topology["geometry"].push_back(coord);
        }
    }

    for (const auto &bond : bonds)
    {
        topology["connectivity"].push_back({bond.atomIndex1, bond.atomIndex2, bond.bondOrder});
    }

    json outputJson = {{"topologies", {topology}}};

    std::ofstream outputFile(outputPath);
    if (!outputFile.is_open())
    {
        std::cerr << "Failed to open " << outputPath << " for writing." << std::endl;
        return;
    }
    outputFile << outputJson.dump(4);
    outputFile.close();
}
