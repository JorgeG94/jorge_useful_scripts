#ifndef MOLECULE_H
#define MOLECULE_H

#include <vector>
#include <string>
#include <utility>
#include <cmath>
#include <map>

// Covalent radii in angstroms
static const std::map<std::string, double> covalentRadii = {
    {"H", 0.31}, {"He", 0.28}, {"Li", 1.28}, {"Be", 0.96}, {"B", 0.84}, 
    {"C", 0.76}, {"N", 0.71}, {"O", 0.66}, {"F", 0.57}, {"Ne", 0.58}, 
    {"Na", 1.66}, {"Mg", 1.41}, {"Al", 1.21}, {"Si", 1.11}, {"P", 1.07}, 
    {"S", 1.05}, {"Cl", 1.02}, {"Ar", 1.06}, {"K", 2.03}, {"Ca", 1.76}, 
    {"Sc", 1.70}, {"Ti", 1.60}, {"V", 1.53}, {"Cr", 1.39}, {"Mn", 1.39}, 
    {"Fe", 1.32}, {"Co", 1.26}, {"Ni", 1.24}, {"Cu", 1.32}, {"Zn", 1.22}, 
    {"Ga", 1.22}, {"Ge", 1.20}, {"As", 1.19}, {"Se", 1.20}, {"Br", 1.20}
};

class Molecule {
public:
    struct Atom {
        std::string element;
        std::vector<double> position; // Assuming a 3D position: x, y, z
    };

    struct Bond {
        int atomIndex1;
        int atomIndex2;
        int bondOrder;
    };

    static double calculateDistance(const Atom& atom1, const Atom& atom2){
    double dx = atom1.position[0] - atom2.position[0];
    double dy = atom1.position[1] - atom2.position[1];
    double dz = atom1.position[2] - atom2.position[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
    }
    static int bondOrderEstimation(const Atom& atom1, const Atom& atom2, double distance){
    double toleranceFactor = 0.4;
    auto findRadii = [](const std::string& element) { return covalentRadii.find(element)->second; };
    
    // Immediate return for hydrogen involvement
    if (atom1.element == "H" || atom2.element == "H") {
        double sumRadii = findRadii(atom1.element) + findRadii(atom2.element);
        return distance <= sumRadii + toleranceFactor ? 1 : 0;
    }

    double singleBondLength = findRadii(atom1.element) + findRadii(atom2.element);
    //if (distance <= singleBondLength * 0.7) return 3; // Triple bond
    //else if (distance <= singleBondLength * 0.9) return 2; // Double bond
    if (distance <= singleBondLength + toleranceFactor) return 1; // Single bond
    else return 0; // No bond
    }
    static std::vector<Bond> calculateBonds(const std::vector<Atom>& atoms){
    std::vector<Bond> bonds;
    for (size_t i = 0; i < atoms.size(); ++i) {
        for (size_t j = i + 1; j < atoms.size(); ++j) {
            double distance = calculateDistance(atoms[i], atoms[j]);
            int bondOrder = bondOrderEstimation(atoms[i], atoms[j], distance);
            if (bondOrder > 0) {
                bonds.push_back({static_cast<int>(i), static_cast<int>(j), bondOrder});
            }
        }
    }
    return bonds;
    }
};

#endif // MOLECULE_H

