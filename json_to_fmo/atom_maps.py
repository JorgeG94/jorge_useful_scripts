def create_symbol_to_element_map():
    # Map atomic symbols to element numbers (as floats)
    symbol_to_element = {
        "H": 1.0, "He": 2.0, "Li": 3.0, "Be": 4.0, "B": 5.0, "C": 6.0,
        "N": 7.0, "O": 8.0, "F": 9.0, "Ne": 10.0, "Na": 11.0, "Mg": 12.0,
        "Al": 13.0, "Si": 14.0, "P": 15.0, "S": 16.0, "Cl": 17.0, "K": 19.0,
        "Ar": 18.0, "Ca": 20.0, "Sc": 21.0, "Ti": 22.0, "V": 23.0, "Cr": 24.0,
        "Mn": 25.0, "Fe": 26.0, "Ni": 28.0, "Co": 27.0, "Cu": 29.0, "Zn": 30.0,
        "Ga": 31.0, "Ge": 32.0, "As": 33.0, "Se": 34.0, "Br": 35.0, "Kr": 36.0,
        "Rb": 37.0, "Sr": 38.0, "Y": 39.0, "Zr": 40.0, "Nb": 41.0, "Mo": 42.0,
        "Tc": 43.0, "Ru": 44.0, "Rh": 45.0, "Pd": 46.0, "Ag": 47.0, "Cd": 48.0,
        "In": 49.0, "Sn": 50.0, "Sb": 51.0, "Te": 52.0, "I": 53.0, "Xe": 54.0,
        "Cs": 55.0, "Ba": 56.0, "La": 57.0, "Ce": 58.0, "Pr": 59.0, "Nd": 60.0,
        "Pm": 61.0, "Sm": 62.0, "Eu": 63.0, "Gd": 64.0, "Tb": 65.0, "Dy": 66.0,
        "Ho": 67.0, "Er": 68.0, "Tm": 69.0, "Yb": 70.0, "Lu": 71.0, "Hf": 72.0,
        "Ta": 73.0, "W": 74.0, "Re": 75.0, "Os": 76.0, "Ir": 77.0, "Pt": 78.0,
        "Au": 79.0, "Hg": 80.0, "Tl": 81.0, "Pb": 82.0, "Bi": 83.0, "Th": 90.0,
        "Pa": 91.0, "U": 92.0, "Np": 93.0, "Pu": 94.0, "Am": 95.0, "Cm": 96.0,
        "Bk": 97.0, "Cf": 98.0, "Es": 99.0, "Fm": 100.0, "Md": 101.0, "No": 102.0,
        "Lr": 103.0, "Rf": 104.0, "Db": 105.0, "Sg": 106.0, "Bh": 107.0, "Hs": 108.0,
        "Mt": 109.0, "Ds": 110.0, "Rg": 111.0, "Cn": 112.0, "Nh": 113.0, "Fl": 114.0,
        "Mc": 115.0, "Lv": 116.0, "Ts": 117.0, "Og": 118.0
    }

    return symbol_to_element

def map_symbols_to_elements(symbols, symbol_to_element):
    # Map atomic symbols to their element numbers
    element_numbers = [symbol_to_element[symbol] for symbol in symbols]

    return element_numbers