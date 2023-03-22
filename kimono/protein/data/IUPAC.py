"""IUPAC data for protein sequences."""

# From biopython 
# https://github.com/biopython/biopython 
# Bio.Data.IUPACData


protein_letters_1to3 = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}

protein_letters_1to3_extended = {
    **protein_letters_1to3,
    **{"B": "Asx", "X": "Xaa", "Z": "Glx", "J": "Xle", "U": "Sec", "O": "Pyl"},
}

protein_letters_3to1 = {value: key for key, value in protein_letters_1to3.items()}
protein_letters_3to1_extended = {
    value: key for key, value in protein_letters_1to3_extended.items()
}