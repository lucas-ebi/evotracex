from Bio.Align import substitution_matrices

ALPHABET: list[str] = [
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y",
]

STEREOCHEMISTRY: dict[str, list[str]] = {
    "Aliphatic":              ["G", "A", "V", "L", "I"],
    "Amide":                  ["N", "Q"],
    "Aromatic":               ["F", "Y", "W"],
    "Basic":                  ["H", "K", "R"],
    "Big":                    ["M", "I", "L", "K", "R"],
    "Hydrophilic":            ["R", "K", "N", "Q", "P", "D"],
    "Median":                 ["E", "V", "Q", "H"],
    "Negatively charged":     ["D", "E"],
    "Non-polar":              ["F", "G", "V", "L", "A", "I", "P", "M", "W"],
    "Polar":                  ["Y", "S", "N", "T", "Q", "C"],
    "Positively charged":     ["K", "R"],
    "Similar (Asn or Asp)":   ["N", "D"],
    "Similar (Gln or Glu)":   ["Q", "E"],
    "Small":                  ["C", "D", "P", "N", "T"],
    "Tiny":                   ["G", "A", "S"],
    "Very hydrophobic":       ["L", "I", "F", "W", "V", "M"],
    "With hydroxyl":          ["S", "T", "Y"],
    "With sulfur":            ["C", "M"],
}

BLOSUM62 = substitution_matrices.load("BLOSUM62")
BLOSUM62_THRESHOLD: float = float(
    sum(BLOSUM62[a, a] for a in ALPHABET) / len(ALPHABET)
)
