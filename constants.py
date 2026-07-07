import json
IUPAC = set(["A","T","C","G","R","Y","S","W","K","M","B","D","H","V","N",".","-"])

PURINE = set(["A","G","R","W","M","K","D","B","H","V","N"])

PYRIMIDINE = set(["C", "T", "U", "Y", "S", "W", "K", "M", "B", "D", "H", "N"])

STANDARD = set(["A","T","G","C"]) #if not in this, then it is ambiguous

INVALID = set(["[SPC]","[INV]"])

STOP_CODONS = set(["TAA","TAG","TGA"])


with open("restriction_enzyme_json","r") as f:
    RESTRICTION_ENZYMES = json.load(f) 

ENZYME_LOOKUP = {name.lower():name for name in RESTRICTION_ENZYMES}

IUPAC_EXPANSION = {
    "A": {"A"},
    "T": {"T"},
    "G": {"G"},
    "C": {"C"},
    "R": {"A", "G"},           
    "Y": {"C", "T"},           
    "M": {"A", "C"},           
    "K": {"G", "T"},           
    "W": {"A", "T"},          
    "S": {"G", "C"},           
    "B": {"C", "G", "T"},      
    "D": {"A", "G", "T"},      
    "H": {"A", "C", "T"},      
    "V": {"A", "C", "G"},      
    "N": {"A", "C", "G", "T"}
}
