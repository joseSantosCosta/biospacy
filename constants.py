IUPAC = set(["A","T","C","G","R","Y","S","W","K","M","B","D","H","V","N",".","-"])

PURINE = set(["A","G","R","W","M","K","D","B","H","V","N"])

PYRIMIDINE = set(["C", "T", "U", "Y", "S", "W", "K", "M", "B", "D", "H", "N"])

STANDARD = set(["A","T","G","C"]) #if not in this, then it is ambiguous

INVALID = set(["[SPC]","[INV]"])

STOP_CODONS = set(["TAA","TAG","TGA"])


RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "NotI": "GCGGCCGC",
    "PstI": "CTGCAG",
    "SmaI": "CCCGGG",
    "XmaI": "CCCGGG",
    "KpnI": "GGTACC",
    "SacI": "GAGCTC",
    "SalI": "GTCGAC",
    "XhoI": "CTCGAG",
    "NheI": "GCTAGC",
    "SpeI": "ACTAGT",
    "AvrII": "CCTAGG",
    "ApaI": "GGGCCC",
    "BglII": "AGATCT",
    "ClaI": "ATCGAT",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "SphI": "GCATGC",
}

