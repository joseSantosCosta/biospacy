import spacy
from spacy.tokens import Token
import constants


#residue level
def is_purine(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return token.text in constants.PURINE

def is_pyrimidine(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return token.text in constants.PYRIMIDINE

def is_standard(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return token.text in constants.STANDARD

def is_ambiguous(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return not token.text in constants.STANDARD and not token.text in constants.INVALID 

def is_invalid(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return token.text in constants.INVALID


Token.set_extension("is_purine",getter=is_purine)
Token.set_extension("is_pyrimidine",getter=is_pyrimidine)
Token.set_extension("is_standard",getter=is_standard)
Token.set_extension("is_ambiguous",getter=is_ambiguous)
Token.set_extension("is_invalid",getter=is_invalid)

#window level
def contains_invalid(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return "[SPC]" in token.text or "[INV]" in token.text

def is_start_codon(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return token.text == "ATG"

def is_stop_codon(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return token.text in constants.STOP_CODONS

def is_homopolymer(token):
    if not isinstance(token,Token):
        raise TypeError("The object passed is not a token object")
    return all(n == token.text[0] for n in token.text)

Token.set_extension("contains_invalid",getter=contains_invalid)
Token.set_extension("is_start_codon",getter=is_start_codon)
Token.set_extension("is_stop_codon",getter=is_stop_codon)
Token.set_extension("is_homopolymer",getter=is_homopolymer)