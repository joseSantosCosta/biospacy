from dataclasses import dataclass
from spacy.tokens import Span, Doc
from spacy.matcher import Matcher
from constants import RESTRICTION_ENZYMES, ENZYME_LOOKUP, IUPAC_EXPANSION, STANDARD,IUPAC
import logging
import logger

logger.setup_logger()
loggerObj = logging.getLogger(__name__)

@dataclass
class MatchResult:
    rule_name:str
    start_position:int
    end_position:int
    matched_sequence:str
    length:int
    span: Span


class BioMatcher:
    def __init__(self,nlp):
        self.nlp = nlp
        self.matcher = Matcher(nlp.vocab)
    
    def add_rule(self,name,pattern):
        if not isinstance(name,str):
            raise TypeError("The name for the pattern has to be a string")
        if not name:
            raise ValueError("The name cannot be an empty string")
        if not isinstance(pattern,list):
            raise TypeError("The pattern has to be a list of dictionaries")
        if not all(isinstance(p,dict) for p in pattern):
            raise TypeError("The pattern has to be a list of dictionaries, at least one element in the provided list is not a dictionary")
        self.matcher.add(name,[pattern])

    def match(self,doc):
        if not isinstance(doc,Doc):
            raise TypeError("The object passed is not a Doc object")
        matches = self.matcher(doc)
        for match_id, start, end in matches:
            rule_name = self.nlp.vocab.strings[match_id]
            start_position = start
            end_position = end
            matched_sequence = doc[start:end].text
            length = len(matched_sequence)
            span = Span(doc,start,end)

            yield MatchResult(
                rule_name=rule_name,
                start_position=start_position,
                end_position=end_position,
                matched_sequence=matched_sequence,
                length=length,
                span = span
            )
    
    def add_orf_rule(self):
        pattern = [
            {"_":{"is_start_codon":True}},
            {"_":{"is_stop_codon":False}, "OP":"*"},
            {"_":{"is_stop_codon":True}}
        ]
        rule_name = "ORF_PATTERN"
        self.matcher.add(rule_name,[pattern])
    
        
    def add_restriction_site_rule(self,enzyme):
        enzyme_lower = enzyme.lower()
        if enzyme_lower not in ENZYME_LOOKUP:
            raise NameError(f"The enzyme name {enzyme} was not found in the enzyme database")
        
        else:
            enzyme = ENZYME_LOOKUP[enzyme_lower]
            enzyme_seq = RESTRICTION_ENZYMES[enzyme]
            pattern = _enzyme_sequence_expansion(enzyme_seq=enzyme_seq)
            rule_name = f"{enzyme}_pattern"
            self.matcher.add(rule_name,[pattern])
    

    #do I want to match for sequences like AAAANAAA
    def add_homopolymer_rule(self,min_length,standard_only = True):
        if min_length <= 1:
            raise ValueError("The minimum length has to be greater than 1")
        list_of_patterns = []
        if standard_only:
            for n in STANDARD:
                 list_of_patterns.append([{"TEXT":n}] * min_length + [{"TEXT":n,"OP":"*"}])
        else:
            for n in IUPAC:
                list_of_patterns.append([{"TEXT":n}] * min_length + [{"TEXT":n,"OP":"*"}])
           
        rule_name = "HOMOPOLYMER_PATTERN"
        self.matcher.add(rule_name,list_of_patterns)
        logging.debug(f"A pattern matching rule named {rule_name} was added together with {len(list_of_patterns)} possible patterns")

#helper function to ambiguous nucleotides in the enzyme recognition site sequence
def _enzyme_sequence_expansion(enzyme_seq,token_mode = "residue"):
        pattern = []
        for n in enzyme_seq:
            if n not in STANDARD:
                pattern.append({"TEXT":{"IN":IUPAC_EXPANSION[n]}})
            else:
                pattern.append({"TEXT":n})
        return pattern

        


