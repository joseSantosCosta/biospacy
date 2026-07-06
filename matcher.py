from dataclasses import dataclass
from spacy.tokens import Span, Doc
from spacy.matcher import Matcher
from constants import RESTRICTION_ENZYMES

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
            matched_sequence = doc[start:end]
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
        if enzyme not in RESTRICTION_ENZYMES:
            raise ValueError("The enzyme selected is unknown for the system")
        
        pattern = [{"TEXT":RESTRICTION_ENZYMES[enzyme]}]
        rule_name = f"{enzyme}_pattern"
        self.matcher.add(rule_name,[pattern])
    
    def add_homopolymer_rule(self,min_length):
        pattern = [
            {"_":{"is_homopolymer":True},"OP":"{min_length}"}
        ]

        rule_name = "HOMOPOLYMER_PATTERN"
        self.matcher.add(rule_name,[pattern])


        


