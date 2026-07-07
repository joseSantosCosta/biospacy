from tokenizer import SequenceTokenizer
from fasta_parser import fasta_parser, FastaRecord
from matcher import MatchResult,  BioMatcher
import spacy
from spacy.tokens import Doc
from dataclasses import dataclass

@dataclass
class PipelineResult:
    record: FastaRecord
    doc: Doc
    match : list


nlp = spacy.blank("xx")

def buil_pipeline(token_mode,k=3,frame=0):
    nlp = spacy.blank("xx")
    tokenizer = SequenceTokenizer(vocab= nlp.vocab, token_mode= token_mode, k=k, frame=frame)
    nlp.tokenizer = tokenizer
    import attributes
    return nlp

def run_pipeline(fasta_path,nlp,bio_matcher = None):
    fasta_records = fasta_parser(fasta_path)

    for record in fasta_records:
        doc = nlp(record.sequence)
        doc._.sequence_id = record.seq_id
        if bio_matcher is not None:
            matches = bio_matcher.match(doc)
        else:
            matches = []
        
        yield PipelineResult(record=record,doc=doc,match=matches)






    

