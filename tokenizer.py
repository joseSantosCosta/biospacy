import spacy
from spacy import util
from spacy.language import Language
from spacy.tokens import Doc, DocBin
from spacy.vocab import Vocab
from constants import IUPAC
from logger import setup_logger
import logging
import fasta_parser
from pathlib import Path

loggerObj = logging.getLogger(__name__)
setup_logger()

Doc.set_extension("sequence_id",default=None)
Doc.set_extension("record_log",default=None)
Doc.set_extension("token_mode",default=None)
Doc.set_extension("k",default=None)
Doc.set_extension("frame",default=None)



class SequenceTokenizer():
    def __init__(self,vocab, token_mode="residue", k=3, frame=0):
        if not isinstance(vocab,Vocab):
            raise TypeError("The vocab has to be a Vocab object from the Spacy library")
        if token_mode not in ["residue","codon","kmer"]:
            raise ValueError("The token_mode has to be residue, codon or kmer")
        
        if (token_mode == "residue" or token_mode == "kmer") and frame != 0:
            raise ValueError("The frame was adjusted, but the token mode selected doesn't account for frame, its better if you leave it as default")
        
        if (token_mode == "residue" or token_mode == "codon") and k != 3:
            raise ValueError("The kmer should only be adjusted when using kmer token mode, leave that value on default")
        
        if token_mode == "kmer" and k < 2:
            raise ValueError("K value cannot be lower than 2 since k=1 is a residue split and having it at 0 or below is not possible")
        
        if token_mode == "codon" and frame < 0 :
            raise ValueError("Frame value cannot be lower than 0")

        self.vocab = vocab
        self.token_mode = token_mode #residue, codon, kmer
        self.k = k
        self.frame = frame #how many nucleotides it skips to start the codon, if frame = 0, then it doesnt skip any

    
    def _replace_invalid_tokens(self,sequence):
        replaced_nucleotides = [] #list of dictionaries where every dict has a position in the sequence as key and a tuple of the nucleotide and what it was replaced for
        if not isinstance(sequence,str):
            raise TypeError ("Please make sure that the sequence passed is a string")
        
        sequence = list(sequence)
        for i, nucleotide in enumerate(sequence):
            if nucleotide not in IUPAC:
                if nucleotide == " ":
                    replaced_nucleotides.append({i : [nucleotide, "[SPC]"]})
                    sequence[i] = "[SPC]"
                else:
                    replaced_nucleotides.append({i: [nucleotide, "[INV]"]})
                    sequence[i] = "[INV]"
        
        for replaced in replaced_nucleotides:
            replaced_keys = replaced.keys()
            replaced_values = list(replaced.values())
            loggerObj.info(f"The nucleotide in position {next(iter(replaced_keys))} was an {replaced_values[0][0]} and it was replace by {replaced_values[0][1]}")
        return sequence
    
    def _split_into_windows(self, token_list, window_size, offset):
        if window_size < 2:
            raise ValueError("The k has to be greater or equal 2")
        if not isinstance(window_size,int):
            raise TypeError("The k selected has to be an integer")
        if not isinstance(offset,int):
            raise TypeError("The frame selected has to be an integer")
        if offset < 0:
            raise ValueError("Frame cannot be lower than 0")
        if offset >= len(token_list):
            raise ValueError("Frame cannot be higher than the length of the sequence")
        if not token_list:
            raise ValueError("The sequence given cannot be null")
        if not isinstance(token_list,list):
            raise TypeError("This function has to receive the list version of the sequence, not the string")
        if len(token_list) - offset < window_size:
            raise ValueError

        stop = offset + ((len(token_list)-offset) // window_size) * window_size   
        windows = range(offset,stop,window_size)
        split_window_sequence = []
        skipped_positions = offset
        truncated_positions = (len(token_list) - offset) % window_size

        if skipped_positions != 0:
            loggerObj.info(f"With the frame set to {offset}, you will skip {skipped_positions} nucleotides")

        if truncated_positions != 0:
            loggerObj.warning(f"Your choice on the kmer size {truncated_positions} nucleotides truncated")

        for i in windows:
            split_window_sequence.append("".join(token_list[i:i+window_size]))
        
        return split_window_sequence
        

    #runs this function when token_mode = residue
    def _split_residue(self,sequence):
        sequence = self._replace_invalid_tokens(sequence)
        return sequence
    
    def _split_codon(self,sequence):
        if self.frame > len(sequence):
            raise ValueError ("The frame defined is higher than the lenght of the sequence")
        
        token_list = self._replace_invalid_tokens(sequence)
        return self._split_into_windows(token_list,3,self.frame)
    
    def _split_kmer(self,sequence):
        if len(sequence) < self.k:
            raise ValueError("The k defined is bigger than the length of the sequence")
        
        token_list = self._replace_invalid_tokens(sequence)
        return self._split_into_windows(token_list,self.k,0)
    
    def __call__(self, sequence):
        if sequence == "":
            raise ValueError("The sequence is empty, cannot tokenize a empty sequence")
        if not isinstance(sequence,str):
            raise TypeError("The offered sequenced is not a string")
        
        if self.token_mode=="residue":
            tokenized_sequence = self._split_residue(sequence=sequence)
            spaces = [False for token in tokenized_sequence]
            doc =  Doc(self.vocab, words=tokenized_sequence,spaces=spaces)
            loggerObj.debug(f"The sequence was tokenized into {len(tokenized_sequence)} residues")
            

        elif self.token_mode == "codon":
            tokenized_sequence = self._split_codon(sequence=sequence)
            spaces = [False for token in tokenized_sequence]
            doc = Doc(self.vocab,words=tokenized_sequence,spaces=spaces)
            loggerObj.debug(f"The sequence was tokenized into {len(tokenized_sequence)} codons")
            

        else:
            tokenized_sequence = self._split_kmer(sequence=sequence)
            spaces = [False for token in tokenized_sequence]
            doc = Doc(self.vocab, words=tokenized_sequence,spaces = spaces)
            loggerObj.debug(f"The sequence was tokenized into {len(tokenized_sequence)} kmers")

        doc._.token_mode = self.token_mode
        doc._.k = self.k
        doc._.frame = self.frame
        
        return doc  
    
    def to_bytes(self):
        serializers = {
            "token_mode" : lambda b: str(self.token_mode).encode("utf8"),
            "k": lambda b: str(self.k).encode("utf8"),
            "frame" : lambda b : str(self.frame).encode("utf8")
        }
        
        return util.to_bytes(serializers)
 
    def from_bytes(self, bytes_data):
        deserializers = {
            "token_mode" : lambda b: setattr(self, "token_mode",b.decode("utf8") ),
            "k": lambda b: setattr(self,"k",int(b.decode("utf8"))),
            "frame" : lambda b : setattr(self,"frame",int(b.decode("utf8")))
        }

        util.from_bytes(bytes_data,deserializers)
        return self
    
    def to_disk(self,path):
        if path == None:
            raise ValueError("The path provided is none")
        path = Path(path)
        if not path.exists():
            raise ValueError("The path provided does not exist")
        with open(path,"wb") as f:
            f.write(self.to_bytes())
    
    def from_disk(self,path):
        if path == None:
            raise ValueError("The path provided is none")
        path = Path(path)
        if not path.exists():
            raise ValueError("The path provided does not exist")
        with open(path,"rb") as f:
            bytes_of_data = f.read()
            self.from_bytes(bytes_of_data)
        return self
    





        

        


                    
                
            

