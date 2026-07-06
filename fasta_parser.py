from dataclasses import dataclass
from pathlib import Path
import os
import logger
import logging

loggerObj = logging.getLogger(__name__)
logger.setup_logger()

@dataclass
class FastaRecord:
    seq_id:str
    description:str
    sequence:str

    def __str__(self):
        gc_content = round(100 * (self.sequence.count("G") + self.sequence.count("C")) / len(self.sequence),2) if len(self.sequence) > 0 else 0
        return f"Sequence id: {self.seq_id} \n Description: {self.description} \n Length of the sequence: {len(self.sequence)} \n GC content(%): {gc_content}"

def fasta_parser(fasta):
    fasta = Path(fasta)
    if not fasta.exists():
        raise FileNotFoundError("The fasta file passed was not found")
    
    if os.path.getsize(fasta) == 0:
        raise ValueError("The file selected is empty")
    
    current_id = None
    current_description = ""
    current_seq = []
    with open(fasta,"r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    logging.info(f"Creating a FastaRecord object instance for the sequence with id {current_id}")
                    yield FastaRecord(seq_id=current_id,description=current_description,seq="".join(current_seq))
                
            
                parts = line[1:].split(maxsplit=1)
                current_id = parts[0]
                current_description = parts[1] if len(parts) > 1 else ""
                current_seq = []
            else:
                current_seq.append(line)
        
        sequence = "".join(current_seq).upper()

        if not sequence:
            logging.warning(f"The sequence with id {current_id} is empty")

        yield FastaRecord(
            seq_id=current_id,
            description=current_description,
            sequence=sequence,
        )


            
    

