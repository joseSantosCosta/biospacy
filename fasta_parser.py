def fasta_to_list(fasta):
    with open(fasta,'r') as f:
        lines = f.readlines()
        lines_striped = [line.strip() for line in lines]
        for line in lines_striped:
            if line == "":
                lines_striped.remove(line)
        
    return lines_striped

def remove_descriptions(lines:list):
    tokenized_fasta = []
    separated_description = [line.split() for line in lines]
    # A list of list, where each list is a line

    for line in separated_description:
        if line[0].startswith(">"):
            while len(line) > 1:
                line.pop()
    
    no_description = [line[0] for line in separated_description]
    return no_description


def extract_id_sequence(fasta):
    id_sequence:dict = {}
    to_extract = remove_descriptions(fasta_to_list(fasta)) #result of a list where it only has ID's and sequences
    id = ""
    sequence = ""
    for string in to_extract:
        if string.startswith(">"):
            id = string
            id_sequence[id] = ""
            sequence = ""
            
        else:
            id_sequence[id] += string

    for record in id_sequence:
        id_sequence[record] = [id_sequence[record]] 
    
    return id_sequence
        
        
        

    

