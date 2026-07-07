import json
from Bio.Restriction import AllEnzymes

db = {enzyme.__name__: enzyme.site for enzyme in AllEnzymes}

with open("restriction_enzyme_json","w") as f:
    json.dump(db,f,indent=4)



