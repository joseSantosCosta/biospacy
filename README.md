# BioSpaCy

A spaCy-based toolkit for biological sequence analysis. BioSpaCy brings the pattern matching capabilities of modern NLP pipelines to DNA sequence analysis, offering a more readable and maintainable alternative to regex for bioinformatics workflows.

This project started as part of an MSc in Bioinformatics and is being actively developed further.

---

## The problem with regex in bioinformatics

If you have spent any time doing sequence analysis, you have probably written something like this:

```python
import re
pattern = re.compile(r'ATG(?:[ATGC]{3})*?(?:TAA|TAG|TGA)')
matches = pattern.findall(sequence)
```

It works. But six months later, you open that script and spend ten minutes figuring out what it does. And that is a simple example — add ambiguity codes, variable length motifs, or physicochemical constraints and it gets ugly fast.

BioSpaCy takes a different approach. Instead of encoding biological rules as character patterns, it treats sequences the way NLP treats text — as a stream of tokens, each carrying biological metadata. Pattern matching then becomes a matter of describing what you are looking for in structured, readable terms:

```python
# Find all open reading frames
bio_matcher.add_orf_rule()
results = bio_matcher.match(doc)
```

The goal is not to replace established tools like BLAST or HMMER for their intended purposes. The goal is to make the kind of exploratory, rule-based sequence analysis that researchers do every day more interpretable and easier to maintain.

---

## Features

- **Three tokenization modes** — single residue, codon (with reading frame support), and non-overlapping k-mers
- **Full IUPAC DNA alphabet support** — ambiguity codes are first-class citizens, not an afterthought
- **Biological token attributes** — each token carries metadata like `is_purine`, `is_start_codon`, `is_stop_codon`, `is_homopolymer`, and more
- **Rule-based motif detection** — write readable pattern rules using token attributes rather than character strings
- **Built-in biological rules** — ORF detection, restriction site matching with IUPAC expansion, homopolymer detection
- **FASTA support** — parse single and multi-sequence FASTA files directly into the pipeline
- **spaCy native** — built on spaCy 3.x, compatible with the broader spaCy ecosystem

---

## Installation

Clone the repository and install the dependencies:

```bash
git clone https://github.com/yourusername/biospacy.git
cd biospacy
pip install -r requirements.txt
```

Dependencies:
- spaCy >= 3.0
- Biopython
- Python >= 3.9

---

## Quick start

### Basic tokenization

```python
import spacy
from tokenizer import SequenceTokenizer

# Build a pipeline in codon mode
nlp = spacy.blank("xx")
nlp.tokenizer = SequenceTokenizer(nlp.vocab, token_mode="codon", frame=0)

# Tokenize a sequence
doc = nlp("ATGCATGCATGA")

for token in doc:
    print(token.text, token._.is_start_codon, token._.is_stop_codon)
```

Output:
```
ATG  True   False
CAT  False  False
GCA  False  False
TGA  False  True
```

### ORF detection

```python
from pipeline import build_pipeline, run_pipeline
from matcher import BioMatcher

# Build a codon mode pipeline
nlp = build_pipeline(token_mode="codon", frame=0)

# Set up the matcher
bio_matcher = BioMatcher(nlp)
bio_matcher.add_orf_rule()

# Run on a FASTA file
for result in run_pipeline("sequences.fasta", nlp, bio_matcher):
    print(f"Sequence: {result.record.seq_id}")
    for match in result.matches:
        print(f"  ORF found at positions {match.start_position}-{match.end_position}")
        print(f"  Sequence: {match.matched_sequence}")
        print(f"  Length: {match.length} nucleotides")
```

### Restriction site detection

```python
nlp = build_pipeline(token_mode="residue")
bio_matcher = BioMatcher(nlp)

# Add a restriction site rule — handles IUPAC ambiguity codes automatically
bio_matcher.add_restriction_site_rule("EcoRI")

for result in run_pipeline("sequences.fasta", nlp, bio_matcher):
    for match in result.matches:
        print(f"EcoRI site at position {match.start_position}")
```

### Custom pattern rules

You can define your own rules using token attributes:

```python
# Find stretches of 5 or more purines
pattern = [
    {"_": {"is_purine": True}, "OP": "{5,}"}
]
bio_matcher.add_rule("PURINE_STRETCH", pattern)
```

This is where BioSpaCy really shines compared to regex — the rule describes what you are looking for biologically, not how to encode it as characters.

---

## Tokenization modes

### Residue mode
Each nucleotide becomes one token. Best for single-nucleotide level analysis — restriction sites, CpG dinucleotides, splice sites.

```python
nlp = build_pipeline(token_mode="residue")
doc = nlp("ATGCAT")
# tokens: ["A", "T", "G", "C", "A", "T"]
```

### Codon mode
Non-overlapping triplets. Supports reading frames 0, 1, and 2. Best for coding sequence analysis — ORF detection, codon usage, start/stop codon identification.

```python
nlp = build_pipeline(token_mode="codon", frame=0)
doc = nlp("ATGCATGCA")
# tokens: ["ATG", "CAT", "GCA"]
```

### K-mer mode
Non-overlapping windows of size k. Best for repeat detection, sequence composition analysis, and pattern matching at a defined resolution.

```python
nlp = build_pipeline(token_mode="kmer", k=4)
doc = nlp("ATGCATGC")
# tokens: ["ATGC", "ATGC"]
```

---

## Token attributes

Every token in a BioSpaCy Doc carries biological metadata accessible via `token._.<attribute>`:

| Attribute | Modes | Description |
|---|---|---|
| `is_purine` | residue | True if the nucleotide is a purine or ambiguous purine |
| `is_pyrimidine` | residue | True if the nucleotide is a pyrimidine or ambiguous pyrimidine |
| `is_standard` | all | True if the nucleotide is A, T, G, or C |
| `is_ambiguous` | all | True if the nucleotide is a valid IUPAC ambiguity code |
| `is_invalid` | all | True if the token is [INV] or [SPC] |
| `is_start_codon` | codon | True if the token is ATG |
| `is_stop_codon` | codon | True if the token is TAA, TAG, or TGA |
| `contains_invalid` | codon, kmer | True if the window contains an invalid token |
| `is_homopolymer` | all | True if all characters in the token are identical |

---

## Handling invalid characters

Any character outside the IUPAC DNA alphabet is replaced with a special `[INV]` token. Whitespace becomes `[SPC]`. Both are treated as single positional units — they occupy one position in the sequence just like any valid nucleotide, preserving positional integrity throughout the pipeline.

```python
doc = nlp("ATG?CA")
# tokens: ["A", "T", "G", "[INV]", "C", "A"]
```

The positions of all replacements are logged at INFO level and accessible via `doc._.replacement_log`.

---

## Project structure

```
bio_spacy/
│
├── constants.py        — IUPAC alphabet, special tokens, restriction enzymes
├── tokenizer.py        — SequenceTokenizer class and factory registration
├── attributes.py       — biological token and Doc extensions
├── fasta_parser.py     — FastaRecord class and parse_fasta generator
├── matcher.py          — BioMatcher class and MatchResult dataclass
├── pipeline.py         — build_pipeline and run_pipeline orchestration
├── logging_config.py   — logging setup utility
└── tests/
      └── test_tokenizer.py
```

---

## Roadmap

This is an active MSc project. Planned additions include:

- Splice site detection (donor and acceptor sites)
- Kozak sequence detection
- CpG dinucleotide analysis
- Protein sequence support (amino acid tokenization with physicochemical attributes)
- k-mer frequency analysis
- Prodigy-style annotation interface for sequence feature labelling
- Benchmarking against regex for interpretability and maintainability

---

## Contributing

This project is in active development. If you are a bioinformatician or Python developer and find it useful or have ideas for biological rules worth implementing, feel free to open an issue or get in touch.

---

## License

MIT

---

## Acknowledgements

Built with [spaCy](https://spacy.io) by Explosion AI. Restriction enzyme data from [REBASE](http://rebase.neb.com) via [Biopython](https://biopython.org).