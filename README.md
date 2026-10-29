# Protein Sequence Explorer

A Python toolkit for reading protein FASTA files, validating sequences, and calculating biochemical features. The parser and features are written from scratch and checked against Biopython.

This is the first in a series of protein bioinformatics projects. The goal is a reliable foundation of hand-built sequence features to reuse later as inputs to machine learning models and as baselines against them, including protein language models (ESM), mutation-effect prediction, and protein design. 

## Status

| Milestone | Status |
|---|---|
| 1. FASTA parsing and validation | Done |
| 2. Sequence features | Done |
| 3. DataFrames, CSV export, `analyze` CLI | Next |
| 4. Visualization | Planned |
| 5. Biopython validation report | Planned |
| 6. Streamlit app | Planned |
| 7. Documentation and demo | Planned |

## Running it

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```bash
uv venv
uv pip install -e ".[dev]"

# Validate a FASTA file and print a summary table
uv run python -m protein_explorer summarize data/sample.fasta

# Run the tests
uv run pytest
```

## Data

`data/sample.fasta` contains 8 reviewed Swiss-Prot entries downloaded from [UniProtKB](https://www.uniprot.org) (CC BY 4.0): lysozyme C, insulin, bovine serum albumin, p53, calmodulin-1, GFP, β2-microglobulin, and hemoglobin β. Accessions and the download command are in [data/README.md](data/README.md).

These are full **precursor** sequences, including signal peptides and propeptides, not the mature proteins.

## Pipeline

1. **Parsing** (`fasta.py`). Reads records without modifying residues. Malformed files (sequence before a header, empty records, duplicate IDs) raise `FastaFormatError` with a line number.
2. **Validation** (`validation.py`). Flags invalid characters with their positions. Counts ambiguous (B/Z/J/X) and nonstandard (U/O) residues, accepts a terminal stop `*`, and rejects sequences that look like DNA or RNA (at least 20 residues long and at least 90% ACGTUN). Invalid sequences are reported rather than dropped.
3. **Features** (`features.py`), described below.

## Features

All features are **estimates calculated from sequence**, not measured values. Ambiguous and nonstandard residues count toward length. A terminal `*` does not.

| Feature | How it's calculated |
|---|---|
| Length | Number of residues |
| Residue counts and frequencies | Count / length for each of the 20 standard residues, always 20 keys in alphabetical order |
| Unique residue types | Number of standard residues present (0–20) |
| Hydrophobic fraction | (A+V+I+L+M+F+W) / length |
| Charged fraction | (D+E+K+R) / length. His is excluded because it is mostly neutral at pH 7.4 |
| Aromatic fraction | (F+W+Y) / length, the same definition as Biopython |
| Acidic and basic counts | D+E and K+R+H |
| GRAVY | Mean Kyte–Doolittle hydropathy over standard residues |
| Molecular weight | Sum of average residue masses + one water (Da). Assumes no PTMs and reduced Cys |
| Net charge at pH | Henderson–Hasselbalch over side chains and both termini, using EMBOSS pKa values |
| Isoelectric point | pH where net charge = 0, found by bisection on pH 0–14 |

## Findings 

- **Molecular weight:** we calculate 69.3 kDa for BSA, against the usual 66.4 kDa. The gap is the 24-residue signal peptide and propeptide in the precursor sequence.
- **Disulfides:** BSA's 17 disulfide bonds make the true mass about 34 Da lower than we calculate, because we assume every Cys is reduced. This is small compared with the precursor effect.
- **Isoelectric point:** our pI values differ from Biopython's by up to 0.58 units (HBB: 7.32 vs 6.74). This comes from using EMBOSS rather than Bjellqvist pKa values.
- **GRAVY and solubility:** insulin has a positive GRAVY (+0.19 for the precursor, +0.22 for the mature A+B chains), yet it is a soluble hormone. GRAVY is a weak predictor of solubility.
- **Sequence length:** insulin is listed as 110 residues because UniProt stores preproinsulin. The active hormone is 51 residues.

## Project layout

```
src/protein_explorer/
  fasta.py         FASTA parser
  validation.py    sequence validation
  features.py      feature calculations
  __main__.py      command-line interface
tests/             pytest suite and malformed-file fixtures
data/              sample FASTA and provenance
```
