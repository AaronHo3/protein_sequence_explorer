# Protein Sequence Explorer

A Python toolkit for reading protein FASTA files, validating sequences, and calculating biochemical features. The parser and features are written from scratch and checked against Biopython.

This is the first in a series of protein bioinformatics projects. The goal is a reliable foundation of hand-built sequence features to reuse later as inputs to machine learning models and as baselines against them, including protein language models (ESM), mutation-effect prediction, and protein design. 

## Status

| Milestone | Status |
|---|---|
| 1. FASTA parsing and validation | Done |
| 2. Sequence features | Done |
| 3. DataFrames, CSV export, `analyze` CLI | Done |
| 4. Visualization | Done |
| 5. Biopython validation report | Done |
| 6. Streamlit app | Next |
| 7. Documentation and demo | Planned |

## Running it

Requires Python 3.11+ and [uv](https://docs.astral.sh/uv/).

```bash
uv venv
uv pip install -e ".[dev]"

# Validate a FASTA file and print a summary table
uv run python -m protein_explorer summarize data/sample.fasta

# Compute all features and write CSV tables to results/
uv run python -m protein_explorer analyze data/sample.fasta --output results/

# Same, with net charge calculated at pH 7.4
uv run python -m protein_explorer analyze data/sample.fasta --output results/ --ph 7.4

# Generate all figures in results/figures/ (PNG by default; also svg or pdf)
uv run python -m protein_explorer plot data/sample.fasta --output results/figures/ --window 9 --format png

# Compare against Biopython and write docs/validation_report.md (needs the dev extra)
uv run python -m protein_explorer.crosscheck data/sample.fasta --output docs/validation_report.md

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
4. **Tables** (`dataset.py`). Builds pandas tables with one row per protein, indexed by FASTA identifier, and extracts the UniProt accession. Invalid sequences stay in the table with `NaN` features. A feature that can't be computed (e.g. molecular weight for a sequence containing `X`) is `NaN`, with the reason in `feature_errors`.
5. **Export** (`export.py`). Writes the tables as CSV files, plus a record of how they were produced.
6. **Figures** (`plots.py`). Static matplotlib figures at 300 dpi. Residues are coloured and grouped by chemical class (nonpolar, aromatic, polar, acidic, basic) using a colour-blind-safe palette.

## Outputs

`analyze` writes four files to the output directory:

| File | Contents |
|---|---|
| `features.csv` | One row per protein: accession, description, length, validation results, all features, warnings, and feature errors |
| `composition.csv` | One row per protein: the 20 amino acid frequencies |
| `summary.csv` | Count, mean, standard deviation, min, quartiles, and max of each numeric feature across valid proteins |
| `run_info.json` | Input file and its SHA-256 checksum, record counts, pH, pKa set, mass type, and software versions |

`features.csv` and `composition.csv` share the identifier column, so joining them gives a feature matrix with one row per protein and about 35 columns, ready for ML.

`plot` writes these figures:

| Figure | What it shows |
|---|---|
| `length_distribution` | Histogram of sequence lengths |
| `gravy_comparison` | GRAVY per protein, sorted; red is hydrophobic, blue is hydrophilic |
| `mw_vs_length` | Molecular weight against length, with a 110 Da per residue reference line |
| `composition_heatmap` | Residue frequencies, proteins × amino acids |
| `composition/<entry>` | One protein's composition, with the dataset mean marked for each residue |
| `hydropathy/<entry>` | Sliding-window Kyte–Doolittle profile along the sequence |

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
| Hydropathy profile | Mean Kyte–Doolittle value in a sliding window (odd width, default 9) plotted at each window's centre residue |

## Validation against Biopython

[docs/validation_report.md](docs/validation_report.md) compares every overlapping feature with Biopython 1.88. Biopython is a dev dependency only, used for checking.

- **Molecular weight, aromaticity, GRAVY:** match to rounding error (≤ 0.00002 Da, exact for the others).
- **pI and net charge:** differ only because of pKa conventions. Given Biopython's pKa set and its terminal-residue adjustments, our code reproduces Biopython's pI with a remaining difference of 0.0000.
- **Tests:** `tests/test_vs_biopython.py` checks each feature for all 8 proteins, with tolerances justified in the file's docstring.

## Findings 

- **Molecular weight:**  69.3 kDa was calculated for BSA, against the usual 66.4 kDa. The gap is the 24-residue signal peptide and propeptide in the precursor sequence.
- **Disulfides:** BSA's 17 disulfide bonds make the true mass about 34 Da lower than we calculate, because it is assumed that every Cys is reduced. This is small compared with the precursor effect.
- **Isoelectric point:** The pI values differ from Biopython's by up to 0.58 units (HBB: 7.32 vs 6.74). Up to 0.51 units comes from using EMBOSS rather than Biopython's pKa values, and up to 0.07 from Biopython adjusting terminal pKa by residue. Both rank the 8 proteins in the same order.
- **GRAVY and solubility:** insulin has a positive GRAVY (+0.19 for the precursor, +0.22 for the mature A+B chains), yet it is a soluble hormone. GRAVY is a weak predictor of solubility.
- **Sequence length:** insulin is listed as 110 residues because UniProt stores preproinsulin. The active hormone is 51 residues.
- **Biopython pI clamping:** Biopython only searches pH 4.05–12, so for poly-Asp (D×20) it returns 4.05, where the net charge is −9.2. Our 0–14 search finds 2.76 with the same pKa values.
- **Residue handling vs Biopython:** Biopython's GRAVY raises an error on X, B, or U, and it counts a terminal `*` as a residue. We skip unknown residues in averages and strip `*`.
- **Transmembrane threshold:** none of the 8 sample proteins sits in a membrane, yet lysozyme and β2-microglobulin exceed the Kyte–Doolittle transmembrane cutoff (1.6 at window 19). In both, the peak falls in the signal peptide, so hydropathy alone can't tell a signal peptide from a transmembrane helix.

## Project layout

```
src/protein_explorer/
  fasta.py         FASTA parser
  validation.py    sequence validation
  features.py      feature calculations
  dataset.py       pandas feature, composition, and summary tables
  export.py        CSV and run_info.json export
  plots.py         matplotlib figures
  crosscheck.py    comparison with Biopython (dev only)
  __main__.py      command-line interface (summarize, analyze, plot)
tests/             pytest suite and malformed-file fixtures
data/              sample FASTA and provenance
docs/              validation report
```
