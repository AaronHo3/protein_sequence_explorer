# Organises parsed records into pandas tables: one row per protein. Invalid records
# are kept as rows with NaN features, never dropped.

import re
from collections.abc import Callable

import pandas as pd

from protein_explorer import features
from protein_explorer.fasta import FastaRecord
from protein_explorer.validation import validate_sequence

# UniProt FASTA identifiers look like sp|P01308|INS_HUMAN or tr|A0A024R161|..., with
# optional isoform suffix (P04637-2)
UNIPROT_ID = re.compile(r"^(?:sp|tr)\|([A-Z0-9]+(?:-\d+)?)\|")

COUNT_COLUMNS = (
    "length", "ambiguous_count", "nonstandard_count",
    "unique_residue_types", "acidic_count", "basic_count",
)

def parse_accession(identifier: str) -> str | None:
    """Extract the UniProt accession from a FASTA identifier, or None if not UniProt."""
    match = UNIPROT_ID.match(identifier)
    return match.group(1) if match else None

# UniProt FASTA descriptions: "<name> OS=<organism> OX=<taxid> GN=<gene> PE=<n> SV=<n>"
UNIPROT_FIELDS = re.compile(r"\s(OS|OX|GN|PE|SV)=")

def parse_uniprot_description(description: str) -> dict[str, str | None]:
    """Protein name, organism, and gene from a UniProt-style description. Non-UniProt
    descriptions return the whole text as the name and None for the other fields."""
    parts = UNIPROT_FIELDS.split(" " + description)
    name = parts[0].strip() or None
    fields = dict(zip(parts[1::2], (v.strip() for v in parts[2::2]), strict=True))
    return {"protein_name": name, "organism": fields.get("OS"), "gene": fields.get("GN")}

def _feature_functions(ph: float) -> dict[str, Callable[[str], float]]:
    """Feature columns in output order. Net charge column is named after its pH."""
    return {
        "unique_residue_types": features.unique_residue_types,
        "molecular_weight": features.molecular_weight,
        "isoelectric_point": features.isoelectric_point,
        f"net_charge_ph{ph:g}": lambda seq: features.net_charge(seq, ph),
        "gravy": features.gravy,
        "hydrophobic_fraction": features.hydrophobic_fraction,
        "charged_fraction": features.charged_fraction,
        "aromatic_fraction": features.aromatic_fraction,
        "acidic_count": features.acidic_count,
        "basic_count": features.basic_count,
    }

def _compute_features(sequence: str, ph: float) -> tuple[dict[str, float], list[str]]:
    """Compute every feature; a feature that raises ValueError becomes NaN with a reason."""
    values: dict[str, float] = {}
    errors: list[str] = []
    for name, func in _feature_functions(ph).items():
        try:
            values[name] = func(sequence)
        except ValueError as error:
            values[name] = float("nan")
            errors.append(f"{name}: {error}")
    return values, errors

def _feature_row(record: FastaRecord, ph: float) -> dict[str, object]:
    """One table row: metadata, validation results, and features (NaN if invalid)."""
    result = validate_sequence(record.sequence)
    if result.is_valid:
        values, errors = _compute_features(record.sequence, ph)
    else:
        values = {name: float("nan") for name in _feature_functions(ph)}
        errors = ["not computed: sequence failed validation"]
    return {
        "identifier": record.identifier,
        "accession": parse_accession(record.identifier),
        "description": record.description,
        **parse_uniprot_description(record.description),
        "length": len(record.sequence.removesuffix("*")),
        "is_valid": result.is_valid,
        "ambiguous_count": result.ambiguous_count,
        "nonstandard_count": result.nonstandard_count,
        "has_terminal_stop": result.has_terminal_stop,
        **values,
        "warnings": "; ".join(result.warnings),
        "feature_errors": "; ".join(errors),
    }

def build_feature_table(records: list[FastaRecord], ph: float = 7.0) -> pd.DataFrame:
    """One row per record, indexed by identifier, in input order."""
    table = pd.DataFrame([_feature_row(record, ph) for record in records])
    table = table.set_index("identifier")
    # Nullable integers keep counts as ints even when some rows are NaN
    return table.astype({column: "Int64" for column in COUNT_COLUMNS})

def build_composition_table(records: list[FastaRecord]) -> pd.DataFrame:
    """Amino acid frequencies: one row per record, 20 residue columns. Invalid rows are NaN."""
    rows = {
        record.identifier: (
            features.residue_frequencies(record.sequence)
            if validate_sequence(record.sequence).is_valid
            else dict.fromkeys(features.AMINO_ACIDS, float("nan"))
        )
        for record in records
    }
    table = pd.DataFrame.from_dict(rows, orient="index", columns=list(features.AMINO_ACIDS))
    table.index.name = "identifier"
    table.columns.name = "residue"
    return table

def summarize_features(table: pd.DataFrame) -> pd.DataFrame:
    """Descriptive statistics across valid proteins: one row per numeric feature."""
    numeric = table.loc[table["is_valid"]].select_dtypes("number")
    summary = numeric.astype(float).describe().T
    summary.index.name = "feature"
    return summary

def filter_features(
    table: pd.DataFrame,
    query: str = "",
    valid_only: bool = False,
    length_range: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """Rows matching a case-insensitive search (identifier or description), validity,
    and an inclusive length range. Returns a new table; the input is not modified."""
    mask = pd.Series(True, index=table.index)
    if query.strip():
        needle = query.strip().lower()
        haystack = table.index.str.lower() + " " + table["description"].fillna("").str.lower()
        mask &= haystack.str.contains(needle, regex=False)
    if valid_only:
        mask &= table["is_valid"].astype(bool)
    if length_range is not None:
        low, high = length_range
        mask &= table["length"].between(low, high).fillna(False).astype(bool)
    return table.loc[mask]


def build_count_table(records: list[FastaRecord]) -> pd.DataFrame:
    """Residue counts: one row per record, 20 residue columns. Invalid rows are NA."""
    rows = {
        record.identifier: (
            features.residue_counts(record.sequence)
            if validate_sequence(record.sequence).is_valid
            else dict.fromkeys(features.AMINO_ACIDS, pd.NA)
        )
        for record in records
    }
    table = pd.DataFrame.from_dict(rows, orient="index", columns=list(features.AMINO_ACIDS))
    table.index.name = "identifier"
    table.columns.name = "residue"
    return table.astype("Int64")

def dataset_composition(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """Dataset-level amino acid frequencies, two ways:

    - residue_weighted: pooled counts / pooled length. Every residue counts equally, so
      long proteins dominate (what you get by concatenating all sequences).
    - protein_weighted: mean of each protein's own frequencies. Every protein counts
      equally, regardless of length.

    Rows with missing counts (invalid sequences) are excluded from both.
    """
    valid = counts.dropna(how="any").astype(float)
    valid_lengths = lengths.loc[valid.index].astype(float)
    if valid.empty:
        return pd.DataFrame(columns=["residue_weighted", "protein_weighted"], dtype=float)
    pooled = valid.sum() / valid_lengths.sum()
    per_protein = valid.div(valid_lengths, axis=0).mean()
    result = pd.DataFrame({"residue_weighted": pooled, "protein_weighted": per_protein})
    result.index.name = "residue"
    return result
