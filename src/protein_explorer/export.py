# Writes analysis tables to CSV, plus a run_info.json recording how they were produced.

import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from protein_explorer.dataset import summarize_features
from protein_explorer.fasta import FastaRecord
from protein_explorer.features import PKA_SET

FEATURES_FILE = "features.csv"
COMPOSITION_FILE = "composition.csv"
SUMMARY_FILE = "summary.csv"
RUN_INFO_FILE = "run_info.json"
# 6 decimals is finer than any estimate here (pI bisection tolerance is 1e-4)
FLOAT_FORMAT = "%.6f"

def build_run_info(input_path: Path, feature_table: pd.DataFrame, ph: float) -> dict[str, object]:
    """Metadata needed to reproduce or interpret a run from a file on disk."""
    return build_run_info_from_bytes(str(input_path), input_path.read_bytes(), feature_table, ph)

def build_run_info_from_bytes(
    input_name: str, input_bytes: bytes, feature_table: pd.DataFrame, ph: float
) -> dict[str, object]:
    """Run metadata for in-memory input (e.g. an uploaded file in the Streamlit app)."""
    return {
        "input_file": input_name,
        "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
        "records": len(feature_table),
        "valid_records": int(feature_table["is_valid"].sum()),
        "net_charge_ph": ph,
        "pka_set": PKA_SET,
        "mass_type": "average",
        "protein_explorer_version": version("protein-explorer"),
        "pandas_version": pd.__version__,
        "python_version": platform.python_version(),
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }

def write_results(
    output_dir: Path,
    feature_table: pd.DataFrame,
    composition_table: pd.DataFrame,
    summary_table: pd.DataFrame,
    run_info: dict[str, object],
) -> list[Path]:
    """Write all outputs to output_dir (created if needed). Returns the written paths."""
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = {
        FEATURES_FILE: feature_table,
        COMPOSITION_FILE: composition_table,
        SUMMARY_FILE: summary_table,
    }
    paths = []
    for filename, table in outputs.items():
        path = output_dir / filename
        table.to_csv(path, float_format=FLOAT_FORMAT)
        paths.append(path)
    run_info_path = output_dir / RUN_INFO_FILE
    run_info_path.write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8")
    paths.append(run_info_path)
    return paths

FASTA_LINE_WIDTH = 60

def records_to_fasta(records: list[FastaRecord], line_width: int = FASTA_LINE_WIDTH) -> str:
    """FASTA text for records, preserving identifiers, descriptions, and residues exactly
    (only line wrapping changes)."""
    if line_width < 1:
        raise ValueError("line_width must be positive.")
    blocks = []
    for record in records:
        header = f">{record.identifier} {record.description}".rstrip()
        lines = [record.sequence[i:i + line_width] for i in range(0, len(record.sequence), line_width)]
        blocks.append("\n".join([header, *lines]))
    return "\n".join(blocks) + "\n"

def _markdown_table(table: pd.DataFrame, float_format: str = "{:.3f}") -> str:
    def cell(value: object) -> str:
        if isinstance(value, float):
            return "" if pd.isna(value) else float_format.format(value)
        text = "" if value is None or value is pd.NA else str(value)
        return text.replace("|", "\\|")  # UniProt IDs contain '|', a Markdown column separator
    header = "| " + " | ".join(str(c) for c in table.columns) + " |"
    divider = "|" + "---|" * len(table.columns)
    rows = ["| " + " | ".join(cell(v) for v in row) + " |" for row in table.itertuples(index=False)]
    return "\n".join([header, divider, *rows])

def build_analysis_report(
    feature_table: pd.DataFrame, run_info: dict[str, object], composition: pd.DataFrame
) -> str:
    """Markdown summary of a run: provenance, per-protein features, dataset statistics,
    dataset composition, and the assumptions behind every estimate."""
    charge = next(c for c in feature_table.columns if c.startswith("net_charge_ph"))
    valid = feature_table[feature_table["is_valid"].astype(bool)]
    per_protein = valid.reset_index()[[
        "identifier", "length", "molecular_weight", "isoelectric_point", charge, "gravy",
        "aromatic_fraction", "hydrophobic_fraction", "charged_fraction",
    ]]
    invalid = feature_table[~feature_table["is_valid"].astype(bool)]
    stats = summarize_features(feature_table).reset_index()[["feature", "count", "mean", "std", "min", "max"]]
    comp = (composition * 100).reset_index()
    sections = [
        "# Protein Sequence Explorer analysis report",
        "",
        f"- Input: `{run_info['input_file']}` (SHA-256 `{run_info['input_sha256']}`)",
        f"- Records: {run_info['records']} ({run_info['valid_records']} valid)",
        f"- Net charge pH: {run_info['net_charge_ph']}; pKa set: {run_info['pka_set']}; "
        f"masses: {run_info['mass_type']}",
        f"- protein_explorer {run_info['protein_explorer_version']}, created {run_info['created_utc']}",
        "",
        "## Per-protein features",
        "",
        _markdown_table(per_protein),
        "",
    ]
    if not invalid.empty:
        sections += ["## Sequences that failed validation", "",
                     _markdown_table(invalid.reset_index()[["identifier", "warnings"]]), ""]
    sections += [
        "## Summary across valid proteins", "", _markdown_table(stats), "",
        "## Dataset composition (%)", "",
        "Residue-weighted pools all residues (long proteins dominate); protein-weighted "
        "averages each protein's own frequencies.", "",
        _markdown_table(comp, "{:.2f}"), "",
        "## Assumptions", "",
        "- All values are estimates from sequence alone, not measurements.",
        "- Molecular weight uses average residue masses plus one water; no modifications, "
        "reduced cysteines, precursor sequences as written.",
        "- pI and net charge use Henderson–Hasselbalch with the pKa set above, including both "
        "termini; structural pKa shifts are ignored.",
        "- GRAVY and hydropathy use the Kyte–Doolittle scale; ambiguous residues are skipped "
        "in averages but count toward length.",
        "",
    ]
    return "\n".join(sections)
