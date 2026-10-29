# Writes analysis tables to CSV, plus a run_info.json recording how they were produced.

import hashlib
import json
import platform
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import pandas as pd

from protein_explorer.features import PKA_SET

FEATURES_FILE = "features.csv"
COMPOSITION_FILE = "composition.csv"
SUMMARY_FILE = "summary.csv"
RUN_INFO_FILE = "run_info.json"
# 6 decimals is finer than any estimate here (pI bisection tolerance is 1e-4)
FLOAT_FORMAT = "%.6f"

def file_sha256(path: Path) -> str:
    """SHA-256 of a file, so results can be traced to the exact input."""
    return hashlib.sha256(path.read_bytes()).hexdigest()

def build_run_info(input_path: Path, feature_table: pd.DataFrame, ph: float) -> dict[str, object]:
    """Metadata needed to reproduce or interpret a run."""
    return {
        "input_file": str(input_path),
        "input_sha256": file_sha256(input_path),
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
