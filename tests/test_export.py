"""Tests for CSV export and run metadata."""

import hashlib
import json
from pathlib import Path

import pandas as pd
import pytest

from protein_explorer.dataset import (
    build_composition_table,
    build_feature_table,
    summarize_features,
)
from protein_explorer.export import build_run_info, write_results
from protein_explorer.fasta import read_fasta

SAMPLE_FASTA = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"


@pytest.fixture
def written(tmp_path) -> dict[str, Path]:
    records = read_fasta(SAMPLE_FASTA)
    features = build_feature_table(records)
    paths = write_results(
        tmp_path / "nested" / "results",
        features,
        build_composition_table(records),
        summarize_features(features),
        build_run_info(SAMPLE_FASTA, features, ph=7.0),
    )
    return {path.name: path for path in paths}


def test_all_output_files_written(written):
    assert set(written) == {"features.csv", "composition.csv", "summary.csv", "run_info.json"}
    assert all(path.exists() for path in written.values())


def test_features_csv_round_trips(written):
    table = pd.read_csv(written["features.csv"], index_col="identifier")
    assert len(table) == 8
    assert table.loc["sp|P02769|ALBU_BOVIN", "molecular_weight"] == pytest.approx(
        69292.68, abs=0.01
    )


def test_composition_csv_has_twenty_residue_columns(written):
    table = pd.read_csv(written["composition.csv"], index_col="identifier")
    assert list(table.columns) == list("ACDEFGHIKLMNPQRSTVWY")


def test_floats_written_with_fixed_precision(written):
    first_row = written["features.csv"].read_text().splitlines()[1]
    assert "16238.493780" in first_row


def test_run_info_records_provenance(written):
    info = json.loads(written["run_info.json"].read_text())
    assert info["input_sha256"] == hashlib.sha256(SAMPLE_FASTA.read_bytes()).hexdigest()
    assert info["records"] == 8
    assert info["valid_records"] == 8
    assert info["net_charge_ph"] == 7.0
    assert "EMBOSS" in info["pka_set"]
    assert info["protein_explorer_version"]
