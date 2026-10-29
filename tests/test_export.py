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


# --- FASTA export, analysis report, figure bundle ----------------------------


def test_fasta_export_round_trips_sample_exactly():
    from protein_explorer.export import records_to_fasta
    assert records_to_fasta(read_fasta(SAMPLE_FASTA)) == SAMPLE_FASTA.read_text()


def test_fasta_export_preserves_residues_and_wraps():
    from protein_explorer.export import records_to_fasta
    from protein_explorer.fasta import FastaRecord
    text = records_to_fasta([FastaRecord("x", "", "mktXAY*")], line_width=3)
    assert text == ">x\nmkt\nXAY\n*\n"
    with pytest.raises(ValueError):
        records_to_fasta([], line_width=0)


def test_analysis_report_contents():
    from protein_explorer.dataset import build_count_table, dataset_composition
    from protein_explorer.export import build_analysis_report
    records = read_fasta(SAMPLE_FASTA)
    features = build_feature_table(records)
    comp = dataset_composition(build_count_table(records), features["length"])
    report = build_analysis_report(features, build_run_info(SAMPLE_FASTA, features, 7.0), comp)
    assert "## Per-protein features" in report and "## Assumptions" in report
    assert r"sp\|P01308\|INS_HUMAN" in report  # pipes escaped for Markdown tables
    assert "SHA-256" in report and "EMBOSS" in report


def test_static_figure_bundle_contains_all_figures():
    import io
    import zipfile
    from protein_explorer.visualization.export import static_figure_bundle
    names = zipfile.ZipFile(io.BytesIO(static_figure_bundle(read_fasta(SAMPLE_FASTA)))).namelist()
    assert len(names) == 20
    assert "hydropathy/INS_HUMAN.png" in names and "mw_vs_length.png" in names
