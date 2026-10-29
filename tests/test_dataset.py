"""Tests for building feature, composition, and summary tables."""

from pathlib import Path

import pandas as pd
import pytest

from protein_explorer.dataset import (
    build_composition_table,
    build_feature_table,
    parse_accession,
    summarize_features,
)
from protein_explorer.fasta import FastaRecord, read_fasta
from protein_explorer.features import AMINO_ACIDS, gravy, molecular_weight

SAMPLE_FASTA = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"

GOOD = FastaRecord("sp|P00001|GOOD_TEST", "valid protein", "MKTAYIAKQRDEHW")
GAPPED = FastaRecord("gapped", "alignment gap", "MKT-AYI")
AMBIGUOUS = FastaRecord("ambiguous", "contains X", "MKTXAYILLV")


@pytest.fixture(scope="module")
def sample_records() -> list[FastaRecord]:
    return read_fasta(SAMPLE_FASTA)


@pytest.fixture(scope="module")
def sample_table(sample_records) -> pd.DataFrame:
    return build_feature_table(sample_records)


# --- Accession parsing -------------------------------------------------------


@pytest.mark.parametrize(
    ("identifier", "expected"),
    [
        ("sp|P01308|INS_HUMAN", "P01308"),
        ("tr|A0A024R161|A0A024R161_HUMAN", "A0A024R161"),
        ("sp|P04637-2|P53_HUMAN", "P04637-2"),  # isoform
        ("P01308", None),
        ("gi|12345|ref|NP_000198.1|", None),
        ("my_protein", None),
    ],
)
def test_parse_accession(identifier, expected):
    assert parse_accession(identifier) == expected


# --- Feature table: shape and content ----------------------------------------


def test_feature_table_has_one_row_per_record_in_order(sample_records, sample_table):
    assert list(sample_table.index) == [r.identifier for r in sample_records]
    assert sample_table.index.name == "identifier"


def test_feature_table_values_match_feature_functions(sample_records, sample_table):
    insulin = next(r for r in sample_records if "INS_HUMAN" in r.identifier)
    row = sample_table.loc[insulin.identifier]
    assert row["accession"] == "P01308"
    assert row["length"] == 110
    assert row["molecular_weight"] == pytest.approx(molecular_weight(insulin.sequence))
    assert row["gravy"] == pytest.approx(gravy(insulin.sequence))


def test_feature_table_all_sample_records_valid(sample_table):
    assert sample_table["is_valid"].all()
    assert (sample_table["feature_errors"] == "").all()
    assert not sample_table.drop(columns=["accession"]).isna().any().any()


def test_count_columns_are_integers(sample_table):
    for column in ["length", "acidic_count", "basic_count", "unique_residue_types"]:
        assert pd.api.types.is_integer_dtype(sample_table[column]), column


def test_net_charge_column_named_after_ph():
    assert "net_charge_ph7" in build_feature_table([GOOD]).columns
    assert "net_charge_ph7.4" in build_feature_table([GOOD], ph=7.4).columns


def test_ph_changes_net_charge():
    at_5 = build_feature_table([GOOD], ph=5.0)["net_charge_ph5"].iloc[0]
    at_9 = build_feature_table([GOOD], ph=9.0)["net_charge_ph9"].iloc[0]
    assert at_5 > at_9


# --- Feature table: invalid and partial records ------------------------------


def test_invalid_record_kept_with_nan_features():
    table = build_feature_table([GOOD, GAPPED])
    assert len(table) == 2
    row = table.loc["gapped"]
    assert row["is_valid"] == False  # noqa: E712 (numpy bool)
    assert pd.isna(row["molecular_weight"])
    assert pd.isna(row["gravy"])
    assert "failed validation" in row["feature_errors"]
    assert "invalid" in row["warnings"]


def test_invalid_record_still_reports_length_and_validation_counts():
    row = build_feature_table([GAPPED]).loc["gapped"]
    assert row["length"] == 7
    assert row["ambiguous_count"] == 0


def test_uncomputable_feature_is_nan_with_reason():
    # X has no mass, but GRAVY and fractions are still defined
    row = build_feature_table([AMBIGUOUS]).loc["ambiguous"]
    assert row["is_valid"] == True  # noqa: E712
    assert pd.isna(row["molecular_weight"])
    assert row["gravy"] == pytest.approx(gravy(AMBIGUOUS.sequence))
    assert "molecular_weight" in row["feature_errors"]
    assert row["ambiguous_count"] == 1


def test_non_uniprot_identifier_has_missing_accession():
    assert pd.isna(build_feature_table([GAPPED]).loc["gapped", "accession"])


# --- Composition table -------------------------------------------------------


def test_composition_table_shape_and_columns(sample_records):
    table = build_composition_table(sample_records)
    assert table.shape == (8, 20)
    assert tuple(table.columns) == AMINO_ACIDS
    assert table.index.name == "identifier"


def test_composition_rows_sum_to_one_for_sample(sample_records):
    sums = build_composition_table(sample_records).sum(axis=1)
    assert sums.to_numpy() == pytest.approx([1.0] * 8)


def test_composition_invalid_record_is_nan_row():
    table = build_composition_table([GOOD, GAPPED])
    assert table.loc["gapped"].isna().all()
    assert not table.loc["sp|P00001|GOOD_TEST"].isna().any()


def test_composition_joins_with_feature_table(sample_records, sample_table):
    # Shared identifier index gives the full ML feature matrix
    matrix = sample_table.join(build_composition_table(sample_records))
    assert matrix.shape == (8, sample_table.shape[1] + 20)


# --- Summary statistics ------------------------------------------------------


def test_summary_has_one_row_per_numeric_feature(sample_table):
    summary = summarize_features(sample_table)
    assert "molecular_weight" in summary.index
    assert "description" not in summary.index
    assert {"mean", "std", "min", "max"} <= set(summary.columns)


def test_summary_values_match_table(sample_table):
    summary = summarize_features(sample_table)
    assert summary.loc["length", "min"] == 110
    assert summary.loc["length", "max"] == 607
    assert summary.loc["gravy", "mean"] == pytest.approx(sample_table["gravy"].mean())


def test_summary_excludes_invalid_records():
    table = build_feature_table([GOOD, GAPPED])
    assert summarize_features(table).loc["length", "count"] == 1
