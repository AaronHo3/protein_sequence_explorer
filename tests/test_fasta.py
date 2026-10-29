"""Tests for the manual FASTA parser.

One test per row of the parser contract, plus integration tests against
real UniProt data in data/sample.fasta.
"""

from pathlib import Path

import pytest

from protein_explorer.fasta import (
    FastaFormatError,
    FastaRecord,
    parse_fasta_lines,
    read_fasta,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SAMPLE_FASTA = PROJECT_ROOT / "data" / "sample.fasta"
FIXTURES = Path(__file__).resolve().parent / "fixtures"


# --- Normal input -----------------------------------------------------------


def test_single_record():
    records = parse_fasta_lines([">seq1 a protein", "MKTAYI"])
    assert records == [FastaRecord("seq1", "a protein", "MKTAYI")]


def test_multiple_records_preserve_file_order():
    lines = [">b", "AAA", ">a", "CCC", ">c", "DDD"]
    assert [r.identifier for r in parse_fasta_lines(lines)] == ["b", "a", "c"]


def test_multiline_sequence_is_concatenated():
    lines = [">seq1 first protein", "MKT", "AYI", ">seq2", "GGG"]
    records = parse_fasta_lines(lines)
    assert records[0].sequence == "MKTAYI"
    assert records[1].sequence == "GGG"


def test_header_without_description_gives_empty_string():
    records = parse_fasta_lines([">P12345", "MKT"])
    assert records[0].identifier == "P12345"
    assert records[0].description == ""


def test_uniprot_identifier_kept_intact():
    header = ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1"
    record = parse_fasta_lines([header, "MALW"])[0]
    assert record.identifier == "sp|P01308|INS_HUMAN"
    assert record.description == "Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1"


# --- Whitespace handling -----------------------------------------------------


def test_blank_lines_are_ignored():
    lines = ["", ">seq1", "", "MKT", "   ", "AYI", "", ">seq2", "GGG", ""]
    records = parse_fasta_lines(lines)
    assert [r.sequence for r in records] == ["MKTAYI", "GGG"]


def test_trailing_whitespace_and_newlines_are_stripped():
    records = parse_fasta_lines([">seq1  \n", "MKT \t\n", "AYI\n"])
    assert records[0].sequence == "MKTAYI"


def test_crlf_line_endings_are_handled():
    records = read_fasta(FIXTURES / "crlf.fasta")
    assert [r.sequence for r in records] == ["MKTAYI", "GGG"]
    assert "\r" not in records[0].description


def test_internal_whitespace_in_sequence_is_removed():
    records = parse_fasta_lines([">seq1", "MKT AYI", "GG\tG"])
    assert records[0].sequence == "MKTAYIGGG"


# --- Biological data is not modified -----------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        "mktayi",  # lowercase
        "MKTXAYI",  # unknown residue
        "MKTAYI*",  # stop codon
        "MKBZJUO",  # ambiguous + nonstandard residues
        "MKT1AYI",  # stray digit
    ],
)
def test_residues_are_preserved_exactly(raw):
    assert parse_fasta_lines([">seq1", raw])[0].sequence == raw


# --- Malformed input ---------------------------------------------------------


def test_sequence_before_header_raises_with_line_number():
    with pytest.raises(FastaFormatError, match=r"(?i)line 1\b"):
        parse_fasta_lines(["MKT", ">seq1", "AAA"])


def test_sequence_before_header_reports_correct_line_after_blank_lines():
    with pytest.raises(FastaFormatError, match=r"(?i)line 3\b"):
        parse_fasta_lines(["", "", "MKT", ">seq1", "AAA"])


@pytest.mark.parametrize("header", [">", ">   "])
def test_header_without_identifier_raises(header):
    with pytest.raises(FastaFormatError, match=r"(?i)line 1\b"):
        parse_fasta_lines([header, "MKT"])


def test_header_with_no_sequence_raises_naming_identifier():
    with pytest.raises(FastaFormatError, match="seq2"):
        parse_fasta_lines([">seq1", "MKT", ">seq2", ">seq3", "GGG"])


def test_final_header_with_no_sequence_raises():
    with pytest.raises(FastaFormatError, match="seq2"):
        parse_fasta_lines([">seq1", "MKT", ">seq2"])


@pytest.mark.parametrize("lines", [[], [""], ["", "   ", "\n"]])
def test_empty_input_raises(lines):
    with pytest.raises(FastaFormatError):
        parse_fasta_lines(lines)


def test_duplicate_identifier_raises_naming_duplicate():
    with pytest.raises(FastaFormatError, match="seq1"):
        parse_fasta_lines([">seq1", "MKT", ">seq2", "GGG", ">seq1 again", "AAA"])


# --- read_fasta (file wrapper) -------------------------------------------------


def test_read_fasta_accepts_str_path():
    assert len(read_fasta(str(SAMPLE_FASTA))) == 8


def test_read_fasta_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_fasta(tmp_path / "does_not_exist.fasta")


def test_read_fasta_orphan_sequence_fixture_raises():
    with pytest.raises(FastaFormatError, match=r"(?i)line 1\b"):
        read_fasta(FIXTURES / "orphan_sequence.fasta")


def test_read_fasta_empty_record_fixture_raises():
    with pytest.raises(FastaFormatError, match="seq2"):
        read_fasta(FIXTURES / "empty_record.fasta")


# --- Integration: real UniProt data ------------------------------------------


@pytest.fixture(scope="module")
def sample_records() -> dict[str, FastaRecord]:
    return {r.identifier: r for r in read_fasta(SAMPLE_FASTA)}


def test_sample_has_eight_records(sample_records):
    assert len(sample_records) == 8


@pytest.mark.parametrize(
    ("identifier", "expected_length"),
    [
        ("sp|P01308|INS_HUMAN", 110),
        ("sp|P68871|HBB_HUMAN", 147),
        ("sp|P02769|ALBU_BOVIN", 607),
        ("sp|P00698|LYSC_CHICK", 147),
        ("sp|P42212|GFP_AEQVI", 238),
    ],
)
def test_sample_lengths_match_uniprot(sample_records, identifier, expected_length):
    assert len(sample_records[identifier].sequence) == expected_length


def test_sample_insulin_sequence_ends(sample_records):
    # Checks that the first and last wrapped lines were joined in the right order
    insulin = sample_records["sp|P01308|INS_HUMAN"].sequence
    assert insulin.startswith("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED")
    assert insulin.endswith("LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")
