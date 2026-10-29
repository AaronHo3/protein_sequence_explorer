"""Tests for protein sequence validation."""

from pathlib import Path

import pytest

from protein_explorer.fasta import read_fasta
from protein_explorer.validation import (
    NUCLEOTIDE_CHECK_MIN_LENGTH,
    ValidationResult,
    validate_sequence,
)

SAMPLE_FASTA = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"


def _has_warning(result: ValidationResult, text: str) -> bool:
    return any(text in warning.lower() for warning in result.warnings)


# --- Valid sequences ---------------------------------------------------------


def test_all_twenty_standard_residues_valid_without_warnings():
    result = validate_sequence("ACDEFGHIKLMNPQRSTVWY")
    assert result == ValidationResult(True, (), 0, 0, False, ())


def test_lowercase_is_valid_with_warning():
    result = validate_sequence("mktayi")
    assert result.is_valid is True
    assert result.invalid_residues == ()
    assert _has_warning(result, "lowercase")


@pytest.mark.parametrize("residue", list("BZJX"))
def test_ambiguous_residues_are_valid_and_counted(residue):
    result = validate_sequence(f"MKT{residue}AY{residue}I")
    assert result.is_valid is True
    assert result.ambiguous_count == 2
    assert _has_warning(result, "ambiguous")


@pytest.mark.parametrize("residue", ["U", "O"])
def test_nonstandard_residues_are_valid_and_counted(residue):
    result = validate_sequence(f"MKT{residue}AYI")
    assert result.is_valid is True
    assert result.nonstandard_count == 1
    assert _has_warning(result, "nonstandard")


def test_terminal_stop_is_valid_with_flag_and_warning():
    result = validate_sequence("MKTAYI*")
    assert result.is_valid is True
    assert result.has_terminal_stop is True
    assert result.invalid_residues == ()
    assert _has_warning(result, "stop")


def test_short_ga_peptide_does_not_trigger_nucleotide_check():
    result = validate_sequence("GAGA")
    assert result.is_valid is True
    assert not _has_warning(result, "nucleic")


# --- Invalid sequences -------------------------------------------------------


def test_internal_stop_is_invalid_with_position():
    result = validate_sequence("MKT*AYI")
    assert result.is_valid is False
    assert result.invalid_residues == ((4, "*"),)


def test_internal_and_terminal_stop_only_internal_reported():
    result = validate_sequence("MK*TAYI*")
    assert result.is_valid is False
    assert result.has_terminal_stop is True
    assert result.invalid_residues == ((3, "*"),)


@pytest.mark.parametrize("char", ["-", ".", "1", "#", "@"])
def test_invalid_characters_are_reported(char):
    result = validate_sequence(f"MK{char}TAYI")
    assert result.is_valid is False
    assert result.invalid_residues == ((3, char),)


def test_all_invalid_positions_reported_not_just_first():
    result = validate_sequence("M-KT.AY1I")
    assert result.invalid_residues == ((2, "-"), (5, "."), (8, "1"))


def test_invalid_character_reported_in_original_case():
    # Positions and characters refer to the input, not an uppercased copy
    result = validate_sequence("mk-tayi")
    assert result.invalid_residues == ((3, "-"),)


def test_empty_sequence_is_invalid():
    result = validate_sequence("")
    assert result.is_valid is False
    assert _has_warning(result, "empty")


def test_stop_codon_only_is_invalid():
    result = validate_sequence("*")
    assert result.is_valid is False


def test_all_ambiguous_sequence_is_invalid():
    result = validate_sequence("XXXXXXXX")
    assert result.is_valid is False
    assert result.ambiguous_count == 8


# --- Nucleotide detection ----------------------------------------------------


def test_dna_sequence_is_flagged_as_nucleic_acid():
    result = validate_sequence("ATGGCTAGCTAGGCTACGATCGATCGTAGCTAGCTAG")
    assert result.is_valid is False
    assert result.invalid_residues == ()  # every char is a legal amino acid code
    assert _has_warning(result, "nucleic")


def test_rna_sequence_is_flagged_as_nucleic_acid():
    result = validate_sequence("AUGGCUAGCUAGGCUACGAUCGAUCGUAGCUAGCUAG")
    assert result.is_valid is False
    assert _has_warning(result, "nucleic")


def test_nucleotide_check_boundary_length():
    just_below = "ACGT" * 5
    assert len(just_below) == NUCLEOTIDE_CHECK_MIN_LENGTH
    assert validate_sequence(just_below[:-1]).is_valid is True
    assert validate_sequence(just_below).is_valid is False


def test_terminal_stop_not_counted_toward_nucleotide_length():
    # 19 nucleotide chars + '*' should stay below the length threshold
    result = validate_sequence("ACGT" * 4 + "ACG" + "*")
    assert result.is_valid is True


def test_ala_gly_rich_protein_below_threshold_is_valid():
    # 80% A/C/G/T, 20% other residues: protein-like, not DNA
    result = validate_sequence("AGAGTKAGAGTKAGCGTKAGAGTK")
    assert result.is_valid is True


# --- Immutability ------------------------------------------------------------


def test_result_is_frozen():
    result = validate_sequence("MKT")
    with pytest.raises(AttributeError):
        result.is_valid = False  # type: ignore[misc]


# --- Integration: real UniProt data ------------------------------------------


def test_all_sample_records_are_valid():
    for record in read_fasta(SAMPLE_FASTA):
        result = validate_sequence(record.sequence)
        assert result.is_valid, f"{record.identifier}: {result.warnings}"
        assert result.invalid_residues == ()
