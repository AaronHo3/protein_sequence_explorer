"""Tests for composition features."""

from pathlib import Path

import pytest

from protein_explorer.fasta import read_fasta
from protein_explorer.features import (
    AMINO_ACIDS,
    residue_counts,
    residue_frequencies,
    sequence_length,
    unique_residue_types,
)

SAMPLE_FASTA = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"
INSULIN_SIGNAL_START = "MALWMRLLPL"  # first 10 residues of preproinsulin (P01308)


# --- Worked example ----------------------------------------------------------


def test_worked_example_counts():
    counts = residue_counts(INSULIN_SIGNAL_START)
    nonzero = {aa: n for aa, n in counts.items() if n}
    assert nonzero == {"A": 1, "L": 4, "M": 2, "P": 1, "R": 1, "W": 1}


def test_worked_example_frequencies():
    freqs = residue_frequencies(INSULIN_SIGNAL_START)
    assert freqs["L"] == pytest.approx(0.4)
    assert freqs["M"] == pytest.approx(0.2)
    assert freqs["C"] == 0.0


def test_worked_example_unique_types():
    assert unique_residue_types(INSULIN_SIGNAL_START) == 6


# --- Fixed shape -------------------------------------------------------------


def test_amino_acids_are_twenty_in_alphabetical_order():
    assert AMINO_ACIDS == tuple("ACDEFGHIKLMNPQRSTVWY")


@pytest.mark.parametrize("func", [residue_counts, residue_frequencies])
def test_all_twenty_keys_present_in_order(func):
    assert tuple(func("AAA")) == AMINO_ACIDS


def test_missing_residues_have_zero_count():
    counts = residue_counts("AAA")
    assert counts["A"] == 3
    assert sum(counts.values()) == 3


# --- Frequencies and denominators --------------------------------------------


def test_frequencies_sum_to_one_for_standard_sequence():
    assert sum(residue_frequencies("ACDEFGHIKLMNPQRSTVWYAAK").values()) == pytest.approx(1.0)


def test_ambiguous_residue_counts_toward_length_but_has_no_key():
    freqs = residue_frequencies("MKTXAYI")
    assert sequence_length("MKTXAYI") == 7
    assert "X" not in freqs
    assert sum(freqs.values()) == pytest.approx(6 / 7)


def test_nonstandard_residue_counts_toward_length():
    assert sequence_length("MKUTAYI") == 7
    assert sum(residue_counts("MKUTAYI").values()) == 6


# --- Stop codon and case -----------------------------------------------------


def test_terminal_stop_excluded_from_length():
    assert sequence_length("MKTAYI*") == 6
    assert residue_frequencies("MKTAYI*") == residue_frequencies("MKTAYI")


def test_lowercase_matches_uppercase():
    assert residue_counts("mktayi") == residue_counts("MKTAYI")
    assert unique_residue_types("mktayi") == 6


# --- Empty input -------------------------------------------------------------


@pytest.mark.parametrize("sequence", ["", "*"])
@pytest.mark.parametrize(
    "func", [sequence_length, residue_counts, residue_frequencies, unique_residue_types]
)
def test_empty_sequence_raises(func, sequence):
    with pytest.raises(ValueError, match="empty"):
        func(sequence)


# --- Integration: real UniProt data ------------------------------------------


@pytest.fixture(scope="module")
def sample() -> dict[str, str]:
    return {r.identifier: r.sequence for r in read_fasta(SAMPLE_FASTA)}


def test_insulin_composition(sample):
    insulin = sample["sp|P01308|INS_HUMAN"]
    counts = residue_counts(insulin)
    assert counts["C"] == 6  # three disulfide bonds
    assert counts["L"] == 20
    assert unique_residue_types(insulin) == 20


def test_bsa_cysteine_count(sample):
    # 17 disulfide bonds plus the free Cys34
    assert residue_counts(sample["sp|P02769|ALBU_BOVIN"])["C"] == 35


def test_sample_frequencies_sum_to_one(sample):
    # Reviewed Swiss-Prot entries contain only standard residues
    for identifier, sequence in sample.items():
        assert sum(residue_frequencies(sequence).values()) == pytest.approx(1.0), identifier
