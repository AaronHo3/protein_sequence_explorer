"""Tests for residue reference data."""

import pytest

from protein_explorer.features import AMINO_ACIDS, KYTE_DOOLITTLE
from protein_explorer.residues import (
    CHEMICAL_CLASS,
    NAMES,
    RESIDUE_ORDER,
    THREE_LETTER,
    residue_info,
    side_chain_charge,
)


def test_every_standard_residue_has_codes_names_and_class():
    for code in AMINO_ACIDS:
        assert code in THREE_LETTER and code in NAMES and code in CHEMICAL_CLASS


def test_residue_order_is_a_permutation_of_the_alphabet():
    assert sorted(RESIDUE_ORDER) == list(AMINO_ACIDS)


@pytest.mark.parametrize(("code", "three", "name"), [
    ("L", "Leu", "Leucine"), ("W", "Trp", "Tryptophan"), ("U", "Sec", "Selenocysteine"),
    ("X", "Xaa", "Unknown residue"),
])
def test_three_letter_codes(code, three, name):
    info = residue_info(code)
    assert (info.three_letter, info.name) == (three, name)


def test_residue_info_is_case_insensitive_and_uses_feature_scale():
    info = residue_info("i")
    assert info.code == "I"
    assert info.hydropathy == KYTE_DOOLITTLE["I"] == 4.5
    assert info.is_standard


def test_nonstandard_and_unknown_symbols_do_not_raise():
    assert residue_info("X").hydropathy is None
    assert not residue_info("X").is_standard
    assert residue_info("-").three_letter == "?"


@pytest.mark.parametrize(("code", "pka", "sign"), [
    ("K", 10.8, 1), ("R", 12.5, 1), ("H", 6.5, 1), ("D", 3.9, -1), ("E", 4.1, -1),
    ("C", 8.5, -1), ("Y", 10.1, -1), ("U", 5.2, -1), ("A", None, 0),
])
def test_side_chain_pka_matches_emboss_set(code, pka, sign):
    info = residue_info(code)
    assert info.side_chain_pka == pka
    assert info.charge_sign == sign


def test_side_chain_charge_half_at_pka_and_zero_for_neutral():
    assert side_chain_charge("K", 10.8) == pytest.approx(0.5)
    assert side_chain_charge("D", 3.9) == pytest.approx(-0.5)
    assert side_chain_charge("A", 7.0) == 0.0


def test_histidine_mostly_neutral_at_physiological_ph():
    assert side_chain_charge("H", 7.4) == pytest.approx(0.112, abs=1e-3)
