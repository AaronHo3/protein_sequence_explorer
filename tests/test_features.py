"""Tests for composition features."""

from pathlib import Path

import pytest

from protein_explorer.fasta import read_fasta
from protein_explorer.features import (
    AMINO_ACIDS,
    KYTE_DOOLITTLE,
    acidic_count,
    aromatic_fraction,
    basic_count,
    charged_fraction,
    gravy,
    hydrophobic_fraction,
    isoelectric_point,
    molecular_weight,
    net_charge,
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


# --- Biochemical fractions and counts ----------------------------------------


def test_worked_example_biochemical_features():
    # M A L W M R L L P L: 8 hydrophobic (M,A,L,W), 1 charged (R), 1 aromatic (W)
    assert hydrophobic_fraction(INSULIN_SIGNAL_START) == pytest.approx(0.8)
    assert charged_fraction(INSULIN_SIGNAL_START) == pytest.approx(0.1)
    assert aromatic_fraction(INSULIN_SIGNAL_START) == pytest.approx(0.1)
    assert acidic_count(INSULIN_SIGNAL_START) == 0
    assert basic_count(INSULIN_SIGNAL_START) == 1


def test_histidine_is_basic_but_not_charged():
    assert basic_count("HHH") == 3
    assert charged_fraction("HHH") == 0.0


def test_histidine_is_not_aromatic():
    # Matches Biopython's aromaticity definition (F, W, Y only)
    assert aromatic_fraction("HFWY") == pytest.approx(0.75)


def test_glycine_and_cysteine_are_not_hydrophobic():
    assert hydrophobic_fraction("GGCC") == 0.0


def test_fractions_use_full_length_as_denominator():
    # X counts toward length, so 3 hydrophobic of 4 residues
    assert hydrophobic_fraction("LLLX") == pytest.approx(0.75)


def test_fractions_ignore_terminal_stop_and_case():
    assert hydrophobic_fraction("llk*") == pytest.approx(2 / 3)


# --- GRAVY -------------------------------------------------------------------


def test_kyte_doolittle_covers_exactly_the_standard_residues():
    assert tuple(sorted(KYTE_DOOLITTLE)) == AMINO_ACIDS


def test_gravy_worked_example():
    # (1.9+1.8+3.8-0.9+1.9-4.5+3.8+3.8-1.6+3.8) / 10
    assert gravy(INSULIN_SIGNAL_START) == pytest.approx(1.38)


def test_gravy_single_residues_match_scale():
    assert gravy("I") == pytest.approx(4.5)
    assert gravy("R") == pytest.approx(-4.5)


def test_gravy_excludes_ambiguous_residues_from_mean():
    assert gravy("IIX") == pytest.approx(4.5)
    assert gravy("IIU") == pytest.approx(4.5)


def test_gravy_all_ambiguous_raises():
    with pytest.raises(ValueError, match="no standard residues"):
        gravy("XXX")


@pytest.mark.parametrize(
    "func",
    [hydrophobic_fraction, charged_fraction, aromatic_fraction, acidic_count, basic_count, gravy],
)
def test_biochemical_features_reject_empty_sequence(func):
    with pytest.raises(ValueError):
        func("*")


# --- Biochemical features on real data ---------------------------------------


def test_calmodulin_is_strongly_acidic(sample):
    # EF-hand Ca2+-binding loops are rich in Asp/Glu
    calmodulin = sample["sp|P0DP23|CALM1_HUMAN"]
    assert acidic_count(calmodulin) == 38
    assert basic_count(calmodulin) == 15


def test_lysozyme_is_basic(sample):
    # Cationic enzyme that binds anionic bacterial cell walls
    lysozyme = sample["sp|P00698|LYSC_CHICK"]
    assert basic_count(lysozyme) > 2 * acidic_count(lysozyme)


@pytest.mark.parametrize(
    ("identifier", "expected_gravy", "expected_aromaticity"),
    [
        # Reference values from Biopython ProteinAnalysis (rounded to 3 dp)
        ("sp|P01308|INS_HUMAN", 0.193, 0.082),
        ("sp|P04637|P53_HUMAN", -0.756, 0.061),
        ("sp|P42212|GFP_AEQVI", -0.521, 0.105),
        ("sp|P68871|HBB_HUMAN", 0.014, 0.088),
    ],
)
def test_sample_gravy_and_aromaticity_match_reference(
    sample, identifier, expected_gravy, expected_aromaticity
):
    assert gravy(sample[identifier]) == pytest.approx(expected_gravy, abs=5e-4)
    assert aromatic_fraction(sample[identifier]) == pytest.approx(expected_aromaticity, abs=5e-4)


# --- Molecular weight --------------------------------------------------------


def test_molecular_weight_diglycine():
    # 2 x Gly residue (57.0513) + one water (18.0153); textbook diglycine = 132.12 Da
    assert molecular_weight("GG") == pytest.approx(132.1179, abs=1e-3)


def test_molecular_weight_single_residue_is_free_amino_acid():
    # Free glycine = 75.07 Da
    assert molecular_weight("G") == pytest.approx(75.067, abs=1e-2)


def test_molecular_weight_includes_selenocysteine():
    assert molecular_weight("GU") - molecular_weight("GC") == pytest.approx(
        150.0379 - 103.1429, abs=1e-3
    )


@pytest.mark.parametrize("residue", list("BZJX"))
def test_molecular_weight_rejects_residues_without_mass(residue):
    with pytest.raises(ValueError, match=residue):
        molecular_weight(f"MK{residue}T")


def test_molecular_weight_ignores_terminal_stop_and_case():
    assert molecular_weight("gg*") == pytest.approx(molecular_weight("GG"))


# --- Net charge --------------------------------------------------------------


def test_net_charge_worked_example():
    # pH 7: Nterm +0.976, K +1.000, H +0.240, Cterm -1.000, D -0.999
    assert net_charge("GDKH", 7.0) == pytest.approx(0.217, abs=1e-3)


def test_net_charge_at_extreme_ph():
    # pH 0: all basic groups and N-terminus protonated, acids neutral
    assert net_charge("KRHD", 0.0) == pytest.approx(4.0, abs=1e-3)
    # pH 14: acids fully deprotonated (-2), but Arg (pKa 12.5) is still ~3% protonated
    assert net_charge("KRHD", 14.0) == pytest.approx(-2.0 + 1 / (1 + 10**1.5), abs=1e-3)


def test_group_is_half_charged_at_its_pka():
    # Lysine side chain contributes +0.5 at pH 10.8; termini nearly cancel to -1
    assert net_charge("GKG", 10.8) - net_charge("GGG", 10.8) == pytest.approx(0.5, abs=1e-6)


def test_net_charge_decreases_with_ph():
    charges = [net_charge("MKTAYIAKQRDEH", ph) for ph in range(15)]
    assert charges == sorted(charges, reverse=True)


def test_net_charge_default_ph_is_seven():
    assert net_charge("GDKH") == net_charge("GDKH", 7.0)


# --- Isoelectric point -------------------------------------------------------


def test_isoelectric_point_diglycine_is_midpoint_of_termini():
    # Only termini ionise, so pI = (8.6 + 3.6) / 2
    assert isoelectric_point("GG") == pytest.approx(6.1, abs=1e-3)


def test_net_charge_is_zero_at_isoelectric_point():
    sequence = "MKTAYIAKQRDEHCY"
    assert net_charge(sequence, isoelectric_point(sequence)) == pytest.approx(0.0, abs=1e-3)


def test_isoelectric_point_rejects_empty_sequence():
    with pytest.raises(ValueError):
        isoelectric_point("*")


# --- Equation-based estimates on real data -----------------------------------


@pytest.mark.parametrize(
    ("identifier", "expected_mw"),
    [
        # Biopython ProteinAnalysis.molecular_weight() (average masses)
        ("sp|P01308|INS_HUMAN", 11980.79),
        ("sp|P02769|ALBU_BOVIN", 69292.68),
        ("sp|P42212|GFP_AEQVI", 26885.99),
    ],
)
def test_sample_molecular_weight_matches_reference(sample, identifier, expected_mw):
    assert molecular_weight(sample[identifier]) == pytest.approx(expected_mw, abs=0.01)


def test_sample_pi_ordering_reflects_biology(sample):
    calmodulin = isoelectric_point(sample["sp|P0DP23|CALM1_HUMAN"])
    lysozyme = isoelectric_point(sample["sp|P00698|LYSC_CHICK"])
    assert calmodulin < 4.5  # acidic Ca2+-binding protein
    assert lysozyme > 9.0  # cationic antibacterial enzyme
    assert net_charge(sample["sp|P00698|LYSC_CHICK"], 7.0) > 0
    assert net_charge(sample["sp|P0DP23|CALM1_HUMAN"], 7.0) < -20
