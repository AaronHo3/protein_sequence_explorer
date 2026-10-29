"""Validation of our feature calculations against Biopython.

Tolerances are justified per feature:
- molecular_weight: our residue masses are Biopython's rounded to 4 dp, so error is at
  most 0.00005 Da per residue -> relative 1e-6 is generous.
- aromaticity, GRAVY: identical definitions and scale -> exact up to float rounding.
- net charge with Biopython's pKa set: identical formula -> exact up to float rounding.
- pI: both bisect to ~1e-4 from different starting intervals -> 1e-3.
"""

from pathlib import Path

import pytest

pytest.importorskip("Bio", reason="biopython is a dev dependency")

from Bio.SeqUtils.ProtParam import ProteinAnalysis  # noqa: E402

from protein_explorer import features  # noqa: E402
from protein_explorer.crosscheck import (  # noqa: E402
    BIOPYTHON_BASE_PKA,
    biopython_pka_set,
    build_report,
    compare_features,
    decompose_pi,
    edge_case_behaviour,
)
from protein_explorer.fasta import read_fasta  # noqa: E402

SAMPLE_FASTA = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"
MW_REL_TOL = 1e-6
FLOAT_TOL = 1e-9
PI_TOL = 1e-3


@pytest.fixture(scope="module")
def records():
    return read_fasta(SAMPLE_FASTA)


def _ids(records):
    return [r.identifier.split("|")[-1] for r in records]


@pytest.fixture(params=range(8), ids=_ids(read_fasta(SAMPLE_FASTA)))
def sequence(request, records) -> str:
    return records[request.param].sequence


# --- Features with identical definitions -------------------------------------


def test_molecular_weight_matches(sequence):
    expected = ProteinAnalysis(sequence).molecular_weight()
    assert features.molecular_weight(sequence) == pytest.approx(expected, rel=MW_REL_TOL)


def test_aromaticity_matches(sequence):
    expected = ProteinAnalysis(sequence).aromaticity()
    assert features.aromatic_fraction(sequence) == pytest.approx(expected, abs=FLOAT_TOL)


def test_gravy_matches(sequence):
    expected = ProteinAnalysis(sequence).gravy()
    assert features.gravy(sequence) == pytest.approx(expected, abs=FLOAT_TOL)


# --- Charge and pI: reproduce Biopython when given its conventions ------------


@pytest.mark.parametrize("ph", [3.0, 7.0, 7.4, 10.0])
def test_net_charge_matches_with_biopython_pka(sequence, ph):
    expected = ProteinAnalysis(sequence).charge_at_pH(ph)
    ours = features.net_charge(sequence, ph, biopython_pka_set(sequence))
    assert ours == pytest.approx(expected, abs=FLOAT_TOL)


def test_isoelectric_point_matches_with_biopython_pka(sequence):
    expected = ProteinAnalysis(sequence).isoelectric_point()
    ours = features.isoelectric_point(sequence, biopython_pka_set(sequence))
    assert ours == pytest.approx(expected, abs=PI_TOL)


def test_emboss_pi_differs_but_by_less_than_one_unit(sequence):
    # Documented difference from pKa choice, not a bug; a 1-unit gap would suggest one
    difference = features.isoelectric_point(sequence) - ProteinAnalysis(sequence).isoelectric_point()
    assert 0 < abs(difference) < 1.0


# --- Terminal pKa adjustment --------------------------------------------------


def test_terminal_adjustment_uses_n_terminal_residue():
    # Biopython: N-terminal Met -> 7.0, Ala -> 7.59; unlisted residue keeps the base 7.5
    assert biopython_pka_set("MKT").positive["Nterm"] == 7.0
    assert biopython_pka_set("AKT").positive["Nterm"] == 7.59
    assert biopython_pka_set("GKT").positive["Nterm"] == BIOPYTHON_BASE_PKA.positive["Nterm"]


def test_terminal_adjustment_uses_c_terminal_residue_ignoring_stop():
    assert biopython_pka_set("MKD").negative["Cterm"] == 4.55
    assert biopython_pka_set("MKE*").negative["Cterm"] == 4.75


def test_pka_sets_are_read_only():
    with pytest.raises(TypeError):
        biopython_pka_set("MKT").positive["K"] = 1.0  # type: ignore[index]


# --- Decomposition and report -------------------------------------------------


def test_pi_decomposition_explains_whole_difference(records):
    pi = decompose_pi(records)
    total = pi["biopython"] - pi["ours_emboss"]
    explained = pi["pka_set_effect"] + pi["terminal_effect"]
    assert (total - explained).abs().max() < PI_TOL
    assert pi["remaining"].abs().max() < PI_TOL


def test_pka_choice_preserves_pi_ranking(records):
    pi = decompose_pi(records)
    assert (pi["ours_emboss"].rank() == pi["biopython"].rank()).all()


def test_compare_features_covers_every_protein_and_feature(records):
    table = compare_features(records)
    assert len(table) == 8 * 5
    agreeing = table[table["feature"].isin(["molecular_weight", "aromaticity", "gravy"])]
    assert agreeing["relative_difference"].abs().max() < MW_REL_TOL


def test_edge_cases_document_known_differences():
    table = edge_case_behaviour().set_index("case")
    assert table.loc["unknown X", "gravy_biopython"] == "error (KeyError)"
    assert table.loc["unknown X", "gravy_ours"] == "-0.7800"
    assert table.loc["terminal stop *", "aromaticity_ours"] == "0.1000"
    assert table.loc["terminal stop *", "aromaticity_biopython"] == "0.0909"
    # U has a pKa in our set only, so our EMBOSS pI moves but Biopython's does not
    assert table.loc["selenocysteine U", "pi_biopython"] == table.loc["standard", "pi_biopython"]
    assert table.loc["selenocysteine U", "pi_ours_emboss"] != table.loc["standard", "pi_ours_emboss"]


def test_biopython_clamps_pi_to_search_range():
    peptide = "D" * 20
    bio = ProteinAnalysis(peptide)
    assert bio.isoelectric_point() == pytest.approx(4.05, abs=PI_TOL)
    assert bio.charge_at_pH(bio.isoelectric_point()) < -5  # not a true pI
    ours = features.isoelectric_point(peptide, biopython_pka_set(peptide))
    assert features.net_charge(peptide, ours, biopython_pka_set(peptide)) == pytest.approx(
        0, abs=1e-2
    )


def test_report_contains_all_sections(records):
    report = build_report(records, SAMPLE_FASTA)
    for heading in ["Agreement by feature", "Breaking down the pI difference",
                    "Mass conventions", "Residues outside the standard 20", "Full comparison"]:
        assert f"## {heading}" in report
    assert "P01308" in report
