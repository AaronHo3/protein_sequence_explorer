"""Tests for residue colouring schemes."""

import pytest

from protein_explorer.features import (
    ACIDIC_RESIDUES,
    AMINO_ACIDS,
    AROMATIC_RESIDUES,
    BASIC_RESIDUES,
    HYDROPHOBIC_RESIDUES,
)
from protein_explorer.visualization.color_schemes import (
    NONSTANDARD_LABEL,
    SCHEME_KEYS,
    color_scheme,
    scheme_labels,
)
from protein_explorer.visualization.theme import DARK, LIGHT


@pytest.mark.parametrize("palette", [LIGHT, DARK], ids=["light", "dark"])
@pytest.mark.parametrize("key", SCHEME_KEYS)
def test_every_scheme_partitions_the_twenty_residues(key, palette):
    scheme = color_scheme(key, palette)
    members = "".join(c.residues for c in scheme.categories)
    assert sorted(members) == list(AMINO_ACIDS)
    assert set(scheme.color_map()) == set(AMINO_ACIDS)


def test_scheme_menu_has_seven_modes_in_order():
    assert list(scheme_labels()) == list(SCHEME_KEYS)
    assert len(SCHEME_KEYS) == 7


def _members(key: str, label_prefix: str) -> set[str]:
    scheme = color_scheme(key, LIGHT)
    return set(next(c.residues for c in scheme.categories if c.label.startswith(label_prefix)))


def test_schemes_match_feature_definitions():
    # Colouring by a property highlights exactly the residues the feature counts
    assert _members("hydropathy", "Hydrophobic") == HYDROPHOBIC_RESIDUES
    assert _members("aromatic", "Aromatic") == AROMATIC_RESIDUES
    assert _members("acid_base", "Acidic") == ACIDIC_RESIDUES
    assert _members("acid_base", "Basic") == BASIC_RESIDUES


def test_charge_scheme_separates_partly_charged_histidine():
    assert _members("charge", "Positive") == {"K", "R"}
    assert _members("charge", "Partly positive") == {"H"}
    assert _members("charge", "Negative") == {"D", "E"}


def test_identity_scheme_gives_each_residue_its_own_colour():
    colors = color_scheme("identity", LIGHT).color_map()
    assert len(set(colors.values())) == 20


def test_category_lookup_handles_nonstandard_codes():
    scheme = color_scheme("chemical_class", LIGHT)
    assert scheme.category_of("l") == "Nonpolar"
    assert scheme.category_of("X") == NONSTANDARD_LABEL


def test_unknown_scheme_raises():
    with pytest.raises(KeyError):
        color_scheme("rainbow", LIGHT)
