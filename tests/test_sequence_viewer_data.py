"""Tests for the sequence viewer data model (the payload sent to the browser)."""

import json

import pytest

from protein_explorer.visualization.color_schemes import color_scheme
from protein_explorer.visualization.sequence_viewer import (
    ResidueAnnotation,
    ResidueRange,
    ViewerOptions,
    build_viewer_data,
    hydropathy_track,
)
from protein_explorer.visualization.theme import DARK, LIGHT

SEQ = "MALWMRLLPLLALLALWGPDPAAA"  # preproinsulin signal peptide, residues 1-24


def _data(sequence=SEQ, **kwargs):
    return build_viewer_data(sequence, color_scheme("chemical_class", LIGHT), LIGHT, **kwargs)


def test_residue_range_is_one_based_inclusive():
    selection = ResidueRange(5, 18)
    assert selection.length == 14
    assert selection.label() == "5–18"
    assert ResidueRange(7, 7).label() == "7"


@pytest.mark.parametrize(("start", "end"), [(0, 3), (5, 4), (-1, 2)])
def test_residue_range_rejects_invalid(start, end):
    with pytest.raises(ValueError):
        ResidueRange(start, end)


def test_range_from_positions():
    assert ResidueRange.from_positions([12, 9, 15]) == ResidueRange(9, 15)
    assert ResidueRange.from_positions([]) is None


def test_payload_is_json_serialisable_and_sends_sequence_once():
    data = _data()
    json.dumps(data)
    assert data["sequence"] == SEQ
    assert set(data["lookup"]) == set(SEQ)


def test_lookup_carries_codes_class_and_hydropathy():
    leucine = _data()["lookup"]["L"]
    assert leucine["three"] == "Leu"
    assert leucine["chemicalClass"] == "Nonpolar"
    assert leucine["hydropathy"] == 3.8


def test_terminal_stop_is_not_displayed():
    assert _data("MKT*")["sequence"] == "MKT"


def test_selection_and_line_length_passed_through():
    data = _data(options=ViewerOptions(line_length=60), selection=ResidueRange(3, 9))
    assert data["selection"] == [3, 9]
    assert data["lineLength"] == 60
    assert _data()["lineLength"] is None  # 0 means fit to width


def test_hydropathy_track_aligned_to_residue_positions():
    track = hydropathy_track(SEQ, 9, LIGHT)
    assert len(track.values) == len(SEQ)
    assert track.values[:4] == (None,) * 4  # no full window before position 5
    assert track.values[4] is not None
    assert track.values[-4:] == (None,) * 4


def test_track_shorter_than_window_is_all_none():
    assert set(hydropathy_track("MKT", 9, LIGHT).values) == {None}


def test_tracks_hidden_when_disabled():
    track = hydropathy_track(SEQ, 9, LIGHT)
    assert len(_data(options=ViewerOptions(tracks=(track,)))["tracks"]) == 1
    assert _data(options=ViewerOptions(tracks=(track,), show_tracks=False))["tracks"] == []


def test_annotations_extension_point():
    data = _data(options=ViewerOptions(annotations=(ResidueAnnotation(1, 24, "Signal peptide", "#000"),)))
    assert data["annotations"] == [{"start": 1, "end": 24, "label": "Signal peptide", "color": "#000"}]


def test_positions_support_gapped_alignment_rows():
    data = build_viewer_data("MK-T", color_scheme("chemical_class", LIGHT), LIGHT,
                             positions=[1, 2, None, 3])
    assert data["positions"] == [1, 2, None, 3]
    assert data["lookup"]["-"]["name"] == "Gap"


def test_positions_length_must_match():
    with pytest.raises(ValueError):
        build_viewer_data("MKT", color_scheme("chemical_class", LIGHT), LIGHT, positions=[1, 2])


def test_scheme_and_palette_change_fills():
    light = _data()["lookup"]["D"]["fill"]
    dark = build_viewer_data(SEQ, color_scheme("chemical_class", DARK), DARK)["lookup"]["D"]["fill"]
    charge = build_viewer_data(SEQ, color_scheme("charge", LIGHT), LIGHT)["lookup"]["D"]["fill"]
    assert len({light, dark, charge}) == 3


def test_nonstandard_residue_gets_neutral_fill():
    data = _data("MXK")
    assert data["lookup"]["X"]["fill"] == LIGHT.neutral_fill
    assert data["lookup"]["X"]["standard"] is False
