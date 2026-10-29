"""Tests for Streamlit-side helpers that can run without a browser."""

import pandas as pd
import pytest

pytest.importorskip("streamlit", reason="streamlit is in the app extra")

from protein_explorer.dataset import build_feature_table  # noqa: E402
from protein_explorer.fasta import FastaRecord  # noqa: E402
from protein_explorer.ui.metric_cards import Metric, metric_cards_html  # noqa: E402
from protein_explorer.ui.pages.import_page import record_status  # noqa: E402
from protein_explorer.ui.residue_inspector import (  # noqa: E402
    range_html,
    range_summary,
    single_residue_html,
)
from protein_explorer.ui.sequence_viewer import parse_selection_event  # noqa: E402
from protein_explorer.visualization.color_schemes import color_scheme  # noqa: E402
from protein_explorer.visualization.sequence_viewer import ResidueRange  # noqa: E402
from protein_explorer.visualization.theme import LIGHT  # noqa: E402

SIGNAL = "MALWMRLLPLLALLALWGPDPAAA"


def test_metric_cards_escape_text():
    html = metric_cards_html([Metric("<b>x</b>", "1 & 2", "kDa", detail="d", help='say "hi"')])
    assert "<b>x</b>" not in html and "&lt;b&gt;" in html
    assert "1 &amp; 2" in html and "&quot;hi&quot;" in html


@pytest.mark.parametrize(("event", "expected"), [
    ({"start": 3, "end": 9}, (True, ResidueRange(3, 9))),
    ({"clear": True}, (True, None)),
    (None, (False, None)),
    ({"start": 0, "end": 2}, (False, None)),
    ({"start": "x"}, (False, None)),
])
def test_parse_selection_event(event, expected):
    assert parse_selection_event(event) == expected


def test_single_residue_panel():
    html = single_residue_html(SIGNAL, 3, color_scheme("chemical_class", LIGHT), LIGHT,
                               window_mean=None, window=9, ph=7.0)
    assert "Leucine" in html and "Leu, L" in html and "3 of 24" in html
    assert "+3.8" in html and "no full window here" in html


def test_ionisable_residue_panel_reports_charge():
    html = single_residue_html("MKD", 3, color_scheme("charge", LIGHT), LIGHT,
                               window_mean=0.5, window=3, ph=7.0)
    assert "3.9 (EMBOSS)" in html and "−1.00" in html


def test_range_summary_values():
    rows = dict(range_summary("DDKK", 7.0))
    assert rows["Acidic / basic"] == "2 / 2"
    assert rows["Length"] == "4 residues"
    assert "Residues 2–5" in range_html("MDDKKA", ResidueRange(2, 5), 7.0)


def test_record_status():
    table = build_feature_table([
        FastaRecord("ok", "", "MKTAYIAKQR"),
        FastaRecord("amb", "", "MKTXAYIAKQR"),
        FastaRecord("bad", "", "MK-TAY"),
    ])
    assert record_status(table).tolist() == ["Valid", "Check", "Failed"]
    assert isinstance(record_status(table), pd.Series)
