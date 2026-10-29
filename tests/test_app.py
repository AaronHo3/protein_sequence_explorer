"""Headless tests for the Streamlit app pages, using Streamlit's AppTest harness."""

from pathlib import Path

import pytest

pytest.importorskip("streamlit", reason="streamlit is in the app extra")

from streamlit.testing.v1 import AppTest  # noqa: E402

from protein_explorer.ui import state  # noqa: E402
from protein_explorer.visualization.sequence_viewer import ResidueRange  # noqa: E402

APP = Path(__file__).resolve().parent.parent / "app" / "streamlit_app.py"
FIXTURES = Path(__file__).resolve().parent / "fixtures"
TIMEOUT = 60
INSULIN = "sp|P01308|INS_HUMAN"


@pytest.fixture
def app() -> AppTest:
    return AppTest.from_file(str(APP), default_timeout=TIMEOUT).run()


def _html(at: AppTest) -> str:
    return " ".join(e.proto.body for e in at.get("html"))


def _open(at: AppTest, page: str) -> AppTest:
    return at.switch_page(f"views/{page}.py").run()


@pytest.mark.parametrize("page", ["import_data", "overview", "protein", "compare", "export"])
def test_every_page_renders_without_errors(app, page):
    _open(app, page)
    assert not app.exception
    assert not app.error


def test_overview_is_default_with_dataset_metrics(app):
    assert app.title[0].value == "Dataset overview"
    html = _html(app)
    assert "Total residues" in html and "1,910" in html


def test_protein_page_shows_selected_protein(app):
    app.session_state[state.PROTEIN] = INSULIN
    _open(app, "protein")
    html = _html(app)
    assert "INS_HUMAN" in html and "11.98" in html and "P01308 on UniProt" in html


def test_residue_selection_is_shared_with_inspector(app):
    app.session_state[state.PROTEIN] = INSULIN
    app.session_state[state.SELECTIONS] = {INSULIN: ResidueRange(3, 3)}
    _open(app, "protein")
    html = _html(app)
    assert "Leucine" in html and "3 of 110" in html


def test_range_selection_offers_region_download(app):
    app.session_state[state.PROTEIN] = INSULIN
    app.session_state[state.SELECTIONS] = {INSULIN: ResidueRange(1, 24)}
    _open(app, "protein")
    assert "Residues 1–24" in _html(app)
    assert any("region" in b.label for b in app.get("download_button"))


def test_search_filter_updates_overview(app):
    app.sidebar.text_input[0].input("homo sapiens").run()
    assert "valid, after filters" in _html(app)
    assert ">5<" in _html(app)


def test_no_matches_shows_guidance(app):
    app.sidebar.text_input[0].input("no such protein").run()
    _open(app, "protein")
    assert not app.exception
    assert any("No proteins match" in i.value for i in app.info)


def test_ph_setting_changes_net_charge_label(app):
    next(s for s in app.sidebar.slider if s.label == "pH for net charge").set_value(5.0).run()
    _open(app, "protein")
    assert "Net charge at pH 5" in _html(app)


def test_colour_mode_persists_across_pages(app):
    _open(app, "protein")
    next(s for s in app.selectbox if s.label == "Colour residues by").set_value("charge").run()
    _open(app, "overview")
    _open(app, "protein")
    assert next(s for s in app.selectbox if s.label == "Colour residues by").value == "charge"


def test_invalid_dataset_marks_records_on_import_page(app):
    data = (FIXTURES / "mixed_validity.fasta").read_bytes()
    app.session_state[state.DATASET] = state.Dataset("mixed_validity.fasta", data, False)
    _open(app, "import_data")
    html = _html(app)
    assert "Failed validation" in html
    assert any("need attention" in w.value for w in app.warning)


def test_compare_page_has_all_sections(app):
    _open(app, "compare")
    headers = [h.value for h in app.main.subheader]
    assert headers == ["Property scatter", "Normalized properties", "Composition heatmap", "Comparison table"]


def test_export_page_offers_all_downloads(app):
    _open(app, "export")
    labels = [b.label for b in app.get("download_button")]
    for expected in ["Features", "Composition", "Residue counts", "Summary statistics",
                     "Dataset composition", "Analysis report (Markdown)", "Run information (JSON)"]:
        assert expected in labels
    assert any(label.startswith("All records") for label in labels)
