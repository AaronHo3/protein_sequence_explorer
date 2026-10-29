# Shared app state: the loaded dataset, cached analysis, filters, and selections. Every
# page reads from here so a selection made on one page is reflected on the others.

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import streamlit as st

from protein_explorer.dataset import (
    build_composition_table,
    build_count_table,
    build_feature_table,
    dataset_composition,
    filter_features,
)
from protein_explorer.fasta import FastaFormatError, FastaRecord, parse_fasta_lines
from protein_explorer.plots import entry_label
from protein_explorer.visualization.sequence_viewer import ResidueRange
from protein_explorer.visualization.theme import ChartPalette, palette_for

SAMPLE_FASTA = Path(__file__).resolve().parents[3] / "data" / "sample.fasta"
SAMPLE_LABEL = "Example: 8 reviewed UniProt proteins"
MAX_UPLOAD_MB = 50

# Session keys. Widget-backed keys are listed in PERSISTENT so their values survive
# navigating to pages where the widget is not drawn.
DATASET = "pe_dataset"
UPLOAD_ERROR = "pe_upload_error"
UPLOADER_VERSION = "pe_uploader_version"
PROTEIN = "pe_protein"
SELECTIONS = "pe_selections"
COMPARE = "pe_compare"
PH = "pe_ph"
WINDOW = "pe_window"
COLOR_MODE = "pe_color_mode"
LINE_LENGTH = "pe_line_length"
SHOW_TRACK = "pe_show_track"
QUERY = "pe_query"
VALID_ONLY = "pe_valid_only"
LENGTH_RANGE = "pe_length_range"
CHART_VERSION = "pe_chart_version"

DEFAULTS: dict[str, object] = {
    PH: 7.0, WINDOW: 9, COLOR_MODE: "chemical_class", LINE_LENGTH: 0, SHOW_TRACK: True,
    QUERY: "", VALID_ONLY: True, UPLOADER_VERSION: 0, CHART_VERSION: 0,
}
PERSISTENT = (PROTEIN, COMPARE, PH, WINDOW, COLOR_MODE, LINE_LENGTH, SHOW_TRACK, QUERY,
              VALID_ONLY, LENGTH_RANGE)

@dataclass(frozen=True)
class Dataset:
    name: str
    data: bytes
    is_example: bool

@dataclass(frozen=True)
class Analysis:
    records: tuple[FastaRecord, ...]
    features: pd.DataFrame
    composition: pd.DataFrame
    counts: pd.DataFrame

    def record(self, identifier: str) -> FastaRecord:
        return next(r for r in self.records if r.identifier == identifier)

    def labels(self) -> pd.Series:
        """Short display label per protein: UniProt entry name if present, else identifier."""
        return pd.Series({i: entry_label(i) for i in self.features.index})

@st.cache_data(show_spinner="Analyzing sequences…", max_entries=8)
def analyze(data: bytes, ph: float) -> Analysis:
    """Parse and compute every table once per (file content, pH)."""
    records = tuple(parse_fasta_lines(data.decode("utf-8").splitlines()))
    return Analysis(
        records=records,
        features=build_feature_table(list(records), ph=ph),
        composition=build_composition_table(list(records)),
        counts=build_count_table(list(records)),
    )

def init_state() -> None:
    """Defaults on first run, then re-save widget values so they persist across pages."""
    for key, value in DEFAULTS.items():
        st.session_state.setdefault(key, value)
    st.session_state.setdefault(SELECTIONS, {})
    if DATASET not in st.session_state:
        st.session_state[DATASET] = Dataset(SAMPLE_LABEL, SAMPLE_FASTA.read_bytes(), True)
    for key in PERSISTENT:
        if key in st.session_state:
            st.session_state[key] = st.session_state[key]

def palette() -> ChartPalette:
    theme = getattr(st.context, "theme", None)
    return palette_for(getattr(theme, "type", None))

def dataset() -> Dataset:
    return st.session_state[DATASET]

def analysis() -> Analysis:
    return analyze(dataset().data, float(st.session_state[PH]))

def _reset_selection_state() -> None:
    for key in (PROTEIN, COMPARE, LENGTH_RANGE):
        st.session_state.pop(key, None)
    st.session_state[SELECTIONS] = {}
    st.session_state[QUERY] = ""

def set_dataset(name: str, data: bytes, is_example: bool = False) -> str | None:
    """Validate and load a dataset. Returns an error message (and keeps the current dataset)
    if the file can't be used."""
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return f"{name} isn't a plain-text file. Upload a FASTA file (.fasta, .fa, .faa)."
    try:
        parse_fasta_lines(text.splitlines())
    except FastaFormatError as error:
        return f"{name} isn't valid FASTA. {error}"
    st.session_state[DATASET] = Dataset(name, data, is_example)
    _reset_selection_state()
    return None

def load_example() -> None:
    set_dataset(SAMPLE_LABEL, SAMPLE_FASTA.read_bytes(), is_example=True)
    st.session_state[UPLOAD_ERROR] = None
    st.session_state[UPLOADER_VERSION] += 1  # clears the file uploader widget

def uploader_key() -> str:
    return f"pe_upload_{st.session_state[UPLOADER_VERSION]}"

def handle_upload() -> None:
    """File uploader callback: load the file or record why it was rejected."""
    uploaded = st.session_state.get(uploader_key())
    if uploaded is None:
        return
    st.session_state[UPLOAD_ERROR] = set_dataset(uploaded.name, uploaded.getvalue())

def filtered_features(features: pd.DataFrame) -> pd.DataFrame:
    length_range = st.session_state.get(LENGTH_RANGE)
    return filter_features(features, st.session_state[QUERY], bool(st.session_state[VALID_ONLY]),
                           tuple(length_range) if length_range else None)

def valid(features: pd.DataFrame) -> pd.DataFrame:
    return features[features["is_valid"].astype(bool)]

def selected_protein(shown: pd.DataFrame) -> str | None:
    """The current protein, falling back to the first shown one if filtered out."""
    current = st.session_state.get(PROTEIN)
    if current in shown.index:
        return current
    return shown.index[0] if len(shown) else None

def select_protein(identifier: str) -> None:
    st.session_state[PROTEIN] = identifier

def residue_selection(identifier: str) -> ResidueRange | None:
    return st.session_state[SELECTIONS].get(identifier)

def set_residue_selection(identifier: str, selection: ResidueRange | None) -> None:
    selections = dict(st.session_state[SELECTIONS])
    if selection is None:
        selections.pop(identifier, None)
    else:
        selections[identifier] = selection
    st.session_state[SELECTIONS] = selections

def chart_key(name: str) -> str:
    """Versioned chart key. Bumping the version after a chart selection is consumed clears
    the chart's own selection box, so the shared selection is the only one shown."""
    return f"{name}_{st.session_state[CHART_VERSION]}"

def bump_chart_version() -> None:
    st.session_state[CHART_VERSION] += 1

def dataset_composition_table(result: Analysis, shown: pd.DataFrame) -> pd.DataFrame:
    ids = valid(shown).index
    return dataset_composition(result.counts.loc[ids], result.features.loc[ids, "length"])
