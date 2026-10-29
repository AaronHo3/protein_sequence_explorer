# Sidebar shared by every page: dataset status, protein picker, pH, and filters.

from html import escape

import pandas as pd
import streamlit as st

from protein_explorer.ui import state
from protein_explorer.ui.protein_selector import protein_selector

def _plural(n: int, noun: str) -> str:
    return f"{n:,} {noun}" + ("" if n == 1 else "s")

def _dataset_status(result: state.Analysis) -> None:
    ds = state.dataset()
    features = result.features
    n_valid = int(features["is_valid"].astype(bool).sum())
    source = "Example dataset" if ds.is_example else "Uploaded file"
    st.html(
        f'<div class="pe-card"><div class="pe-card-label">{source}</div>'
        f'<div style="font-weight:600;margin:2px 0 4px;word-break:break-word">{escape(ds.name)}</div>'
        f'<div class="pe-card-detail">{_plural(len(features), "protein")}, {n_valid} valid</div></div>'
    )

def _length_filter(features: pd.DataFrame) -> None:
    lengths = features["length"].dropna().astype(int)
    low, high = int(lengths.min()), int(lengths.max())
    if low >= high:
        st.session_state.pop(state.LENGTH_RANGE, None)
        return
    current = st.session_state.get(state.LENGTH_RANGE)
    if not current or current[0] < low or current[1] > high:
        st.session_state[state.LENGTH_RANGE] = (low, high)
    st.slider("Length (residues)", low, high, key=state.LENGTH_RANGE)

def render_sidebar(result: state.Analysis) -> pd.DataFrame:
    """Draw the sidebar and return the filtered feature table."""
    with st.sidebar:
        _dataset_status(result)
        shown = state.filtered_features(result.features)
        protein_selector(shown)
        st.caption(f"{len(shown)} of {_plural(len(result.features), 'protein')} match the filters.")

        st.subheader("Filter proteins")
        st.text_input("Search", key=state.QUERY, placeholder="Name, accession, organism")
        st.toggle("Valid sequences only", key=state.VALID_ONLY,
                  help="Hide sequences that failed validation (they stay in exports).")
        _length_filter(result.features)

        st.subheader("Settings")
        st.slider("pH for net charge", 0.0, 14.0, step=0.1, key=state.PH,
                  help="Used for net charge and side-chain charges. pI does not depend on it.")
    return shown
