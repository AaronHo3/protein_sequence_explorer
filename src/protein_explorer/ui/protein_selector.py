# Searchable protein picker (type to filter), shared across pages via session state.

import pandas as pd
import streamlit as st

from protein_explorer.ui import state

def protein_option_label(features: pd.DataFrame, identifier: str) -> str:
    name = features.loc[identifier, "protein_name"] if "protein_name" in features else None
    label = state.entry_label(identifier)
    return f"{label}, {name}" if isinstance(name, str) and name and name != label else label

def protein_selector(shown: pd.DataFrame, container=st.sidebar) -> str | None:
    """Select one protein from the shown (filtered) table. Keeps the current choice when it
    is still shown; otherwise falls back to the first protein."""
    if shown.empty:
        container.caption("No proteins match the filters.")
        return None
    current = state.selected_protein(shown)
    st.session_state[state.PROTEIN] = current
    return container.selectbox(
        "Protein", list(shown.index), key=state.PROTEIN,
        format_func=lambda i: protein_option_label(shown, i),
        help="Type to search by name or accession.",
    )
