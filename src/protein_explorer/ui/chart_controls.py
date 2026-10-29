# Small, consistent control widgets used above charts.

from collections.abc import Mapping

import streamlit as st

from protein_explorer.visualization.theme import plotly_config

def choice(label: str, options: Mapping[str, str], key: str, default: str,
           help: str | None = None) -> str:
    """Segmented control over {value: label}. Never returns None: clicking the active
    option again keeps the current choice."""
    st.session_state.setdefault(key, default)
    if st.session_state[key] not in options:
        st.session_state[key] = default
    value = st.segmented_control(label, list(options), format_func=options.__getitem__,
                                 key=key, help=help)
    return value if value is not None else default

def show_chart(fig, filename: str, key: str | None = None, on_select=None,
               selection_mode=("points", "box")):
    """Plotly chart with the shared toolbar config (PNG export at 3x)."""
    kwargs = {"width": "stretch", "theme": None, "config": plotly_config(filename), "key": key}
    if on_select is not None:
        kwargs.update(on_select=on_select, selection_mode=selection_mode)
    return st.plotly_chart(fig, **kwargs)
