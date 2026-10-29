# App shell: page config, shared state, styles, sidebar, and top navigation following the
# workflow Import -> Overview -> Protein -> Compare -> Export.

import streamlit as st

from protein_explorer.ui import state
from protein_explorer.ui.sidebar import render_sidebar
from protein_explorer.ui.styles import inject_styles

PAGES = [
    ("views/import_data.py", "Import", ":material/upload_file:"),
    ("views/overview.py", "Overview", ":material/insights:"),
    ("views/protein.py", "Protein", ":material/biotech:"),
    ("views/compare.py", "Compare", ":material/compare_arrows:"),
    ("views/export.py", "Export", ":material/download:"),
]
DEFAULT_PAGE = "views/overview.py"

def run() -> None:
    st.set_page_config(page_title="Protein Sequence Explorer", page_icon="🧬", layout="wide")
    state.init_state()
    inject_styles(state.palette())
    navigation = st.navigation(
        [st.Page(path, title=title, icon=icon, default=path == DEFAULT_PAGE)
         for path, title, icon in PAGES],
        position="top",
    )
    with st.sidebar:
        st.markdown("### Protein Sequence Explorer")
    render_sidebar(state.analysis())  # datasets are validated in state.set_dataset
    navigation.run()
