# Import: upload or load the example dataset, then review validation results per record.

import pandas as pd
import streamlit as st

from protein_explorer.ui import state
from protein_explorer.ui.metric_cards import Metric, metric_cards
from protein_explorer.ui.styles import page_header

FILE_TYPES = ["fasta", "fa", "faa", "fas", "txt"]

def record_status(features: pd.DataFrame) -> pd.Series:
    """Failed (invalid), Check (valid but with warnings such as X, U, or a stop), or Valid."""
    valid = features["is_valid"].astype(bool)
    has_warnings = features["warnings"].fillna("") != ""
    return pd.Series(
        ["Failed" if not v else "Check" if w else "Valid" for v, w in zip(valid, has_warnings, strict=True)],
        index=features.index,
    )

def _reset_workspace() -> None:
    """Clear filters, selections, and settings; keep the loaded dataset."""
    for key in (*state.PERSISTENT, state.SELECTIONS):
        st.session_state.pop(key, None)
    for key, value in state.DEFAULTS.items():
        if key not in (state.UPLOADER_VERSION, state.CHART_VERSION):
            st.session_state[key] = value
    st.session_state[state.SELECTIONS] = {}

def _upload_panel() -> None:
    st.file_uploader(
        "Drag and drop a protein FASTA file, or browse", type=FILE_TYPES, key=state.uploader_key(),
        on_change=state.handle_upload,
        help=f"Drag and drop or browse. Up to {state.MAX_UPLOAD_MB} MB.",
    )
    error = st.session_state.get(state.UPLOAD_ERROR)
    if error:
        st.error(f"{error} The previous dataset is still loaded.", icon=":material/error:")
    left, right = st.columns(2)
    left.button("Load example dataset", on_click=state.load_example, width="stretch",
                icon=":material/science:", disabled=state.dataset().is_example,
                help="8 reviewed UniProt proteins (data/sample.fasta).")
    right.button("Reset workspace", on_click=_reset_workspace, width="stretch",
                 icon=":material/restart_alt:",
                 help="Clear filters, selections, and settings. Keeps the loaded dataset.")

def _format_panel() -> None:
    with st.container(border=True):
        st.markdown("**What the importer accepts**")
        st.markdown(
            "- Protein FASTA: a `>` header line, then one or more sequence lines\n"
            "- One or many records per file; line breaks and spaces inside sequences are ignored\n"
            "- File types: .fasta, .fa, .faa, .fas, .txt\n"
            "- Sequences are never modified. Invalid characters, DNA/RNA, and ambiguous "
            "residues are reported, not removed"
        )

def render() -> None:
    page_header("Import sequences",
                "Load protein sequences in FASTA format. Every record is validated on import, "
                "and problems are reported rather than silently fixed.")
    left, right = st.columns([3, 2], gap="large")
    with left:
        _upload_panel()
    with right:
        _format_panel()

    result = state.analysis()
    features = result.features
    status = record_status(features)
    valid = features[features["is_valid"].astype(bool)]
    st.subheader("Current dataset")
    metric_cards([
        Metric("Records", str(len(features))),
        Metric("Valid", str(int((status == "Valid").sum()))),
        Metric("Valid with warnings", str(int((status == "Check").sum())),
               help="Contains ambiguous (B/Z/J/X) or nonstandard (U/O) residues, lowercase, or a terminal stop."),
        Metric("Failed validation", str(int((status == "Failed").sum())),
               help="Invalid characters or a likely DNA/RNA sequence. Kept in exports, excluded from charts."),
        Metric("Total residues", f"{int(valid['length'].sum()):,}", detail="valid sequences"),
    ])
    if (status != "Valid").any():
        st.warning("Some records need attention. Warnings are listed in the table below.",
                   icon=":material/warning:")
    table = pd.DataFrame({
        "Status": status,
        "Protein": features["protein_name"],
        "Organism": features["organism"],
        "Length": features["length"],
        "Warnings": features["warnings"].fillna(""),
    })
    st.dataframe(
        table, width="stretch",
        column_config={
            "Status": st.column_config.TextColumn(width="small"),
            "Length": st.column_config.NumberColumn(format="%d aa", width="small"),
            "Warnings": st.column_config.TextColumn(width="large"),
        },
    )
    st.page_link("views/overview.py", label="Continue to the dataset overview",
                 icon=":material/arrow_forward:")
