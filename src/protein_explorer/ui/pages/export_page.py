# Export: tables, sequences, reproducibility records, and figures.

import json

import streamlit as st

from protein_explorer.dataset import summarize_features
from protein_explorer.export import build_analysis_report, build_run_info_from_bytes, records_to_fasta
from protein_explorer.ui import state
from protein_explorer.ui.styles import page_header
from protein_explorer.visualization.export import static_figure_bundle

BUNDLE = "pe_figure_bundle"

def _download(column, label: str, data, filename: str, mime: str, help: str) -> None:
    column.download_button(label, data, filename, mime, width="stretch", help=help,
                           icon=":material/download:")

def render() -> None:
    result = state.analysis()
    dataset = state.dataset()
    shown = state.filtered_features(result.features)
    stem = "example" if dataset.is_example else dataset.name.rsplit(".", 1)[0]
    ph = float(st.session_state[state.PH])
    run_info = build_run_info_from_bytes(dataset.name, dataset.data, result.features, ph)
    composition = state.dataset_composition_table(result, result.features)
    page_header("Export results",
                "Download the calculated tables, sequences, and figures. Tables include every "
                "record in the file, with failed sequences marked rather than dropped.")

    with st.container(border=True):
        st.markdown("**Tables (CSV)**")
        cols = st.columns(4)
        _download(cols[0], "Features", result.features.to_csv(), f"{stem}_features.csv", "text/csv",
                  "One row per protein: metadata, validation, and every feature.")
        _download(cols[1], "Composition", result.composition.to_csv(), f"{stem}_composition.csv",
                  "text/csv", "Amino acid frequencies (fractions), one row per protein.")
        _download(cols[2], "Residue counts", result.counts.to_csv(), f"{stem}_counts.csv",
                  "text/csv", "Residue counts, one row per protein.")
        _download(cols[3], "Summary statistics", summarize_features(result.features).to_csv(),
                  f"{stem}_summary.csv", "text/csv", "Mean, SD, and quartiles across valid proteins.")
        _download(st, "Dataset composition", composition.to_csv(), f"{stem}_dataset_composition.csv",
                  "text/csv", "Residue-weighted and protein-weighted dataset frequencies.")

    with st.container(border=True):
        st.markdown("**Sequences (FASTA)**")
        records = list(result.records)
        valid_ids = set(state.valid(result.features).index)
        cols = st.columns(3)
        _download(cols[0], f"All records ({len(records)})", records_to_fasta(records),
                  f"{stem}.fasta", "text/x-fasta", "Identifiers, descriptions, and residues unchanged.")
        valid_records = [r for r in records if r.identifier in valid_ids]
        _download(cols[1], f"Valid only ({len(valid_records)})", records_to_fasta(valid_records),
                  f"{stem}_valid.fasta", "text/x-fasta", "Records that passed validation.")
        shown_records = [r for r in records if r.identifier in set(shown.index)]
        if shown_records:
            _download(cols[2], f"Matching filters ({len(shown_records)})",
                      records_to_fasta(shown_records), f"{stem}_filtered.fasta", "text/x-fasta",
                      "Records matching the sidebar filters.")

    with st.container(border=True):
        st.markdown("**Reproducibility**")
        cols = st.columns(2)
        _download(cols[0], "Analysis report (Markdown)",
                  build_analysis_report(result.features, run_info, composition),
                  f"{stem}_report.md", "text/markdown",
                  "Provenance, per-protein features, summary statistics, and assumptions.")
        _download(cols[1], "Run information (JSON)", json.dumps(run_info, indent=2),
                  f"{stem}_run_info.json", "application/json",
                  "Input checksum, pH, pKa set, and software versions.")

    with st.container(border=True):
        st.markdown("**Figures**")
        st.caption("Every interactive chart has a camera button in its toolbar that saves a "
                   "PNG at 3× resolution. The bundle below contains the static 300 dpi figures "
                   "for every valid protein.")
        window = int(st.session_state[state.WINDOW])
        if st.button("Prepare figure bundle", icon=":material/photo_library:"):
            with st.spinner("Drawing figures…"):
                st.session_state[BUNDLE] = (dataset.data, window, static_figure_bundle(records, window))
        bundle = st.session_state.get(BUNDLE)
        if bundle and bundle[0] == dataset.data and bundle[1] == window:
            st.download_button("Download figures (ZIP)", bundle[2], f"{stem}_figures.zip",
                               "application/zip", icon=":material/download:")
