# Overview: dataset-level summary of the proteins that pass the sidebar filters.

import streamlit as st

from protein_explorer.ui import state
from protein_explorer.ui.chart_controls import choice, show_chart
from protein_explorer.ui.metric_cards import Metric, metric_cards
from protein_explorer.ui.styles import page_header
from protein_explorer.visualization.composition import SORT_ORDERS, dataset_composition_chart
from protein_explorer.visualization.distributions import (
    default_bins,
    distribution_histogram,
    mw_vs_length,
)

def _summary_cards(valid) -> None:
    lengths, mw, pi = valid["length"].astype(float), valid["molecular_weight"], valid["isoelectric_point"]
    metric_cards([
        Metric("Proteins", str(len(valid)), detail="valid, after filters"),
        Metric("Total residues", f"{int(lengths.sum()):,}"),
        Metric("Mean length", f"{lengths.mean():.0f}", "aa"),
        Metric("Median length", f"{lengths.median():.0f}", "aa"),
        Metric("Mean molecular weight", f"{mw.mean() / 1000:.1f}", "kDa",
               help="Average residue masses; precursor sequences as written, no modifications."),
        Metric("Estimated pI", f"{pi.median():.2f}", detail=f"median, range {pi.min():.2f}–{pi.max():.2f}",
               help="Henderson–Hasselbalch with EMBOSS pKa values."),
    ])

def render() -> None:
    result = state.analysis()
    palette = state.palette()
    shown = state.filtered_features(result.features)
    valid = state.valid(shown)
    labels = result.labels()
    page_header("Dataset overview",
                f"{state.dataset().name}. Summaries use the {len(valid)} valid proteins that "
                "match the sidebar filters.")
    if valid.empty:
        st.info("No valid proteins match the current filters. Adjust them in the sidebar.",
                icon=":material/filter_alt_off:")
        return
    _summary_cards(valid)

    left, right = st.columns(2, gap="large")
    with left:
        st.subheader("Length distribution")
        bins = st.slider("Bins", 1, 40, default_bins(len(valid)), key="pe_length_bins")
        show_chart(distribution_histogram(valid["length"], palette, xlabel="Length (residues)",
                                          bins=bins, labels=labels, value_format=".0f"),
                   "length_distribution")
    with right:
        st.subheader("Estimated isoelectric point")
        bins = st.slider("Bins", 1, 40, default_bins(len(valid)), key="pe_pi_bins")
        show_chart(distribution_histogram(valid["isoelectric_point"], palette,
                                          xlabel="Isoelectric point (pH, estimated)",
                                          bins=bins, labels=labels, value_format=".2f"),
                   "pi_distribution")

    left, right = st.columns(2, gap="large")
    with left:
        st.subheader("Molecular weight vs. length")
        fit = st.toggle("Show linear fit", value=True, key="pe_mw_fit",
                        help="Least squares. The slope is the average mass per residue.")
        show_chart(mw_vs_length(valid, palette, labels=labels,
                                highlight=state.selected_protein(shown), show_fit=fit),
                   "mw_vs_length")
    with right:
        st.subheader("Dataset composition")
        order = choice("Order", SORT_ORDERS, "pe_dataset_order", "class")
        composition = state.dataset_composition_table(result, shown)
        show_chart(dataset_composition_chart(composition, palette, order=order),
                   "dataset_composition")
        st.caption("Residue-weighted pools every residue, so long proteins dominate. "
                   "Protein-weighted averages each protein's own composition, so every "
                   "protein counts equally.")
