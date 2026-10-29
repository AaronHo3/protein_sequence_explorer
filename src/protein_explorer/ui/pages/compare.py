# Comparison workspace: property scatter, normalized property heatmap, composition
# heatmap, and a sortable table, all for one chosen set of proteins.

from functools import partial

import pandas as pd
import streamlit as st

from protein_explorer.residues import RESIDUE_ORDER
from protein_explorer.ui import state
from protein_explorer.ui.chart_controls import show_chart
from protein_explorer.ui.styles import page_header
from protein_explorer.visualization.comparison import (
    comparison_features,
    comparison_scatter,
    feature_label,
    feature_name,
    normalized_feature_heatmap,
)
from protein_explorer.visualization.heatmaps import ROW_ORDERS, composition_heatmap

def _on_scatter_select(key: str, options: list[str]) -> None:
    """One point selects that protein; several (box) become the comparison set."""
    event = st.session_state.get(key)
    points = (event or {}).get("selection", {}).get("points", [])
    ids = [p["customdata"][0] for p in points if p.get("customdata")]
    ids = [i for i in dict.fromkeys(ids) if i in options]
    if len(ids) == 1:
        state.select_protein(ids[0])
        st.toast(f"Selected {state.entry_label(ids[0])}. The Protein page now shows it.",
                 icon=":material/biotech:")
    elif len(ids) > 1:
        st.session_state[state.COMPARE] = ids
        st.toast(f"Comparing {len(ids)} proteins from the selection.", icon=":material/compare_arrows:")
    if ids:
        state.bump_chart_version()

def _chosen_proteins(options: list[str], labels: pd.Series) -> list[str]:
    current = st.session_state.get(state.COMPARE)
    st.session_state[state.COMPARE] = [i for i in (current or options) if i in options] or options
    return st.multiselect("Proteins to compare", options, key=state.COMPARE,
                          format_func=lambda i: labels[i],
                          help="Defaults to every valid protein that matches the sidebar filters.")

def _scatter_section(data: pd.DataFrame, labels: pd.Series, options: list[str]) -> None:
    st.subheader("Property scatter")
    columns = comparison_features(data)
    c1, c2, c3, c4, c5 = st.columns([1.3, 1.3, 1.4, 0.8, 1], vertical_alignment="bottom")
    x = c1.selectbox("Horizontal axis", columns, index=columns.index("isoelectric_point"),
                     format_func=feature_label, key="pe_cmp_x")
    y = c2.selectbox("Vertical axis", columns, index=columns.index("gravy"),
                     format_func=feature_label, key="pe_cmp_y")
    query = c3.text_input("Highlight", placeholder="Search proteins", key="pe_cmp_highlight")
    show_labels = c4.toggle("Labels", value=True, key="pe_cmp_labels")
    by_organism = c5.toggle("By organism", value=False, key="pe_cmp_group",
                            help="Colour by organism; up to 3 groups, the rest as Other.")
    needle = query.strip().lower()
    highlighted = [i for i in data.index if needle and (
        needle in i.lower() or needle in str(data.loc[i, "description"]).lower())] if needle else []
    if needle and not highlighted:
        st.caption(f"No proteins match '{query}'.")
    key = state.chart_key("pe_cmp_scatter")
    fig = comparison_scatter(data, x, y, state.palette(), labels=labels, highlighted=highlighted,
                             selected=st.session_state.get(state.PROTEIN), show_labels=show_labels,
                             group_by="organism" if by_organism else None)
    show_chart(fig, f"{x}_vs_{y}", key=key, on_select=partial(_on_scatter_select, key, options))
    st.caption("Click a point to open that protein; drag a box to compare only those proteins. "
               "The ring marks the protein selected in the sidebar.")

def _normalized_section(data: pd.DataFrame, labels: pd.Series) -> None:
    st.subheader("Normalized properties")
    columns = comparison_features(data)
    chosen = st.multiselect("Properties", columns, default=columns, format_func=feature_name,
                            key="pe_cmp_props")
    if len(data) < 2 or not chosen:
        st.info("Choose at least two proteins and one property.", icon=":material/info:")
        return
    show_chart(normalized_feature_heatmap(data, chosen, state.palette(), labels=labels),
               "normalized_properties")
    st.caption("Each column is z-scored across the proteins shown: (value − mean) / standard "
               "deviation. Colour compares a protein with the others in this set, not with "
               "proteins in general, and never compares properties with different units "
               "directly. Hover for the raw value.")

def _composition_section(result: state.Analysis, ids: list[str], labels: pd.Series) -> None:
    st.subheader("Composition heatmap")
    orders = {**ROW_ORDERS, **{r: f"By {r} frequency" for r in RESIDUE_ORDER}}
    order = st.selectbox("Row order", list(orders), format_func=orders.__getitem__,
                         key="pe_cmp_order",
                         help="Clustering uses average linkage on Euclidean distance between "
                              "composition vectors.")
    show_chart(composition_heatmap(result.composition.loc[ids], state.palette(), order=order,
                                   labels=labels), "composition_heatmap")

def _table_section(data: pd.DataFrame, labels: pd.Series) -> None:
    st.subheader("Comparison table")
    columns = comparison_features(data)
    table = data[["protein_name", "organism", *columns]].copy()
    table.insert(0, "entry", labels.reindex(table.index))
    query = st.text_input("Filter rows", placeholder="Type to filter", key="pe_cmp_table_query")
    if query.strip():
        text = table.astype(str).agg(" ".join, axis=1).str.lower()
        table = table[text.str.contains(query.strip().lower(), regex=False)]
    config = {c: st.column_config.NumberColumn(feature_name(c), format="%.3f") for c in columns}
    config.update({
        "length": st.column_config.NumberColumn("Length", format="%d"),
        "molecular_weight": st.column_config.NumberColumn("Molecular weight (Da)", format="%.1f"),
        "entry": st.column_config.TextColumn("Entry"),
        "protein_name": st.column_config.TextColumn("Protein"),
        "organism": st.column_config.TextColumn("Organism"),
    })
    st.dataframe(table, width="stretch", column_config=config, hide_index=True)
    st.download_button("Download table as CSV", table.to_csv(), "comparison.csv", "text/csv",
                       icon=":material/download:")

def render() -> None:
    result = state.analysis()
    shown = state.filtered_features(result.features)
    labels = result.labels()
    page_header("Compare proteins",
                "Investigate how the selected proteins differ in their biochemical properties "
                "and composition.")
    options = list(state.valid(shown).index)
    if not options:
        st.info("No valid proteins match the current filters.", icon=":material/filter_alt_off:")
        return
    ids = _chosen_proteins(options, labels)
    if not ids:
        st.info("Choose proteins to compare.", icon=":material/info:")
        return
    data = result.features.loc[ids]
    _scatter_section(data, labels, options)
    _normalized_section(data, labels)
    _composition_section(result, ids, labels)
    _table_section(data, labels)
