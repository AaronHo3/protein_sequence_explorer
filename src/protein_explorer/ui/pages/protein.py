# Protein explorer: one protein's metrics, the interactive sequence viewer, residue
# inspector, hydropathy profile, and composition. Selections are shared between the
# viewer, the hydropathy plot, and the inspector through ui.state.

from functools import partial
from html import escape

import pandas as pd
import streamlit as st

from protein_explorer.export import records_to_fasta
from protein_explorer.features import residue_counts
from protein_explorer.glossary import FEATURE_NOTES
from protein_explorer.ui import state
from protein_explorer.ui.chart_controls import choice, show_chart
from protein_explorer.ui.metric_cards import Metric, metric_cards
from protein_explorer.ui.residue_inspector import residue_inspector
from protein_explorer.ui.sequence_viewer import sequence_viewer
from protein_explorer.visualization.color_schemes import SCHEME_KEYS, color_scheme, scheme_labels
from protein_explorer.visualization.composition import SORT_ORDERS, composition_bar_chart
from protein_explorer.visualization.hydropathy import hydropathy_figure, selection_from_points
from protein_explorer.visualization.sequence_viewer import (
    LINE_LENGTH_OPTIONS,
    ViewerOptions,
    build_viewer_data,
    hydropathy_track,
)

WINDOWS = [5, 7, 9, 11, 13, 15, 17, 19, 21]
MAX_RESIDUE_BARS = 3000  # above this, per-residue bars are unreadable and slow to draw
UNIPROT_URL = "https://www.uniprot.org/uniprotkb/{}/entry"

def _header(record, row: pd.Series) -> None:
    name = row.get("protein_name") or record.description or record.identifier
    parts = [escape(str(name))]
    if isinstance(row.get("organism"), str):
        parts.append(f"<em>{escape(row['organism'])}</em>")
    if isinstance(row.get("gene"), str):
        parts.append(f"gene {escape(row['gene'])}")
    if isinstance(row.get("accession"), str):
        url = UNIPROT_URL.format(escape(row["accession"]))
        parts.append(f'<a href="{url}" target="_blank" rel="noopener">{escape(row["accession"])} on UniProt</a>')
    st.html(
        f'<div class="pe-protein-title"><h1>{escape(state.entry_label(record.identifier))}</h1>'
        f'<span class="pe-protein-sub">{", ".join(parts)}</span></div>'
    )

def _metrics(row: pd.Series, ph: float) -> None:
    charge = next(c for c in row.index if c.startswith("net_charge_ph"))
    metric_cards([
        Metric("Length", f"{row['length']:,}", "aa"),
        Metric("Molecular weight", f"{row['molecular_weight'] / 1000:.2f}", "kDa",
               help="Average residue masses plus one water. Sequence as written; no modifications."),
        Metric("Isoelectric point", f"{row['isoelectric_point']:.2f}", detail="estimated",
               help="pH at which the Henderson–Hasselbalch net charge is zero (EMBOSS pKa values)."),
        Metric(f"Net charge at pH {ph:g}", f"{row[charge]:+.2f}", detail="estimated",
               help="Side chains plus both termini; ignores structural pKa shifts."),
        Metric("GRAVY", f"{row['gravy']:+.3f}", detail="Kyte–Doolittle mean",
               help="Positive is hydrophobic overall, negative is hydrophilic."),
        Metric("Aromatic", f"{row['aromatic_fraction']:.1%}", detail="F, W, Y"),
        Metric("Hydrophobic", f"{row['hydrophobic_fraction']:.1%}", detail="A, V, I, L, M, F, W"),
        Metric("Charged", f"{row['charged_fraction']:.1%}", detail="D, E, K, R"),
    ])

def _on_hydropathy_select(identifier: str, key: str) -> None:
    event = st.session_state.get(key)
    points = (event or {}).get("selection", {}).get("points", [])
    selection = selection_from_points(points)
    if selection is not None:
        state.set_residue_selection(identifier, selection)
        state.bump_chart_version()

def _viewer_section(record, palette, window: int, selection) -> None:
    st.subheader("Sequence")
    labels = scheme_labels()
    c1, c2, c3 = st.columns([1.6, 1.2, 1.4], vertical_alignment="bottom")
    mode = c1.selectbox("Colour residues by", list(SCHEME_KEYS), format_func=labels.__getitem__,
                        key=state.COLOR_MODE)
    line = c2.selectbox("Line length", list(LINE_LENGTH_OPTIONS), key=state.LINE_LENGTH,
                        format_func=lambda n: "Fit to width" if not n else f"{n} residues")
    show_track = c3.toggle("Hydropathy strip", key=state.SHOW_TRACK,
                           help=f"Window-mean hydropathy (window {window}) under each residue.")
    scheme = color_scheme(mode, palette)
    st.caption(scheme.description)
    track = hydropathy_track(record.sequence, window, palette)
    data = build_viewer_data(record.sequence, scheme, palette,
                             ViewerOptions(line_length=line, show_tracks=show_track, tracks=(track,)),
                             selection=selection)
    viewer, inspector = st.columns([2.5, 1], gap="large")
    with viewer:
        sequence_viewer(data, key=f"pe_viewer_{record.identifier}",
                        on_select=partial(state.set_residue_selection, record.identifier))
        st.download_button("Download sequence as FASTA", records_to_fasta([record]),
                           file_name=f"{state.entry_label(record.identifier)}.fasta",
                           mime="text/x-fasta", icon=":material/download:")
    with inspector:
        st.markdown("**Residue inspector**")
        residue_inspector(record.sequence, selection, scheme, palette, window=window,
                          window_means=list(track.values), ph=float(st.session_state[state.PH]),
                          record_id=record.identifier, description=record.description)

def _hydropathy_section(record, palette, window: int, selection) -> None:
    st.subheader("Hydropathy profile")
    c1, c2, c3 = st.columns([1.6, 1.2, 1.2], vertical_alignment="bottom")
    c1.select_slider("Window (residues)", WINDOWS, key=state.WINDOW,
                     help="Odd window centred on each residue. 9 shows surface patches; "
                          "19 matches the length of a membrane-spanning helix.")
    too_long = len(record.sequence) > MAX_RESIDUE_BARS
    residues_on = c2.toggle(
        "Residue values", value=True, key="pe_show_residue_values", disabled=too_long,
        help=(f"Hidden for sequences over {MAX_RESIDUE_BARS:,} residues: the bars merge "
              "into noise at that density and slow the chart." if too_long
              else "Kyte–Doolittle value of each individual residue."),
    ) and not too_long
    reference = c3.toggle("1.6 reference", value=False, key="pe_show_reference",
                          help="Kyte & Doolittle's cutoff for window 19. A peak above it is "
                               "not evidence of a transmembrane helix on its own.")
    key = state.chart_key(f"pe_hydro_{record.identifier}")
    show_chart(hydropathy_figure(record.sequence, palette, window, selection=selection,
                                 show_residue_values=residues_on, show_reference=reference),
               "hydropathy_profile", key=key,
               on_select=partial(_on_hydropathy_select, record.identifier, key))
    st.caption(
        f"Bars are each residue's Kyte–Doolittle value; the line is the mean of a "
        f"{window}-residue window centred on each position. The first and last {window // 2} "
        "residues have no full window, so the line starts and ends inside the sequence. "
        "Drag across the plot to select residues; the selection is shared with the sequence "
        "viewer and inspector above."
    )

def _composition_section(record, row, palette, reference) -> None:
    st.subheader("Amino acid composition")
    c1, c2, c3 = st.columns([1.2, 2, 1.4], vertical_alignment="bottom")
    with c1:
        unit = choice("Show", {"percent": "Percent", "count": "Count"}, "pe_comp_unit", "percent")
    with c2:
        order = choice("Order", SORT_ORDERS, "pe_comp_order", "class")
    show_ref = c3.toggle("Dataset reference", value=True, key="pe_comp_ref",
                         help="Protein-weighted mean composition of the valid proteins shown.")
    counts = pd.Series(residue_counts(record.sequence))
    show_chart(composition_bar_chart(counts, int(row["length"]), palette, as_percent=unit == "percent",
                                     order=order, reference=reference if show_ref else None),
               "composition")

def render() -> None:
    result = state.analysis()
    palette = state.palette()
    shown = state.filtered_features(result.features)
    identifier = state.selected_protein(shown)
    if identifier is None:
        st.title("Protein")
        st.info("No proteins match the current filters. Adjust them in the sidebar.",
                icon=":material/filter_alt_off:")
        return
    record = result.record(identifier)
    row = result.features.loc[identifier]
    window = int(st.session_state[state.WINDOW])
    selection = state.residue_selection(identifier)
    if selection is not None and selection.end > row["length"]:
        selection = None

    _header(record, row)
    if row["warnings"]:
        (st.error if not row["is_valid"] else st.warning)(row["warnings"], icon=":material/warning:")
    if row["is_valid"]:
        _metrics(row, float(st.session_state[state.PH]))
    _viewer_section(record, palette, window, selection)
    if not row["is_valid"]:
        st.info("Features and charts are not calculated for sequences that fail validation.",
                icon=":material/info:")
        return
    _hydropathy_section(record, palette, window, selection)
    reference = state.dataset_composition_table(result, shown)["protein_weighted"]
    _composition_section(record, row, palette, reference)
    with st.expander("How these values are calculated"):
        for note in FEATURE_NOTES:
            st.markdown(f"**{note.name}.** {note.meaning}\n\n*Calculation:* {note.calculation}\n\n"
                        f"*Limitations:* {note.caveat}")
