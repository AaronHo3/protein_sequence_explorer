# Residue inspection panel for the current selection: one residue in detail, or summary
# statistics for a range. Values come from residues.py and features.py.

from collections import Counter
from html import escape

import streamlit as st

from protein_explorer.features import (
    ACIDIC_RESIDUES,
    AROMATIC_RESIDUES,
    BASIC_RESIDUES,
    HYDROPHOBIC_RESIDUES,
    KYTE_DOOLITTLE,
)
from protein_explorer.residues import residue_info, side_chain_charge
from protein_explorer.visualization.color_schemes import ColorScheme
from protein_explorer.visualization.sequence_viewer import ResidueRange
from protein_explorer.visualization.theme import ChartPalette, readable_text

def _signed(value: float, digits: int = 2) -> str:
    return f"{value:+.{digits}f}".replace("-", "−")

def _rows(rows: list[tuple[str, str]]) -> str:
    return "<dl>" + "".join(f"<dt>{escape(k)}</dt><dd>{escape(v)}</dd>" for k, v in rows) + "</dl>"

def single_residue_html(
    sequence: str, position: int, scheme: ColorScheme, palette: ChartPalette,
    window_mean: float | None, window: int, ph: float,
) -> str:
    symbol = sequence[position - 1]
    info = residue_info(symbol)
    fill = scheme.color_map().get(info.code, scheme.nonstandard_color)
    rows = [
        ("Position", f"{position} of {len(sequence)}"),
        ("Class", info.chemical_class),
        ("Hydropathy (residue)", _signed(info.hydropathy, 1) if info.hydropathy is not None else "not defined"),
        (f"Window mean (w={window})", _signed(window_mean) if window_mean is not None else "no full window here"),
    ]
    if info.side_chain_pka is not None:
        rows.append(("Side-chain pKa", f"{info.side_chain_pka:g} (EMBOSS)"))
        rows.append((f"Side-chain charge at pH {ph:g}", _signed(side_chain_charge(info.code, ph))))
    else:
        rows.append(("Side-chain charge", "none (not ionisable)"))
    if info.code == "C":
        rows.append(("Note", "Assumes a free thiol; a disulfide-bonded Cys does not ionise"))
    if not info.is_standard:
        rows.append(("Note", "Ambiguous or nonstandard code; excluded from averages"))
    return (
        '<div class="pe-inspector"><div class="pe-inspector-head">'
        f'<div class="pe-inspector-tile" style="background:{fill};color:{readable_text(fill, palette)}">{escape(symbol)}</div>'
        f'<div><div class="pe-inspector-name">{escape(info.name)}</div>'
        f'<div class="pe-inspector-codes">{escape(info.three_letter)}, {escape(info.code)}</div></div></div>'
        f"{_rows(rows)}</div>"
    )

def range_summary(segment: str, ph: float) -> list[tuple[str, str]]:
    counts = Counter(segment.upper())
    kd = [KYTE_DOOLITTLE[r] for r in segment.upper() if r in KYTE_DOOLITTLE]
    def n(group: frozenset[str]) -> int:
        return sum(counts[r] for r in group)
    charge = sum(side_chain_charge(r, ph) for r in segment)
    return [
        ("Length", f"{len(segment)} residues"),
        ("Mean hydropathy", _signed(sum(kd) / len(kd)) if kd else "not defined"),
        ("Hydrophobic", f"{n(HYDROPHOBIC_RESIDUES)} ({n(HYDROPHOBIC_RESIDUES) / len(segment):.0%})"),
        ("Acidic / basic", f"{n(ACIDIC_RESIDUES)} / {n(BASIC_RESIDUES)}"),
        ("Aromatic", str(n(AROMATIC_RESIDUES))),
        (f"Side-chain charge at pH {ph:g}", _signed(charge)),
    ]

def range_html(sequence: str, selection: ResidueRange, ph: float) -> str:
    segment = sequence[selection.start - 1:selection.end]
    shown = segment if len(segment) <= 120 else segment[:117] + "…"
    return (
        '<div class="pe-inspector"><div class="pe-inspector-head"><div>'
        f'<div class="pe-inspector-name">Residues {selection.start}–{selection.end}</div>'
        '<div class="pe-inspector-codes">Side-chain charges only; termini excluded</div>'
        f"</div></div>{_rows(range_summary(segment, ph))}"
        f'<div class="pe-inspector-seq">{escape(shown)}</div></div>'
    )

def residue_inspector(
    sequence: str, selection: ResidueRange | None, scheme: ColorScheme, palette: ChartPalette,
    *, window: int, window_means: list[float | None], ph: float, record_id: str,
    description: str,
) -> None:
    """Render the inspector; offers the selected region as FASTA when a range is selected."""
    residues = sequence.removesuffix("*")
    if selection is None or selection.end > len(residues):
        st.html('<div class="pe-inspector"><p class="pe-empty">Select a residue in the sequence '
                "to inspect it, shift-click for a range, or drag across the hydropathy plot."
                "</p></div>")
        return
    if selection.start == selection.end:
        mean = window_means[selection.start - 1] if window_means else None
        st.html(single_residue_html(residues, selection.start, scheme, palette, mean, window, ph))
        return
    st.html(range_html(residues, selection, ph))
    segment = residues[selection.start - 1:selection.end]
    header = f">{record_id}/{selection.start}-{selection.end} {description}".rstrip()
    st.download_button(
        "Download region as FASTA", f"{header}\n{segment}\n",
        file_name=f"{record_id.split('|')[-1]}_{selection.start}-{selection.end}.fasta",
        mime="text/x-fasta", width="stretch", icon=":material/download:",
    )
