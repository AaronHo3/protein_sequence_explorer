# Interactive hydropathy profile. Distinguishes per-residue Kyte-Doolittle values (bars)
# from the sliding-window mean (line), and supports linked residue selection.

import math

import plotly.graph_objects as go

from protein_explorer.features import KYTE_DOOLITTLE, hydropathy_profile
from protein_explorer.residues import THREE_LETTER
from protein_explorer.visualization.sequence_viewer import ResidueRange
from protein_explorer.visualization.theme import ChartPalette, mix, style_figure

KD_REFERENCE = 1.6  # Kyte & Doolittle's cutoff, defined for window 19; a reference, not a call
WINDOW_TRACE = "Window mean"

def window_bounds(length: int, window: int) -> tuple[int, int] | None:
    """First and last residue positions that have a full window, or None if none do."""
    if length < window:
        return None
    half = window // 2
    return half + 1, length - half

def hydropathy_figure(
    sequence: str,
    palette: ChartPalette,
    window: int = 9,
    *,
    show_residue_values: bool = True,
    show_reference: bool = False,
    selection: ResidueRange | None = None,
    height: int = 360,
) -> go.Figure:
    """Window-mean line with positive (hydrophobic) and negative (hydrophilic) regions
    shaded, optional per-residue bars, and the current residue selection marked.

    Window means exist only for positions with a full window; the first and last
    window // 2 residues have none, and the line is not extrapolated into them.
    """
    residues = sequence.upper().removesuffix("*")
    fig = go.Figure()

    if show_residue_values:
        values = [KYTE_DOOLITTLE.get(r) for r in residues]
        fig.add_bar(
            x=list(range(1, len(residues) + 1)), y=values, name="Residue value",
            marker={"color": mix(palette.muted, palette.surface, 0.62), "line": {"width": 0}},
            opacity=0.75,
            customdata=[THREE_LETTER.get(r, r) for r in residues],
            hovertemplate="%{customdata} %{x}: %{y:+.1f} (residue)<extra></extra>",
        )

    profile = hydropathy_profile(residues, window) if len(residues) >= window else []
    if profile:
        x = [p for p, _ in profile]
        y = [None if math.isnan(v) else v for _, v in profile]
        positive = [None if v is None else max(v, 0) for v in y]
        negative = [None if v is None else min(v, 0) for v in y]
        fig.add_scatter(x=x, y=positive, fill="tozeroy", mode="none", hoverinfo="skip",
                        fillcolor=mix(palette.hydrophobic, palette.surface, 0.72),
                        name="Hydrophobic (above 0)")
        fig.add_scatter(x=x, y=negative, fill="tozeroy", mode="none", hoverinfo="skip",
                        fillcolor=mix(palette.hydrophilic, palette.surface, 0.72),
                        name="Hydrophilic (below 0)")
        half = window // 2
        fig.add_scatter(
            x=x, y=y, mode="lines+markers", name=f"{WINDOW_TRACE} (window {window})",
            line={"color": palette.ink, "width": 1.8},
            marker={"size": 4, "color": palette.ink, "opacity": 0.35},
            customdata=[[p - half, p + half] for p in x],
            hovertemplate=("Residues %{customdata[0]}–%{customdata[1]}, centre %{x}<br>"
                           "Window mean %{y:+.2f}<extra></extra>"),
        )
    else:
        fig.add_annotation(
            text=f"Sequence is shorter than the {window}-residue window; "
                 "only per-residue values are shown.",
            xref="paper", yref="paper", x=0.5, y=1.0, showarrow=False,
            font={"color": palette.ink_secondary, "size": 12},
        )

    fig.add_hline(y=0, line={"color": palette.baseline, "width": 1.2})
    if show_reference:
        fig.add_hline(y=KD_REFERENCE, line={"color": palette.muted, "width": 1, "dash": "dash"},
                      annotation_text=f"KD reference {KD_REFERENCE} (defined for window 19)",
                      annotation_position="top left",
                      annotation_font={"color": palette.ink_secondary, "size": 11})
    if selection is not None:
        fig.add_vrect(x0=selection.start - 0.5, x1=selection.end + 0.5, line_width=0,
                      fillcolor=palette.selection, opacity=0.16, layer="below")
        fig.add_vline(x=selection.start, line={"color": palette.selection, "width": 1.5})
        if selection.end != selection.start:
            fig.add_vline(x=selection.end, line={"color": palette.selection, "width": 1.5})

    fig.update_layout(
        dragmode="select", selectdirection="h", hovermode="closest", bargap=0.05,
        clickmode="event+select",
    )
    fig.update_xaxes(range=[0.5, len(residues) + 0.5], showgrid=False)
    fig.update_yaxes(range=[-4.8, 4.8], zeroline=False)
    return style_figure(fig, palette, "Residue position",
                        "Kyte–Doolittle hydropathy", height=height)

def selection_from_points(points: list[dict]) -> ResidueRange | None:
    """Residue range covered by points selected on the hydropathy chart (any trace)."""
    xs = [int(round(p["x"])) for p in points if p.get("x") is not None]
    return ResidueRange.from_positions(xs) if xs else None
