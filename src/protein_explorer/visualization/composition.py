# Amino acid composition charts: one protein (counts or percentages) and dataset-level
# composition (residue- vs protein-weighted).

import pandas as pd
import plotly.graph_objects as go

from protein_explorer.residues import CHEMICAL_CLASS, NAMES, RESIDUE_CLASSES, RESIDUE_ORDER
from protein_explorer.visualization.theme import ChartPalette, style_figure

SORT_ORDERS = {"class": "Chemical class", "alphabetical": "Alphabetical", "frequency": "Frequency"}

def residue_order(values: pd.Series, order: str) -> list[str]:
    """Residue display order: grouped by chemical class, alphabetical, or by descending value
    (ties broken alphabetically so the order is stable)."""
    if order == "class":
        return list(RESIDUE_ORDER)
    if order == "alphabetical":
        return sorted(RESIDUE_ORDER)
    if order == "frequency":
        frame = pd.DataFrame({"r": list(RESIDUE_ORDER), "v": values.reindex(list(RESIDUE_ORDER)).fillna(0).to_numpy()})
        return frame.sort_values(["v", "r"], ascending=[False, True])["r"].tolist()
    raise ValueError(f"Unknown order '{order}'. Choose from {', '.join(SORT_ORDERS)}.")

def _class_colors(palette: ChartPalette) -> dict[str, str]:
    return dict(zip(RESIDUE_CLASSES, palette.series, strict=False))

def composition_bar_chart(
    counts: pd.Series,
    length: int,
    palette: ChartPalette,
    *,
    as_percent: bool = True,
    order: str = "class",
    reference: pd.Series | None = None,
    height: int = 340,
) -> go.Figure:
    """Bars for the 20 standard residues, coloured by chemical class.

    counts: residue -> count. length: sequence length used as the percentage denominator
    (it includes ambiguous residues, matching residue_frequencies).
    reference: residue -> fraction (e.g. dataset composition), drawn as ticks in the same
    units as the bars.
    """
    values = counts.reindex(list(RESIDUE_ORDER)).fillna(0).astype(float)
    shown = values / length * 100 if as_percent else values
    axis_order = residue_order(shown, order)
    colors = _class_colors(palette)
    fig = go.Figure()
    for class_name, members in RESIDUE_CLASSES.items():
        residues = [r for r in axis_order if r in members]
        fig.add_bar(
            x=residues, y=[shown[r] for r in residues], name=class_name,
            marker={"color": colors[class_name], "line": {"width": 0}},
            customdata=[[NAMES[r], int(values[r]), values[r] / length * 100] for r in residues],
            hovertemplate=("%{customdata[0]} (%{x})<br>%{customdata[1]} residues, "
                           "%{customdata[2]:.1f}%<extra>" + class_name + "</extra>"),
        )
    if reference is not None:
        ref = reference.reindex(axis_order).astype(float)
        ref_values = ref * 100 if as_percent else ref * length
        fig.add_scatter(
            x=axis_order, y=ref_values, mode="markers", name="Dataset reference",
            marker={"symbol": "line-ew-open", "size": 20, "color": palette.ink,
                    "line": {"width": 2, "color": palette.ink}},
            hovertemplate=("Dataset reference %{x}: %{y:.1f}" + ("%" if as_percent else " residues")
                           + "<extra></extra>"),
        )
    fig.update_layout(bargap=0.18)
    fig.update_xaxes(categoryorder="array", categoryarray=axis_order, showgrid=False)
    ylabel = "Frequency (%)" if as_percent else "Count (residues)"
    return style_figure(fig, palette, "Amino acid", ylabel, height=height)

def dataset_composition_chart(
    composition: pd.DataFrame, palette: ChartPalette, *, order: str = "class", height: int = 340
) -> go.Figure:
    """Grouped bars: residue-weighted vs protein-weighted dataset frequencies (%)."""
    data = composition.reindex(list(RESIDUE_ORDER)) * 100
    axis_order = residue_order(data["residue_weighted"], order)
    fig = go.Figure()
    styles = {
        "residue_weighted": ("Residue-weighted (pooled)", palette.series[0]),
        "protein_weighted": ("Protein-weighted (mean of proteins)", palette.series[1]),
    }
    for column, (label, color) in styles.items():
        fig.add_bar(
            x=axis_order, y=data.loc[axis_order, column], name=label, marker={"color": color},
            customdata=[[NAMES[r], CHEMICAL_CLASS[r]] for r in axis_order],
            hovertemplate=("%{customdata[0]} (%{x}), %{customdata[1]}<br>%{y:.2f}%<extra>"
                           + label + "</extra>"),
        )
    fig.update_layout(barmode="group", bargap=0.22, bargroupgap=0.06)
    fig.update_xaxes(categoryorder="array", categoryarray=axis_order, showgrid=False)
    return style_figure(fig, palette, "Amino acid", "Frequency (%)", height=height)
