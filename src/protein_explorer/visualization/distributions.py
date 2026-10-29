# Distribution charts: histograms with mean/median and a rug of individual proteins (so
# small datasets stay readable), and molecular weight vs length with a linear fit.

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from protein_explorer.visualization.theme import ChartPalette, style_figure

def default_bins(n: int) -> int:
    """Sturges' rule, capped so tiny datasets still get a handful of bins."""
    return max(1, min(int(np.ceil(np.log2(max(n, 1)) + 1)), 50))

def distribution_histogram(
    values: pd.Series,
    palette: ChartPalette,
    *,
    xlabel: str,
    bins: int | None = None,
    labels: pd.Series | None = None,
    value_format: str = ".1f",
    height: int = 320,
) -> go.Figure:
    """Histogram with mean and median lines and a rug marking each protein.

    The rug matters for small datasets, where a histogram alone hides how few points
    there are. labels: per-value names for the rug hover.
    """
    data = values.dropna().astype(float)
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, row_heights=[0.86, 0.14],
                        vertical_spacing=0.03)
    if data.empty:
        fig.add_annotation(text="No values to show.", showarrow=False, xref="paper",
                           yref="paper", x=0.5, y=0.5)
        return style_figure(fig, palette, xlabel, "Number of proteins", height=height)

    edges = np.histogram_bin_edges(data, bins=bins or default_bins(len(data)))
    if edges[0] == edges[-1]:  # every value identical
        edges = np.array([edges[0] - 0.5, edges[0] + 0.5])
    counts, _ = np.histogram(data, bins=edges)
    centers = (edges[:-1] + edges[1:]) / 2
    fig.add_bar(
        x=centers, y=counts, width=np.diff(edges) * 0.94, name="Proteins",
        marker={"color": palette.series[0]},
        customdata=np.column_stack([edges[:-1], edges[1:]]),
        hovertemplate=(f"%{{customdata[0]:{value_format}}} to %{{customdata[1]:{value_format}}}"
                       "<br>%{y} proteins<extra></extra>"),
        row=1, col=1,
    )
    names = labels.reindex(data.index).fillna("").tolist() if labels is not None else [""] * len(data)
    fig.add_scatter(
        x=data, y=[0] * len(data), mode="markers", name="Individual proteins",
        marker={"symbol": "line-ns-open", "size": 14, "color": palette.ink_secondary,
                "line": {"width": 1.5, "color": palette.ink_secondary}},
        text=names, hovertemplate=f"%{{text}}<br>%{{x:{value_format}}}<extra></extra>",
        row=2, col=1,
    )
    mean, median = data.mean(), data.median()
    # Put each label on the side away from the other line so close values don't collide
    for stat, value, dash, side in (
        ("Mean", mean, "dash", "top right" if mean >= median else "top left"),
        ("Median", median, "dot", "top left" if mean >= median else "top right"),
    ):
        fig.add_vline(x=value, line={"color": palette.ink, "width": 1.2, "dash": dash},
                      annotation_text=f"{stat} {value:{value_format}}",
                      annotation_position=side, annotation_font={"size": 11, "color": palette.ink_secondary},
                      row=1, col=1)
    fig.update_layout(showlegend=False, bargap=0)
    style_figure(fig, palette, "", "", height=height)
    fig.update_yaxes(title_text="Number of proteins", rangemode="tozero", row=1, col=1)
    fig.update_yaxes(visible=False, row=2, col=1)
    fig.update_xaxes(title_text=xlabel, row=2, col=1)
    fig.update_xaxes(showgrid=False)
    return fig

def linear_fit(x: pd.Series, y: pd.Series) -> tuple[float, float, float] | None:
    """Least-squares slope, intercept, and R² (None with fewer than 3 points)."""
    if len(x) < 3:
        return None
    slope, intercept = np.polyfit(x.astype(float), y.astype(float), 1)
    predicted = slope * x + intercept
    ss_res = float(((y - predicted) ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r_squared = 1 - ss_res / ss_tot if ss_tot else 1.0
    return float(slope), float(intercept), r_squared

def mw_vs_length(
    features: pd.DataFrame,
    palette: ChartPalette,
    *,
    labels: pd.Series | None = None,
    highlight: str | None = None,
    show_fit: bool = True,
    height: int = 360,
) -> go.Figure:
    """Molecular weight (kDa) vs length (residues), with an optional least-squares fit."""
    data = features[["length", "molecular_weight"]].dropna().astype(float)
    names = (labels.reindex(data.index) if labels is not None else pd.Series(data.index, index=data.index)).tolist()
    kda = data["molecular_weight"] / 1000
    fig = go.Figure()
    fit = linear_fit(data["length"], data["molecular_weight"]) if show_fit else None
    if fit is not None:
        slope, intercept, r2 = fit
        xs = np.array([data["length"].min(), data["length"].max()])
        fig.add_scatter(
            x=xs, y=(slope * xs + intercept) / 1000, mode="lines",
            name=f"Linear fit: {slope:.1f} Da per residue (R² = {r2:.4f})",
            line={"color": palette.muted, "dash": "dash", "width": 1.2}, hoverinfo="skip",
        )
    is_highlight = data.index == highlight
    fig.add_scatter(
        x=data["length"][~is_highlight], y=kda[~is_highlight], mode="markers", name="Proteins",
        text=[n for n, h in zip(names, is_highlight, strict=True) if not h],
        marker={"size": 10, "color": palette.series[0], "line": {"width": 1.5, "color": palette.surface}},
        hovertemplate="%{text}<br>%{x} residues<br>%{y:.2f} kDa<extra></extra>",
    )
    if is_highlight.any():
        fig.add_scatter(
            x=data["length"][is_highlight], y=kda[is_highlight], mode="markers+text",
            name="Selected protein", text=[n for n, h in zip(names, is_highlight, strict=True) if h],
            textposition="top center", textfont={"color": palette.ink},
            marker={"size": 13, "color": palette.accent, "line": {"width": 2, "color": palette.surface}},
            hovertemplate="%{text}<br>%{x} residues<br>%{y:.2f} kDa<extra></extra>",
        )
    return style_figure(fig, palette, "Length (residues)", "Molecular weight (kDa)", height=height)
