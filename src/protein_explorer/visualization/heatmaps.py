# Amino acid composition heatmap (proteins x residues) with row ordering options,
# including hierarchical clustering.

import pandas as pd
import plotly.graph_objects as go
from scipy.cluster.hierarchy import leaves_list, linkage

from protein_explorer.residues import CHEMICAL_CLASS, NAMES, RESIDUE_ORDER
from protein_explorer.visualization.theme import ChartPalette, style_figure

ROW_ORDERS = {
    "input": "File order",
    "name": "Name",
    "cluster": "Clustered by composition",
}

def order_rows(composition: pd.DataFrame, method: str, labels: pd.Series | None = None) -> list:
    """Row order for the heatmap.

    input: as in the file. name: by display label. cluster: average-linkage hierarchical
    clustering on Euclidean distance between composition vectors, so proteins with similar
    composition sit together (needs at least 3 proteins to be informative). A one-letter
    residue code sorts by that residue's frequency, highest first.
    """
    rows = list(composition.index)
    if method == "input" or len(rows) < 2:
        return rows
    if method == "name":
        keys = labels.reindex(rows) if labels is not None else pd.Series(rows, index=rows)
        return sorted(rows, key=lambda r: str(keys[r]).lower())
    if method == "cluster":
        if len(rows) < 3:
            return rows
        tree = linkage(composition.loc[rows].to_numpy(dtype=float), method="average",
                       metric="euclidean")
        return [rows[i] for i in leaves_list(tree)]
    if method in RESIDUE_ORDER:
        return composition[method].sort_values(ascending=False, kind="stable").index.tolist()
    raise ValueError(f"Unknown row order '{method}'.")

def composition_heatmap(
    composition: pd.DataFrame,
    palette: ChartPalette,
    *,
    order: str = "input",
    labels: pd.Series | None = None,
    height: int | None = None,
) -> go.Figure:
    """Residue frequencies (%) per protein; columns grouped by chemical class. Rows with
    missing values (invalid sequences) are left out."""
    data = composition.dropna().reindex(columns=list(RESIDUE_ORDER)) * 100
    rows = order_rows(data, order, labels)
    data = data.loc[rows]
    names = (labels.reindex(rows) if labels is not None else pd.Series(rows, index=rows)).tolist()
    steps = len(palette.sequential) - 1
    hover = [[f"{NAMES[r]}, {CHEMICAL_CLASS[r].lower()}" for r in RESIDUE_ORDER]] * len(data)
    fig = go.Figure(go.Heatmap(
        z=data.to_numpy(), x=list(RESIDUE_ORDER), y=names, zmin=0,
        colorscale=[[i / steps, c] for i, c in enumerate(palette.sequential)],
        customdata=hover, xgap=1, ygap=1,
        hovertemplate="%{y}<br>%{x} (%{customdata}): %{z:.1f}%<extra></extra>",
        colorbar={"title": {"text": "Frequency (%)", "font": {"color": palette.ink_secondary}},
                  "tickfont": {"color": palette.ink_secondary}, "outlinewidth": 0},
    ))
    fig.update_yaxes(autorange="reversed", showgrid=False)
    fig.update_xaxes(showgrid=False)
    return style_figure(fig, palette, "Amino acid (grouped by chemical class)", "",
                        height=height or max(260, 30 * len(data) + 120))
