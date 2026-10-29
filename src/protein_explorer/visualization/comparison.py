# Protein comparison: property scatter (with highlighting and grouping) and a normalized
# property heatmap. Features with different units are compared only after z-scoring.

from collections.abc import Sequence

import numpy as np
import pandas as pd
import plotly.graph_objects as go

from protein_explorer.visualization.theme import ChartPalette, style_figure

FEATURE_LABELS: dict[str, str] = {
    "length": "Length (residues)",
    "molecular_weight": "Molecular weight (Da)",
    "isoelectric_point": "Isoelectric point (estimated)",
    "gravy": "GRAVY (Kyte–Doolittle)",
    "aromatic_fraction": "Aromatic fraction",
    "hydrophobic_fraction": "Hydrophobic fraction",
    "charged_fraction": "Charged fraction",
    "acidic_count": "Acidic residues (D+E)",
    "basic_count": "Basic residues (K+R+H)",
    "unique_residue_types": "Distinct residue types",
}
MAX_GROUPS = 3  # categorical colours validated for scatter (all-pairs); the rest fold to Other
OTHER = "Other"

def net_charge_column(features: pd.DataFrame) -> str | None:
    return next((c for c in features.columns if c.startswith("net_charge_ph")), None)

def comparison_features(features: pd.DataFrame) -> list[str]:
    """The comparable properties present in a feature table, in a fixed order."""
    order = ["length", "molecular_weight", "isoelectric_point", "gravy", "aromatic_fraction",
             "hydrophobic_fraction", "charged_fraction"]
    charge = net_charge_column(features)
    return [c for c in order if c in features.columns] + ([charge] if charge else [])

def feature_label(column: str) -> str:
    if column.startswith("net_charge_ph"):
        return f"Net charge at pH {column.removeprefix('net_charge_ph')} (estimated)"
    return FEATURE_LABELS.get(column, column)

def feature_name(column: str) -> str:
    """Label without unit or qualifier, for headings."""
    return feature_label(column).split(" (")[0]

def group_labels(values: pd.Series, max_groups: int = MAX_GROUPS) -> pd.Series:
    """Keep the most common groups and fold the rest (and missing) into 'Other', so a
    scatter never needs more categorical colours than were validated."""
    filled = values.fillna("Unknown")
    keep = filled.value_counts().index[:max_groups]
    return filled.where(filled.isin(keep), OTHER)

def comparison_scatter(
    features: pd.DataFrame,
    x: str,
    y: str,
    palette: ChartPalette,
    *,
    labels: pd.Series | None = None,
    highlighted: Sequence[str] = (),
    selected: str | None = None,
    group_by: str | None = None,
    show_labels: bool = True,
    height: int = 420,
) -> go.Figure:
    """Scatter of two properties. Highlighted proteins (e.g. search matches) are emphasised;
    the selected protein gets the accent colour. customdata[0] is the identifier, so a
    click can be mapped back to a protein."""
    data = features[[x, y]].dropna()
    names = (labels.reindex(data.index) if labels is not None else pd.Series(data.index, index=data.index))
    groups = (group_labels(features.loc[data.index, group_by]) if group_by
              else pd.Series("Proteins", index=data.index))
    dim = bool(highlighted)
    fig = go.Figure()
    for slot, group in enumerate(sorted(groups.unique(), key=lambda g: (g == OTHER, g))):
        members = data.index[groups == group]
        color = palette.muted if group == OTHER else palette.series[slot % MAX_GROUPS]
        opacity = [1.0 if (not dim or i in highlighted) else 0.25 for i in members]
        fig.add_scatter(
            x=data.loc[members, x], y=data.loc[members, y], name=group,
            mode="markers+text" if show_labels else "markers",
            text=names.loc[members] if show_labels else None, textposition="top center",
            textfont={"size": 11, "color": palette.ink_secondary},
            marker={"size": 11, "color": color, "opacity": opacity,
                    "line": {"width": 1.5, "color": palette.surface}},
            customdata=np.column_stack([members, names.loc[members]]),
            hovertemplate=(f"%{{customdata[1]}}<br>{feature_name(x)}: %{{x:,.3~f}}"
                           f"<br>{feature_name(y)}: %{{y:,.3~f}}<extra></extra>"),
        )
    if selected is not None and selected in data.index:
        fig.add_scatter(
            x=[data.loc[selected, x]], y=[data.loc[selected, y]], mode="markers",
            name="Selected protein", hoverinfo="skip",
            marker={"size": 20, "color": "rgba(0,0,0,0)",
                    "line": {"width": 2.5, "color": palette.accent}},
        )
    fig.update_layout(showlegend=bool(group_by) or selected is not None, clickmode="event+select")
    return style_figure(fig, palette, feature_label(x), feature_label(y), height=height)

def _wrap(label: str) -> str:
    """Two-line axis label so long property names stay horizontal."""
    words = label.split()
    if len(words) < 2:
        return label
    middle = (len(words) + 1) // 2
    return " ".join(words[:middle]) + "<br>" + " ".join(words[middle:])

def zscore_table(features: pd.DataFrame, columns: Sequence[str]) -> pd.DataFrame:
    """Z-score each property across the given proteins: (value - mean) / standard deviation
    (sample SD). A property with no spread gets 0. This puts properties with different
    units on one scale; it is relative to the proteins compared, not an absolute score."""
    data = features[list(columns)].astype(float)
    sd = data.std(ddof=1).replace(0, np.nan)
    return ((data - data.mean()) / sd).fillna(0.0)

def normalized_feature_heatmap(
    features: pd.DataFrame,
    columns: Sequence[str],
    palette: ChartPalette,
    *,
    labels: pd.Series | None = None,
    height: int | None = None,
) -> go.Figure:
    """Proteins x properties, coloured by z-score (diverging around the mean), with raw
    values in the hover so units are never lost."""
    data = features[list(columns)].dropna()
    z = zscore_table(data, columns)
    names = (labels.reindex(data.index) if labels is not None else pd.Series(data.index, index=data.index))
    limit = max(2.0, float(np.nanmax(np.abs(z.to_numpy()))) if len(z) else 2.0)
    raw = data.to_numpy()
    fig = go.Figure(go.Heatmap(
        z=z.to_numpy(), x=[_wrap(feature_name(c)) for c in columns], y=names.tolist(),
        zmin=-limit, zmax=limit, zmid=0,
        colorscale=[[0, palette.diverging_low], [0.5, palette.diverging_mid], [1, palette.diverging_high]],
        customdata=raw, xgap=2, ygap=2,
        hovertemplate="%{y}<br>%{x}: %{customdata:,.3~f}<br>z = %{z:+.2f}<extra></extra>",
        colorbar={"title": {"text": "z-score", "font": {"color": palette.ink_secondary}},
                  "tickfont": {"color": palette.ink_secondary}, "outlinewidth": 0},
    ))
    fig.update_yaxes(autorange="reversed", showgrid=False)
    fig.update_xaxes(showgrid=False, side="top", tickangle=0)
    return style_figure(fig, palette, "", "", height=height or max(260, 34 * len(data) + 110))
