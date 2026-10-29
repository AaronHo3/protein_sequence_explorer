# Static, publication-quality figures built with matplotlib's object-oriented API
# (no pyplot), so no global figure state and no GUI backend is needed.

import re
from pathlib import Path

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.patches import Patch

from protein_explorer.dataset import build_composition_table, build_feature_table
from protein_explorer.fasta import FastaRecord
from protein_explorer.features import hydropathy_profile
from protein_explorer.residues import RESIDUE_CLASSES, RESIDUE_ORDER  # noqa: F401 (re-exported)
from protein_explorer.visualization.theme import LIGHT

# Static figures use the light palette (publication / print)
SURFACE = LIGHT.surface
INK = LIGHT.ink
INK_SECONDARY = LIGHT.ink_secondary
MUTED = LIGHT.muted
GRIDLINE = LIGHT.gridline
BASELINE = LIGHT.baseline
SERIES = LIGHT.series
NEGATIVE, POSITIVE = LIGHT.diverging_low, LIGHT.diverging_high
SEQUENTIAL = LinearSegmentedColormap.from_list("seq_blue", list(LIGHT.sequential))
CLASS_COLORS = dict(zip(RESIDUE_CLASSES, LIGHT.series, strict=True))
RESIDUE_COLORS = {r: CLASS_COLORS[c] for c, members in RESIDUE_CLASSES.items() for r in members}

MEAN_RESIDUE_MASS = 110.0  # Da, common rule of thumb for MW ≈ 110 × length
TM_THRESHOLD = 1.6  # Kyte-Doolittle transmembrane cutoff (window 19)
DPI = 300

def entry_label(identifier: str) -> str:
    """Short label: UniProt entry name (INS_HUMAN) if present, else the identifier."""
    parts = identifier.split("|")
    return parts[2] if len(parts) >= 3 and parts[0] in {"sp", "tr"} and parts[2] else identifier

def _new_axes(width: float, height: float) -> tuple[Figure, Axes]:
    """Figure with one axes in the shared style: recessive spines and gridlines."""
    fig = Figure(figsize=(width, height), facecolor=SURFACE, layout="constrained")
    ax = fig.add_subplot()
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(BASELINE)
    ax.tick_params(colors=INK_SECONDARY, labelsize=9)
    ax.grid(color=GRIDLINE, linewidth=0.6)
    ax.set_axisbelow(True)
    return fig, ax

def _label(ax: Axes, title: str, xlabel: str, ylabel: str) -> None:
    ax.set_title(title, color=INK, fontsize=11, loc="left", fontweight="bold")
    ax.set_xlabel(xlabel, color=INK_SECONDARY, fontsize=9)
    ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize=9)

def _class_legend(ax: Axes) -> None:
    handles = [Patch(color=color, label=name) for name, color in CLASS_COLORS.items()]
    ax.legend(handles=handles, frameon=False, fontsize=8, labelcolor=INK_SECONDARY,
              ncols=len(handles), loc="upper left", bbox_to_anchor=(0, 1.0))

def plot_composition(
    frequencies: pd.Series, title: str, reference: pd.Series | None = None
) -> Figure:
    """Bar chart of one protein's residue frequencies, grouped by chemical class.

    If reference is given (e.g. the dataset mean), it is drawn as a short tick per residue.
    """
    fig, ax = _new_axes(8, 3.6)
    values = frequencies.reindex(list(RESIDUE_ORDER)) * 100
    positions = range(len(RESIDUE_ORDER))
    ax.bar(positions, values, width=0.8, color=[RESIDUE_COLORS[r] for r in RESIDUE_ORDER],
           edgecolor=SURFACE, linewidth=1)
    if reference is not None:
        ref = reference.reindex(list(RESIDUE_ORDER)) * 100
        ax.scatter(positions, ref, marker="_", s=220, color=INK, linewidths=1.5,
                   label="Dataset mean", zorder=3)
    ax.set_xticks(list(positions), RESIDUE_ORDER)
    ax.grid(axis="x", visible=False)
    ax.set_ylim(0, max(values.max(), 0 if reference is None else ref.max()) * 1.25)
    _class_legend(ax)
    if reference is not None:
        ax.add_artist(ax.get_legend())
        ax.legend(frameon=False, fontsize=8, labelcolor=INK_SECONDARY, loc="upper right")
    _label(ax, title, "Amino acid", "Frequency (%)")
    return fig

def plot_length_distribution(features: pd.DataFrame, bins: int | str = "auto") -> Figure:
    """Histogram of sequence lengths across proteins."""
    fig, ax = _new_axes(6, 3.6)
    ax.hist(features["length"].dropna().astype(int), bins=bins, color=SERIES[0],
            edgecolor=SURFACE, linewidth=2)
    ax.grid(axis="x", visible=False)
    _label(ax, "Sequence length distribution", "Length (residues)", "Number of proteins")
    return fig

def plot_gravy_comparison(features: pd.DataFrame) -> Figure:
    """Horizontal bars of GRAVY per protein, sorted; red = hydrophobic, blue = hydrophilic."""
    data = features["gravy"].dropna().sort_values()
    fig, ax = _new_axes(6, 0.4 * len(data) + 1.2)
    colors = [POSITIVE if value > 0 else NEGATIVE for value in data]
    ax.barh([entry_label(i) for i in data.index], data, color=colors, height=0.7)
    ax.axvline(0, color=BASELINE, linewidth=1)
    ax.grid(axis="y", visible=False)
    ax.legend(handles=[Patch(color=POSITIVE, label="Hydrophobic (> 0)"),
                       Patch(color=NEGATIVE, label="Hydrophilic (< 0)")],
              frameon=False, fontsize=8, labelcolor=INK_SECONDARY, loc="upper left")
    _label(ax, "Hydropathy by protein (GRAVY)", "GRAVY (Kyte–Doolittle, mean)", "")
    return fig

def plot_mw_vs_length(features: pd.DataFrame) -> Figure:
    """Scatter of molecular weight against length, with a 110 Da/residue reference line."""
    data = features[["length", "molecular_weight"]].dropna()
    kda = data["molecular_weight"] / 1000
    fig, ax = _new_axes(6, 4.2)
    max_length = data["length"].max() * 1.08
    ax.plot([0, max_length], [0, max_length * MEAN_RESIDUE_MASS / 1000], color=MUTED,
            linestyle="--", linewidth=1, label=f"{MEAN_RESIDUE_MASS:.0f} Da per residue")
    ax.scatter(data["length"], kda, s=40, color=SERIES[0], edgecolor=SURFACE,
               linewidth=1.5, zorder=3)
    # Proteins of similar size overlap (e.g. three ~148-residue proteins), so labels are
    # stacked upward with a minimum gap and joined to their point by a leader line.
    min_gap = kda.max() * 0.05
    label_x_offset = max_length * 0.04
    previous_label_y = float("-inf")
    for identifier in kda.sort_values().index:
        length, mass = data.loc[identifier, "length"], kda[identifier]
        label_y = max(mass, previous_label_y + min_gap)
        previous_label_y = label_y
        ax.annotate(
            entry_label(identifier), (length, mass), xytext=(length + label_x_offset, label_y),
            textcoords="data", fontsize=7.5, color=INK_SECONDARY, va="center",
            arrowprops={"arrowstyle": "-", "color": BASELINE, "linewidth": 0.6,
                        "shrinkA": 0, "shrinkB": 3},
        )
    ax.set_xlim(0, max_length)
    ax.set_ylim(0, None)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_SECONDARY, loc="upper left")
    _label(ax, "Molecular weight vs. length", "Length (residues)", "Molecular weight (kDa)")
    return fig

def plot_composition_heatmap(composition: pd.DataFrame) -> Figure:
    """Heatmap of residue frequencies: proteins × residues, residues grouped by class."""
    data = composition.dropna().reindex(columns=list(RESIDUE_ORDER)) * 100
    fig, ax = _new_axes(9, 0.38 * len(data) + 1.6)
    ax.grid(False)
    image = ax.imshow(data.to_numpy(), cmap=SEQUENTIAL, aspect="auto", vmin=0)
    ax.set_xticks(range(len(RESIDUE_ORDER)), RESIDUE_ORDER)
    ax.set_yticks(range(len(data)), [entry_label(i) for i in data.index])
    for residue, tick in zip(RESIDUE_ORDER, ax.get_xticklabels(), strict=True):
        tick.set_color(RESIDUE_COLORS[residue])
        tick.set_fontweight("bold")
    colorbar = fig.colorbar(image, ax=ax, shrink=0.9)
    colorbar.set_label("Frequency (%)", color=INK_SECONDARY, fontsize=9)
    colorbar.outline.set_visible(False)
    colorbar.ax.tick_params(colors=INK_SECONDARY, labelsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    _label(ax, "Amino acid composition", "Amino acid (grouped by class)", "")
    return fig

def plot_hydropathy_profile(sequence: str, title: str, window: int = 9) -> Figure:
    """Sliding-window Kyte-Doolittle profile along the sequence."""
    positions, values = zip(*hydropathy_profile(sequence, window), strict=True)
    fig, ax = _new_axes(8, 3.6)
    ax.axhline(0, color=BASELINE, linewidth=1)
    ax.axhline(TM_THRESHOLD, color=MUTED, linestyle="--", linewidth=1,
               label=f"Transmembrane threshold ({TM_THRESHOLD}, window 19)")
    ax.plot(positions, values, color=SERIES[0], linewidth=1.6)
    ax.set_xlim(1, len(sequence.removesuffix("*")))
    ax.legend(frameon=False, fontsize=8, labelcolor=INK_SECONDARY, loc="lower right")
    _label(ax, title, "Residue position",
           f"Hydropathy (Kyte–Doolittle, window {window})")
    return fig

def save_figure(fig: Figure, path: Path) -> Path:
    """Save a figure at publication resolution; format comes from the file suffix."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI, facecolor=SURFACE)
    return path

def _safe_filename(label: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", label)

def write_all_figures(
    records: list[FastaRecord], output_dir: Path, window: int = 9, fmt: str = "png"
) -> list[Path]:
    """Dataset-level figures plus a composition chart and hydropathy profile per valid protein."""
    features = build_feature_table(records)
    composition = build_composition_table(records)
    dataset_figures = {
        "length_distribution": plot_length_distribution(features),
        "gravy_comparison": plot_gravy_comparison(features),
        "mw_vs_length": plot_mw_vs_length(features),
        "composition_heatmap": plot_composition_heatmap(composition),
    }
    paths = [save_figure(fig, output_dir / f"{name}.{fmt}") for name, fig in dataset_figures.items()]

    mean_composition = composition.mean()
    for record in records:
        if not features.loc[record.identifier, "is_valid"]:
            continue
        label = entry_label(record.identifier)
        name = _safe_filename(label)
        paths.append(save_figure(
            plot_composition(composition.loc[record.identifier], f"{label} composition",
                             reference=mean_composition),
            output_dir / "composition" / f"{name}.{fmt}",
        ))
        if len(record.sequence.removesuffix("*")) >= window:
            paths.append(save_figure(
                plot_hydropathy_profile(record.sequence, f"{label} hydropathy profile", window),
                output_dir / "hydropathy" / f"{name}.{fmt}",
            ))
    return paths
