"""Tests for figure generation. Checks structure (data, labels, files), not pixels."""

from pathlib import Path

import pandas as pd
import pytest
from matplotlib.figure import Figure

from protein_explorer.dataset import build_composition_table, build_feature_table
from protein_explorer.fasta import FastaRecord, read_fasta
from protein_explorer.features import AMINO_ACIDS
from protein_explorer.plots import (
    RESIDUE_CLASSES,
    RESIDUE_COLORS,
    RESIDUE_ORDER,
    entry_label,
    plot_composition,
    plot_composition_heatmap,
    plot_gravy_comparison,
    plot_hydropathy_profile,
    plot_length_distribution,
    plot_mw_vs_length,
    save_figure,
    write_all_figures,
)

SAMPLE_FASTA = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"


@pytest.fixture(scope="module")
def records() -> list[FastaRecord]:
    return read_fasta(SAMPLE_FASTA)


@pytest.fixture(scope="module")
def features(records) -> pd.DataFrame:
    return build_feature_table(records)


@pytest.fixture(scope="module")
def composition(records) -> pd.DataFrame:
    return build_composition_table(records)


def _axes(fig: Figure):
    return fig.axes[0]


# --- Residue classes ---------------------------------------------------------


def test_residue_classes_partition_the_twenty_amino_acids():
    assert sorted(RESIDUE_ORDER) == list(AMINO_ACIDS)
    assert sum(len(m) for m in RESIDUE_CLASSES.values()) == 20


def test_every_residue_has_a_colour():
    assert set(RESIDUE_COLORS) == set(AMINO_ACIDS)


@pytest.mark.parametrize(
    ("identifier", "expected"),
    [("sp|P01308|INS_HUMAN", "INS_HUMAN"), ("tr|A0A0|A0A0_HUMAN", "A0A0_HUMAN"),
     ("my_protein", "my_protein"), ("a|b|c", "a|b|c")],
)
def test_entry_label(identifier, expected):
    assert entry_label(identifier) == expected


# --- Individual figures ------------------------------------------------------


def test_composition_bars_follow_class_order(composition):
    insulin = composition.loc["sp|P01308|INS_HUMAN"]
    ax = _axes(plot_composition(insulin, "INS_HUMAN composition"))
    assert [t.get_text() for t in ax.get_xticklabels()] == list(RESIDUE_ORDER)
    heights = [bar.get_height() for bar in ax.patches]
    assert heights[RESIDUE_ORDER.index("L")] == pytest.approx(insulin["L"] * 100)
    assert ax.get_ylabel() == "Frequency (%)"


def test_composition_with_reference_has_both_legends(composition):
    fig = plot_composition(composition.iloc[0], "t", reference=composition.mean())
    texts = [t.get_text() for legend in _axes(fig).findobj(lambda a: hasattr(a, "get_texts"))
             for t in legend.get_texts()]
    assert "Dataset mean" in texts
    assert "Nonpolar" in texts


def test_length_distribution_counts_all_proteins(features):
    ax = _axes(plot_length_distribution(features))
    assert sum(bar.get_height() for bar in ax.patches) == 8
    assert "residues" in ax.get_xlabel()


def test_gravy_comparison_sorted_and_coloured_by_sign(features):
    ax = _axes(plot_gravy_comparison(features))
    labels = [t.get_text() for t in ax.get_yticklabels()]
    assert labels[0] == "P53_HUMAN"  # most hydrophilic at the bottom
    assert labels[-1] == "INS_HUMAN"  # most hydrophobic at the top
    widths = [bar.get_width() for bar in ax.patches[:8]]
    assert widths == sorted(widths)


def test_mw_vs_length_labels_every_point(features):
    ax = _axes(plot_mw_vs_length(features))
    annotations = {t.get_text() for t in ax.texts}
    assert {"INS_HUMAN", "ALBU_BOVIN", "CALM1_HUMAN"} <= annotations
    assert ax.get_ylabel() == "Molecular weight (kDa)"


def test_mw_vs_length_labels_do_not_overlap(features):
    ax = _axes(plot_mw_vs_length(features))
    label_ys = sorted(t.xyann[1] for t in ax.texts)
    gaps = [b - a for a, b in zip(label_ys, label_ys[1:])]
    assert min(gaps) > 0


def test_composition_heatmap_shape(composition):
    fig = plot_composition_heatmap(composition)
    image = _axes(fig).images[0]
    assert image.get_array().shape == (8, 20)
    assert len(fig.axes) == 2  # heatmap + colourbar


def test_hydropathy_profile_plot(records):
    insulin = next(r for r in records if "INS_HUMAN" in r.identifier)
    ax = _axes(plot_hydropathy_profile(insulin.sequence, "t", window=19))
    line = ax.lines[-1]
    assert len(line.get_xdata()) == 110 - 19 + 1
    assert "window 19" in ax.get_ylabel()


def test_invalid_records_excluded_from_dataset_plots():
    records = [FastaRecord("good", "", "MKTAYIAKQR"), FastaRecord("bad", "", "MK-TAYI")]
    ax = _axes(plot_gravy_comparison(build_feature_table(records)))
    assert [t.get_text() for t in ax.get_yticklabels()] == ["good"]


# --- Saving ------------------------------------------------------------------


@pytest.mark.parametrize("suffix", ["png", "svg", "pdf"])
def test_save_figure_formats(tmp_path, features, suffix):
    path = save_figure(plot_length_distribution(features), tmp_path / "sub" / f"len.{suffix}")
    assert path.exists() and path.stat().st_size > 0


def test_write_all_figures(tmp_path, records):
    paths = write_all_figures(records, tmp_path, window=9)
    # 4 dataset figures + composition and hydropathy for each of 8 proteins
    assert len(paths) == 4 + 2 * 8
    assert (tmp_path / "composition" / "ALBU_BOVIN.png").exists()


def test_write_all_figures_skips_invalid_and_short(tmp_path):
    records = [
        FastaRecord("long", "", "MKTAYIAKQRDEHW"),
        FastaRecord("short", "", "MKT"),  # shorter than window: composition only
        FastaRecord("bad", "", "MK-TAYI"),  # invalid: nothing per-protein
    ]
    names = {p.relative_to(tmp_path).as_posix() for p in write_all_figures(records, tmp_path)}
    assert "hydropathy/long.png" in names
    assert "composition/short.png" in names
    assert "hydropathy/short.png" not in names
    assert not any("bad" in name for name in names)
