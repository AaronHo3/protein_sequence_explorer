"""Tests for Plotly figure builders: the data and labels behind each chart."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from protein_explorer.dataset import (
    build_composition_table,
    build_count_table,
    build_feature_table,
    dataset_composition,
)
from protein_explorer.fasta import FastaRecord, read_fasta
from protein_explorer.features import residue_counts
from protein_explorer.residues import RESIDUE_ORDER
from protein_explorer.visualization.comparison import (
    comparison_features,
    comparison_scatter,
    feature_label,
    group_labels,
    normalized_feature_heatmap,
    zscore_table,
)
from protein_explorer.visualization.composition import (
    composition_bar_chart,
    dataset_composition_chart,
    residue_order,
)
from protein_explorer.visualization.distributions import (
    default_bins,
    distribution_histogram,
    linear_fit,
    mw_vs_length,
)
from protein_explorer.visualization.heatmaps import composition_heatmap, order_rows
from protein_explorer.visualization.hydropathy import (
    hydropathy_figure,
    selection_from_points,
    window_bounds,
)
from protein_explorer.visualization.sequence_viewer import ResidueRange
from protein_explorer.visualization.theme import DARK, LIGHT

SAMPLE = Path(__file__).resolve().parent.parent / "data" / "sample.fasta"
INSULIN = "sp|P01308|INS_HUMAN"


@pytest.fixture(scope="module")
def records():
    return read_fasta(SAMPLE)


@pytest.fixture(scope="module")
def features(records):
    return build_feature_table(records)


@pytest.fixture(scope="module")
def composition(records):
    return build_composition_table(records)


def _sequence(records, identifier=INSULIN):
    return next(r.sequence for r in records if r.identifier == identifier)


# --- Hydropathy ----------------------------------------------------------------


def test_window_bounds():
    assert window_bounds(110, 9) == (5, 106)
    assert window_bounds(5, 9) is None


def test_hydropathy_distinguishes_residue_values_from_window_mean(records):
    fig = hydropathy_figure(_sequence(records), LIGHT, 9)
    names = [t.name for t in fig.data]
    assert "Residue value" in names and "Window mean (window 9)" in names
    bars = next(t for t in fig.data if t.name == "Residue value")
    line = next(t for t in fig.data if t.name.startswith("Window mean"))
    assert len(bars.x) == 110
    assert line.x[0] == 5 and line.x[-1] == 106  # not extrapolated to the ends


def test_hydropathy_shades_positive_and_negative_separately(records):
    fig = hydropathy_figure(_sequence(records), LIGHT, 9)
    positive = next(t for t in fig.data if t.name.startswith("Hydrophobic"))
    negative = next(t for t in fig.data if t.name.startswith("Hydrophilic"))
    assert min(v for v in positive.y if v is not None) >= 0
    assert max(v for v in negative.y if v is not None) <= 0


def test_hydropathy_reference_line_off_by_default(records):
    # Peaks are not labelled as transmembrane unless the user asks for the reference
    seq = _sequence(records)
    def has_reference(fig):
        return any(s.y0 == 1.6 for s in fig.layout.shapes)
    assert not has_reference(hydropathy_figure(seq, LIGHT))
    assert has_reference(hydropathy_figure(seq, LIGHT, show_reference=True))


def test_hydropathy_marks_selection(records):
    fig = hydropathy_figure(_sequence(records), LIGHT, selection=ResidueRange(5, 18))
    rects = [s for s in fig.layout.shapes if s.type == "rect"]
    assert rects and rects[0].x0 == 4.5 and rects[0].x1 == 18.5


def test_hydropathy_short_sequence_has_no_window_line():
    fig = hydropathy_figure("MKTAY", LIGHT, 9)
    assert not any((t.name or "").startswith("Window mean") for t in fig.data)
    assert "shorter than" in fig.layout.annotations[0].text


def test_selection_from_points():
    assert selection_from_points([{"x": 12}, {"x": 9.0}, {"x": 15}]) == ResidueRange(9, 15)
    assert selection_from_points([]) is None


# --- Composition ---------------------------------------------------------------


def test_residue_orders():
    values = pd.Series({r: 0.0 for r in RESIDUE_ORDER} | {"L": 5.0, "A": 3.0})
    assert residue_order(values, "class") == list(RESIDUE_ORDER)
    assert residue_order(values, "alphabetical")[:3] == ["A", "C", "D"]
    assert residue_order(values, "frequency")[:2] == ["L", "A"]
    with pytest.raises(ValueError):
        residue_order(values, "random")


def test_composition_counts_vs_percent(records):
    seq = _sequence(records)
    counts = pd.Series(residue_counts(seq))
    percent = composition_bar_chart(counts, 110, LIGHT, as_percent=True)
    raw = composition_bar_chart(counts, 110, LIGHT, as_percent=False)
    leucine = lambda fig: dict(zip(fig.data[0].x, fig.data[0].y))["L"]  # noqa: E731
    assert leucine(raw) == 20
    assert leucine(percent) == pytest.approx(20 / 110 * 100)
    assert raw.layout.yaxis.title.text == "Count (residues)"


def test_composition_reference_in_bar_units(records):
    counts = pd.Series(residue_counts(_sequence(records)))
    reference = pd.Series({r: 0.05 for r in RESIDUE_ORDER})
    fig = composition_bar_chart(counts, 110, LIGHT, as_percent=False, reference=reference)
    ticks = next(t for t in fig.data if t.name == "Dataset reference")
    assert ticks.y[0] == pytest.approx(0.05 * 110)
    assert ticks.marker.color == LIGHT.ink


def test_dataset_composition_chart_two_weightings(records, features):
    comp = dataset_composition(build_count_table(records), features["length"])
    fig = dataset_composition_chart(comp, LIGHT)
    assert [t.name.split(" ")[0] for t in fig.data] == ["Residue-weighted", "Protein-weighted"]


# --- Distributions -------------------------------------------------------------


def test_default_bins_sturges():
    assert default_bins(8) == 4
    assert default_bins(1) == 1


def test_histogram_counts_and_reference_lines(features):
    fig = distribution_histogram(features["length"], LIGHT, xlabel="Length", bins=4)
    assert sum(fig.data[0].y) == 8
    rug = fig.data[1]
    assert len(rug.x) == 8 and rug.marker.color == LIGHT.ink_secondary
    labels = [a.text for a in fig.layout.annotations]
    assert any(t.startswith("Mean") for t in labels) and any(t.startswith("Median") for t in labels)


def test_histogram_identical_values_and_empty():
    assert sum(distribution_histogram(pd.Series([5.0, 5.0]), LIGHT, xlabel="x").data[0].y) == 2
    assert "No values" in distribution_histogram(pd.Series([], dtype=float), LIGHT, xlabel="x").layout.annotations[0].text


def test_linear_fit_recovers_known_slope():
    x = pd.Series([1.0, 2.0, 3.0, 4.0])
    slope, intercept, r2 = linear_fit(x, 110 * x + 18)
    assert (slope, intercept, r2) == pytest.approx((110, 18, 1))
    assert linear_fit(x[:2], x[:2]) is None


def test_mw_vs_length_fit_and_highlight(features):
    fig = mw_vs_length(features, LIGHT, highlight=INSULIN)
    assert any("Da per residue" in t.name for t in fig.data)
    selected = next(t for t in fig.data if t.name == "Selected protein")
    assert selected.x[0] == 110


# --- Comparison ----------------------------------------------------------------


def test_comparison_features_include_net_charge(features):
    columns = comparison_features(features)
    assert columns[0] == "length" and columns[-1] == "net_charge_ph7"
    assert len(columns) == 8
    assert feature_label("net_charge_ph7.4") == "Net charge at pH 7.4 (estimated)"


def test_scatter_customdata_maps_points_to_identifiers(features):
    fig = comparison_scatter(features, "isoelectric_point", "gravy", LIGHT, selected=INSULIN)
    ids = [row[0] for row in fig.data[0].customdata]
    assert set(ids) == set(features.index)
    assert any(t.name == "Selected protein" for t in fig.data)


def test_scatter_highlight_dims_other_points(features):
    fig = comparison_scatter(features, "length", "gravy", LIGHT, highlighted=[INSULIN])
    opacity = dict(zip((r[0] for r in fig.data[0].customdata), fig.data[0].marker.opacity))
    assert opacity[INSULIN] == 1.0 and min(opacity.values()) < 1.0


def test_group_labels_fold_extra_groups():
    # Three most common groups keep their names; the rest (and missing) become Other
    groups = group_labels(pd.Series(["a"] * 3 + ["b"] * 2 + ["c"] * 2 + ["d", None]))
    assert set(groups) == {"a", "b", "c", "Other"}
    assert groups.iloc[-2:].tolist() == ["Other", "Other"]


def test_zscore_is_mean_zero_unit_sd(features):
    z = zscore_table(features, ["length", "gravy"])
    assert z.mean().to_numpy() == pytest.approx([0, 0], abs=1e-12)
    assert z.std(ddof=1).to_numpy() == pytest.approx([1, 1])


def test_zscore_constant_column_is_zero():
    z = zscore_table(pd.DataFrame({"a": [3.0, 3.0, 3.0]}), ["a"])
    assert (z["a"] == 0).all()


def test_normalized_heatmap_keeps_raw_values_in_hover(features):
    fig = normalized_feature_heatmap(features, ["length", "gravy"], DARK)
    assert np.asarray(fig.data[0].customdata)[:, 0].tolist() == features["length"].tolist()
    assert fig.data[0].zmid == 0


# --- Composition heatmap -------------------------------------------------------


def test_order_rows(composition):
    assert order_rows(composition, "input") == list(composition.index)
    labels = pd.Series({i: i.split("|")[-1] for i in composition.index})
    assert order_rows(composition, "name", labels)[0] == "sp|P02769|ALBU_BOVIN"
    assert order_rows(composition, "L")[0] == INSULIN  # highest Leu frequency
    clustered = order_rows(composition, "cluster")
    assert sorted(clustered) == sorted(composition.index) and clustered != []
    with pytest.raises(ValueError):
        order_rows(composition, "zzz")


def test_heatmap_drops_invalid_rows_and_keeps_class_order():
    records = [FastaRecord("ok", "", "MKTAYIAKQR"), FastaRecord("bad", "", "MK-TAY")]
    fig = composition_heatmap(build_composition_table(records), LIGHT)
    assert list(fig.data[0].y) == ["ok"]
    assert list(fig.data[0].x) == list(RESIDUE_ORDER)
