"""Tests for design tokens and colour helpers."""

import pytest

from protein_explorer.visualization.theme import (
    DARK,
    LIGHT,
    contrast_ratio,
    mix,
    palette_for,
    plotly_config,
    readable_text,
)


def test_palette_for_theme_type():
    assert palette_for("dark") is DARK
    assert palette_for("light") is LIGHT
    assert palette_for(None) is LIGHT


def test_mix_endpoints_and_midpoint():
    assert mix("#000000", "#ffffff", 0) == "#000000"
    assert mix("#000000", "#ffffff", 1) == "#ffffff"
    assert mix("#000000", "#ffffff", 0.5) == "#808080"


def test_contrast_ratio_black_white_is_21():
    assert contrast_ratio("#000000", "#ffffff") == pytest.approx(21.0)


@pytest.mark.parametrize("fill", ["#2a78d6", "#eda100", "#ecebe6", "#0d366b", "#f0a9c4"])
def test_readable_text_meets_aa_large_text_contrast(fill):
    assert contrast_ratio(fill, readable_text(fill, LIGHT)) >= 3.0


def test_dark_mode_is_its_own_palette():
    assert LIGHT.series != DARK.series
    assert LIGHT.surface != DARK.surface
    assert DARK.sequential[0] == "#0d366b"  # low end recedes toward the dark surface


def test_accent_distinct_from_data_series():
    for palette in (LIGHT, DARK):
        assert palette.accent not in palette.series


def test_plotly_config_exports_png_at_3x():
    config = plotly_config("chart")
    assert config["toImageButtonOptions"] == {"format": "png", "scale": 3, "filename": "chart"}
    assert config["displaylogo"] is False


def test_palettes_are_immutable():
    with pytest.raises(AttributeError):
        LIGHT.accent = "#000000"  # type: ignore[misc]
