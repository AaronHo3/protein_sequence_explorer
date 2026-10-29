# Design tokens shared by every view: chart palettes (light and dark), UI accents, and
# chart styling. Categorical slots are validated for colour-vision deficiency separately
# for each surface; dark is its own set of steps, not an automatic inversion.
# Kept free of plotly/streamlit imports so static (matplotlib) figures can use it too.

from dataclasses import dataclass
from typing import Any

FONT_FAMILY = "IBM Plex Sans, system-ui, sans-serif"
MONO_FAMILY = "IBM Plex Mono, ui-monospace, monospace"

@dataclass(frozen=True)
class ChartPalette:
    """Colour roles for one surface. Data colours never double as UI chrome."""
    name: str
    surface: str
    raised: str  # card and panel background
    ink: str
    ink_secondary: str
    muted: str
    gridline: str
    baseline: str
    series: tuple[str, ...]  # categorical slots in fixed (validated) order
    sequential: tuple[str, ...]  # low -> high; the low end recedes toward the surface
    diverging_low: str  # values below zero
    diverging_high: str  # values above zero
    diverging_mid: str  # neutral midpoint (grey, never a hue)
    hydrophobic: str
    hydrophilic: str
    charge_positive: str
    charge_negative: str
    charge_partial: str  # His: only partly protonated near pH 7
    aromatic: str
    aromatic_partial: str  # His: aromatic ring, but not counted in aromaticity
    neutral_fill: str  # residues not highlighted by the current scheme
    accent: str  # interactive UI elements (indigo), distinct from data colours
    accent_soft: str
    selection: str  # outline for selected residues and regions

LIGHT = ChartPalette(
    name="light",
    surface="#fcfcfb",
    raised="#ffffff",
    ink="#0b0b0b",
    ink_secondary="#52514e",
    muted="#898781",
    gridline="#e1e0d9",
    baseline="#c3c2b7",
    series=("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"),
    sequential=("#cde2fb", "#86b6ef", "#3987e5", "#256abf", "#184f95", "#0d366b"),
    diverging_low="#2a78d6",
    diverging_high="#e34948",
    diverging_mid="#f0efec",
    hydrophobic="#eb6834",
    hydrophilic="#2a78d6",
    charge_positive="#2a78d6",
    charge_negative="#e34948",
    charge_partial="#1baf7a",
    aromatic="#6250d6",
    aromatic_partial="#b4a9f0",
    neutral_fill="#ecebe6",
    accent="#4c55c9",
    accent_soft="#e8e9fb",
    selection="#4c55c9",
)

DARK = ChartPalette(
    name="dark",
    surface="#1a1a19",
    raised="#222220",
    ink="#f0efec",
    ink_secondary="#c3c2b7",
    muted="#898781",
    gridline="#2c2c2a",
    baseline="#383835",
    series=("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181"),
    sequential=("#0d366b", "#184f95", "#256abf", "#3987e5", "#6da7ec", "#b7d3f6"),
    diverging_low="#3987e5",
    diverging_high="#e66767",
    diverging_mid="#383835",
    hydrophobic="#d95926",
    hydrophilic="#3987e5",
    charge_positive="#3987e5",
    charge_negative="#e66767",
    charge_partial="#199e70",
    aromatic="#9085e9",
    aromatic_partial="#5d52c4",
    neutral_fill="#2f2f2c",
    accent="#8c93ec",
    accent_soft="#2b2d4a",
    selection="#a9aff5",
)

def palette_for(theme_type: str | None) -> ChartPalette:
    """Palette for Streamlit's active theme ('light', 'dark', or unknown -> light)."""
    return DARK if theme_type == "dark" else LIGHT

def _rgb(hex_color: str) -> tuple[int, int, int]:
    value = hex_color.lstrip("#")
    return int(value[0:2], 16), int(value[2:4], 16), int(value[4:6], 16)

def mix(color: str, other: str, amount: float) -> str:
    """Blend two hex colours: amount 0 -> color, 1 -> other."""
    a, b = _rgb(color), _rgb(other)
    blended = (round(x + (y - x) * amount) for x, y in zip(a, b, strict=True))
    return "#" + "".join(f"{c:02x}" for c in blended)

def relative_luminance(hex_color: str) -> float:
    """WCAG relative luminance (0 = black, 1 = white)."""
    def channel(c: int) -> float:
        s = c / 255
        return s / 12.92 if s <= 0.03928 else ((s + 0.055) / 1.055) ** 2.4
    r, g, b = (channel(c) for c in _rgb(hex_color))
    return 0.2126 * r + 0.7152 * g + 0.0722 * b

def contrast_ratio(a: str, b: str) -> float:
    high, low = sorted((relative_luminance(a), relative_luminance(b)), reverse=True)
    return (high + 0.05) / (low + 0.05)

def readable_text(fill: str, palette: ChartPalette) -> str:
    """Near-black or white label, whichever contrasts more with the fill."""
    dark, light = "#0b0b0b", "#ffffff"
    return dark if contrast_ratio(fill, dark) >= contrast_ratio(fill, light) else light

def style_figure(
    fig: Any, palette: ChartPalette, xlabel: str = "", ylabel: str = "", height: int | None = None
) -> Any:
    """Shared Plotly look: transparent background, recessive grid, palette ink, readable hover.
    Typed as Any so this module does not import plotly."""
    axis = {
        "gridcolor": palette.gridline, "linecolor": palette.baseline, "zeroline": False,
        "automargin": True, "tickfont": {"color": palette.ink_secondary, "size": 12},
        "title_font": {"color": palette.ink_secondary, "size": 13},
    }
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        font={"family": FONT_FAMILY, "color": palette.ink, "size": 13},
        margin={"l": 56, "r": 16, "t": 64, "b": 44},  # top room: legend below the toolbar
        hoverlabel={"bgcolor": palette.raised, "bordercolor": palette.baseline,
                    "font": {"family": FONT_FAMILY, "color": palette.ink, "size": 12}},
        legend={"orientation": "h", "yanchor": "bottom", "y": 1.0, "x": 0,
                "font": {"color": palette.ink_secondary, "size": 12}},
        modebar={"bgcolor": "rgba(0,0,0,0)", "color": palette.muted,
                 "activecolor": palette.accent},
    )
    if height is not None:
        fig.update_layout(height=height)
    fig.update_xaxes(title_text=xlabel, **axis)
    fig.update_yaxes(title_text=ylabel, **axis)
    return fig

def plotly_config(filename: str) -> dict[str, object]:
    """Toolbar settings: PNG export at 3x (slides, papers), no logo, no clutter."""
    return {
        "displaylogo": False,
        "modeBarButtonsToRemove": ["lasso2d", "autoScale2d", "toggleSpikelines"],
        "toImageButtonOptions": {"format": "png", "scale": 3, "filename": filename},
    }
