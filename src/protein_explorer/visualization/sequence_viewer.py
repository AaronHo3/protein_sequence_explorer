# Data model for the interactive sequence viewer. Builds a JSON-serialisable payload that
# the browser component renders; no Streamlit here, so it is testable and reusable.
#
# Extension points (designed in, used as noted):
#   tracks      per-residue numeric rows under the sequence. Used now for hydropathy;
#               later conservation scores or experimental measurements.
#   annotations labelled residue ranges drawn above the sequence (domains, mutations).
#   positions   column -> residue number mapping, so gapped alignment rows can reuse the
#               viewer with correct numbering (None marks a gap column).

import math
from collections.abc import Sequence
from dataclasses import dataclass, field

from protein_explorer.features import hydropathy_profile
from protein_explorer.residues import residue_info
from protein_explorer.visualization.color_schemes import ColorScheme
from protein_explorer.visualization.theme import ChartPalette, readable_text

BLOCK_SIZE = 10
LINE_LENGTH_OPTIONS: tuple[int, ...] = (0, 30, 50, 60, 80, 100)  # 0 = fit to width
GAP = "-"

@dataclass(frozen=True)
class ResidueRange:
    """Inclusive 1-based residue range, the standard biological numbering convention."""
    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 1 or self.end < self.start:
            raise ValueError(f"Invalid residue range {self.start}-{self.end}.")

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def label(self) -> str:
        return str(self.start) if self.start == self.end else f"{self.start}–{self.end}"

    @classmethod
    def from_positions(cls, positions: Sequence[int]) -> "ResidueRange | None":
        valid = [int(p) for p in positions]
        return cls(min(valid), max(valid)) if valid else None

@dataclass(frozen=True)
class ResidueTrack:
    """A numeric value per residue (None where undefined), drawn as a colour strip."""
    name: str
    values: tuple[float | None, ...]
    low: float
    high: float
    low_color: str
    mid_color: str
    high_color: str
    description: str = ""

@dataclass(frozen=True)
class ResidueAnnotation:
    """A labelled residue range, e.g. a domain or a mutation."""
    start: int
    end: int
    label: str
    color: str

@dataclass(frozen=True)
class ViewerOptions:
    line_length: int = 0  # 0: as many 10-residue blocks as fit the width
    show_tracks: bool = True
    tracks: tuple[ResidueTrack, ...] = field(default_factory=tuple)
    annotations: tuple[ResidueAnnotation, ...] = field(default_factory=tuple)

def hydropathy_track(sequence: str, window: int, palette: ChartPalette) -> ResidueTrack:
    """Window-mean hydropathy aligned to residue positions. Positions within window//2 of
    either end have no full window, so they are None rather than extrapolated."""
    residues = sequence.removesuffix("*")
    values: list[float | None] = [None] * len(residues)
    if len(residues) >= window:
        for position, value in hydropathy_profile(residues, window):
            values[position - 1] = None if math.isnan(value) else round(value, 3)
    return ResidueTrack(
        name=f"Hydropathy, window {window}",
        values=tuple(values),
        low=-4.5, high=4.5,
        low_color=palette.hydrophilic, mid_color=palette.diverging_mid,
        high_color=palette.hydrophobic,
        description="Kyte–Doolittle sliding-window mean centred on each residue.",
    )

def _residue_lookup(sequence: str, scheme: ColorScheme, palette: ChartPalette) -> dict:
    """Per distinct symbol: fill, text colour, and the reference data shown on hover."""
    colors = scheme.color_map()
    lookup = {}
    for symbol in sorted(set(sequence)):
        info = residue_info(symbol)
        if symbol == GAP:
            fill, name = "transparent", "Gap"
        else:
            fill, name = colors.get(info.code, scheme.nonstandard_color), info.name
        lookup[symbol] = {
            "fill": fill,
            "text": palette.muted if symbol == GAP else readable_text(
                fill if fill != "transparent" else palette.surface, palette),
            "three": info.three_letter,
            "name": name,
            "category": scheme.category_of(info.code),
            "chemicalClass": info.chemical_class,
            "hydropathy": info.hydropathy,
            "standard": info.is_standard,
        }
    return lookup

def build_viewer_data(
    sequence: str,
    scheme: ColorScheme,
    palette: ChartPalette,
    options: ViewerOptions = ViewerOptions(),
    selection: ResidueRange | None = None,
    positions: Sequence[int | None] | None = None,
) -> dict:
    """JSON payload for the browser viewer. The sequence is sent once as a string with a
    small per-symbol lookup, so payload size grows with length, not with residue metadata."""
    residues = sequence.removesuffix("*")
    if positions is not None and len(positions) != len(residues):
        raise ValueError("positions must have one entry per sequence column.")
    return {
        "sequence": residues,
        "positions": list(positions) if positions is not None else None,
        "lookup": _residue_lookup(residues, scheme, palette),
        "legend": [{"label": c.label, "color": c.color, "residues": c.residues}
                   for c in scheme.categories],
        "blockSize": BLOCK_SIZE,
        "lineLength": options.line_length or None,
        "selection": [selection.start, selection.end] if selection else None,
        "tracks": [
            {"name": t.name, "values": list(t.values), "low": t.low, "high": t.high,
             "lowColor": t.low_color, "midColor": t.mid_color, "highColor": t.high_color,
             "description": t.description}
            for t in (options.tracks if options.show_tracks else ())
        ],
        "annotations": [
            {"start": a.start, "end": a.end, "label": a.label, "color": a.color}
            for a in options.annotations
        ],
        "colors": {"selection": palette.selection, "muted": palette.muted,
                   "surface": palette.surface, "raised": palette.raised},
    }
