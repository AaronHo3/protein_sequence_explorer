# Residue colouring schemes for the sequence viewer and legends. Category memberships
# reuse the residue sets in features.py, so colouring by "hydrophobic" highlights exactly
# the residues counted by hydrophobic_fraction.

from dataclasses import dataclass

from protein_explorer.features import (
    ACIDIC_RESIDUES,
    AROMATIC_RESIDUES,
    BASIC_RESIDUES,
    HYDROPHOBIC_RESIDUES,
)
from protein_explorer.residues import RESIDUE_CLASSES, RESIDUE_ORDER
from protein_explorer.visualization.theme import ChartPalette, mix

NONSTANDARD_LABEL = "Ambiguous or nonstandard"

@dataclass(frozen=True)
class Category:
    label: str
    residues: str  # one-letter codes in this category
    color: str
    note: str = ""

@dataclass(frozen=True)
class ColorScheme:
    key: str
    label: str
    description: str
    categories: tuple[Category, ...]
    nonstandard_color: str

    def color_map(self) -> dict[str, str]:
        """One-letter code -> fill colour for every standard residue."""
        return {r: c.color for c in self.categories for r in c.residues}

    def category_of(self, residue: str) -> str:
        upper = residue.upper()
        return next((c.label for c in self.categories if upper in c.residues), NONSTANDARD_LABEL)

def _standard_except(*groups: frozenset[str] | str) -> str:
    excluded = set().union(*(set(g) for g in groups))
    return "".join(r for r in RESIDUE_ORDER if r not in excluded)

def _ordered(residues: frozenset[str] | str) -> str:
    return "".join(r for r in RESIDUE_ORDER if r in residues)

def _identity(palette: ChartPalette) -> tuple[Category, ...]:
    """One colour per amino acid, built as tints/shades of its chemical-class hue so related
    residues read as a family. The one-letter code is always shown, so colour is never the
    only identifier."""
    categories = []
    for slot, (class_name, members) in enumerate(RESIDUE_CLASSES.items()):
        base = palette.series[slot]
        for i, residue in enumerate(members):
            # Spread from a darker shade to a lighter tint across the class
            t = (i / max(len(members) - 1, 1)) * 0.6 - 0.25
            color = mix(base, "#000000", -t) if t < 0 else mix(base, "#ffffff", t)
            categories.append(Category(residue, residue, color, class_name))
    return tuple(categories)

SCHEME_KEYS: tuple[str, ...] = (
    "chemical_class", "identity", "hydropathy", "charge", "polarity", "aromatic", "acid_base",
)

def color_scheme(key: str, palette: ChartPalette) -> ColorScheme:
    """Build a colouring scheme for a palette. Raises KeyError for unknown keys."""
    neutral = palette.neutral_fill
    builders = {
        "chemical_class": lambda: ColorScheme(
            "chemical_class", "Chemical class",
            "Five disjoint side-chain classes. Context-dependent residues are placed by their "
            "dominant property: Cys and Tyr can be weakly polar, His is basic but mostly "
            "neutral at pH 7.4.",
            tuple(Category(name, members, palette.series[i])
                  for i, (name, members) in enumerate(RESIDUE_CLASSES.items())),
            neutral,
        ),
        "identity": lambda: ColorScheme(
            "identity", "Amino acid identity",
            "Each amino acid has its own colour, shaded within its chemical-class family.",
            _identity(palette), neutral,
        ),
        "hydropathy": lambda: ColorScheme(
            "hydropathy", "Hydrophobic vs hydrophilic",
            "Hydrophobic = A, V, I, L, M, F, W (the set used for hydrophobic fraction). "
            "Gly, Cys, Pro, and Tyr are borderline and counted as hydrophilic here.",
            (Category("Hydrophobic", _ordered(HYDROPHOBIC_RESIDUES), palette.hydrophobic),
             Category("Hydrophilic", _standard_except(HYDROPHOBIC_RESIDUES), palette.hydrophilic)),
            neutral,
        ),
        "charge": lambda: ColorScheme(
            "charge", "Positive vs negative charge",
            "Side-chain charge near pH 7.4. His (pKa about 6.5) is only about 10% protonated, "
            "so it is shown separately. Termini are not shown.",
            (Category("Positive (K, R)", "KR", palette.charge_positive),
             Category("Partly positive (H)", "H", palette.charge_partial),
             Category("Negative (D, E)", "DE", palette.charge_negative),
             Category("Uncharged", _standard_except("KRHDE"), neutral)),
            neutral,
        ),
        "polarity": lambda: ColorScheme(
            "polarity", "Polar vs nonpolar",
            "Nonpolar = G, A, V, L, I, M, P, F, W. Polar includes charged side chains and "
            "Tyr (hydroxyl). Classifications of Cys, Trp, and Tyr vary between textbooks.",
            (Category("Nonpolar", "GAVLIMPFW", palette.hydrophobic),
             Category("Polar", "STCNQYDEKRH", palette.hydrophilic)),
            neutral,
        ),
        "aromatic": lambda: ColorScheme(
            "aromatic", "Aromatic residues",
            "Aromatic = F, W, Y (the set used for aromaticity). His has an aromatic imidazole "
            "ring but is not counted, matching Biopython.",
            (Category("Aromatic (F, W, Y)", _ordered(AROMATIC_RESIDUES), palette.aromatic),
             Category("His (not counted)", "H", palette.aromatic_partial),
             Category("Other", _standard_except(AROMATIC_RESIDUES, "H"), neutral)),
            neutral,
        ),
        "acid_base": lambda: ColorScheme(
            "acid_base", "Acidic, basic, neutral",
            "Acidic = D, E. Basic = K, R, H (the sets used for acidic and basic counts).",
            (Category("Acidic", _ordered(ACIDIC_RESIDUES), palette.charge_negative),
             Category("Basic", _ordered(BASIC_RESIDUES), palette.charge_positive),
             Category("Neutral", _standard_except(ACIDIC_RESIDUES, BASIC_RESIDUES), neutral)),
            neutral,
        ),
    }
    return builders[key]()

def scheme_labels() -> dict[str, str]:
    """Scheme key -> display label, in menu order (labels do not depend on the palette)."""
    from protein_explorer.visualization.theme import LIGHT
    return {key: color_scheme(key, LIGHT).label for key in SCHEME_KEYS}
