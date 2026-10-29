# Reference data for individual residues: names, codes, chemical class, hydropathy, and
# side-chain ionisation. Values come from features.py so every view uses the same
# definitions as the calculated features. No UI code here.

from dataclasses import dataclass

from protein_explorer.features import (
    KYTE_DOOLITTLE,
    NEGATIVE_PKA,
    POSITIVE_PKA,
)

THREE_LETTER: dict[str, str] = {
    "A": "Ala", "C": "Cys", "D": "Asp", "E": "Glu", "F": "Phe", "G": "Gly", "H": "His",
    "I": "Ile", "K": "Lys", "L": "Leu", "M": "Met", "N": "Asn", "P": "Pro", "Q": "Gln",
    "R": "Arg", "S": "Ser", "T": "Thr", "V": "Val", "W": "Trp", "Y": "Tyr",
    # Ambiguous and nonstandard codes (IUPAC)
    "B": "Asx", "Z": "Glx", "J": "Xle", "X": "Xaa", "U": "Sec", "O": "Pyl",
}

NAMES: dict[str, str] = {
    "A": "Alanine", "C": "Cysteine", "D": "Aspartate", "E": "Glutamate",
    "F": "Phenylalanine", "G": "Glycine", "H": "Histidine", "I": "Isoleucine",
    "K": "Lysine", "L": "Leucine", "M": "Methionine", "N": "Asparagine",
    "P": "Proline", "Q": "Glutamine", "R": "Arginine", "S": "Serine",
    "T": "Threonine", "V": "Valine", "W": "Tryptophan", "Y": "Tyrosine",
    "B": "Aspartate or asparagine", "Z": "Glutamate or glutamine",
    "J": "Leucine or isoleucine", "X": "Unknown residue",
    "U": "Selenocysteine", "O": "Pyrrolysine",
}

# Disjoint chemical classes. Display order is also the default axis order, so each
# class forms one contiguous block in charts.
RESIDUE_CLASSES: dict[str, str] = {
    "Nonpolar": "GAVLIMP",
    "Aromatic": "FWY",
    "Polar": "STCNQ",
    "Acidic": "DE",
    "Basic": "KRH",
}
RESIDUE_ORDER: tuple[str, ...] = tuple("".join(RESIDUE_CLASSES.values()))
CHEMICAL_CLASS: dict[str, str] = {r: name for name, members in RESIDUE_CLASSES.items() for r in members}

@dataclass(frozen=True)
class ResidueInfo:
    """Everything the interface shows about one residue type."""
    code: str
    three_letter: str
    name: str
    chemical_class: str
    hydropathy: float | None  # Kyte-Doolittle; None for ambiguous/nonstandard codes
    side_chain_pka: float | None
    charge_sign: int  # +1 basic, -1 acidic, 0 not ionisable in the pKa set
    is_standard: bool

def residue_info(code: str) -> ResidueInfo:
    """Reference data for a one-letter code (case-insensitive). Unknown symbols are reported
    as nonstandard rather than raising, so the viewer can show any sequence."""
    upper = code.upper()
    pka, sign = _side_chain_pka(upper)
    return ResidueInfo(
        code=upper,
        three_letter=THREE_LETTER.get(upper, "?"),
        name=NAMES.get(upper, f"Unrecognised symbol '{code}'"),
        chemical_class=CHEMICAL_CLASS.get(upper, "Ambiguous or nonstandard"),
        hydropathy=KYTE_DOOLITTLE.get(upper),
        side_chain_pka=pka,
        charge_sign=sign,
        is_standard=upper in CHEMICAL_CLASS,
    )

def _side_chain_pka(code: str) -> tuple[float | None, int]:
    if code in POSITIVE_PKA:
        return POSITIVE_PKA[code], 1
    if code in NEGATIVE_PKA:
        return NEGATIVE_PKA[code], -1
    return None, 0

def side_chain_charge(code: str, ph: float) -> float:
    """Average side-chain charge at a pH (Henderson-Hasselbalch, EMBOSS pKa values).
    Returns 0 for side chains that are not ionisable in the pKa set."""
    pka, sign = _side_chain_pka(code.upper())
    if pka is None:
        return 0.0
    if sign > 0:
        return 1 / (1 + 10 ** (ph - pka))
    return -1 / (1 + 10 ** (pka - ph))
