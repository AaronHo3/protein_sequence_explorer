# Sequence-level feature calculations. Functions assume the sequence has already
# passed validation (see validation.py) and do not re-validate it.

from collections import Counter

from protein_explorer.validation import STANDARD_RESIDUES, STOP

# Fixed column order for feature vectors: same keys, same order, every run
AMINO_ACIDS: tuple[str, ...] = tuple(sorted(STANDARD_RESIDUES))

# Residue groupings. These are conventions, not physical constants; other sources differ.
# Aliphatic + large nonpolar side chains; excludes G (no side chain), C, P, Y
HYDROPHOBIC_RESIDUES = frozenset("AVILMFW")
ACIDIC_RESIDUES = frozenset("DE")
# H included: its side chain is basic, though mostly neutral at pH 7.4 (pKa ~6.0)
BASIC_RESIDUES = frozenset("KRH")
# Charged at physiological pH; H excluded for the reason above
CHARGED_RESIDUES = frozenset("DEKR")
# Same definition as Biopython's aromaticity (Lobry & Gautier, 1994); H excluded
AROMATIC_RESIDUES = frozenset("FWY")

# Kyte & Doolittle (1982) hydropathy index. Positive = hydrophobic.
KYTE_DOOLITTLE: dict[str, float] = {
    "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5, "F": 2.8,
    "G": -0.4, "H": -3.2, "I": 4.5, "K": -3.9, "L": 3.8,
    "M": 1.9, "N": -3.5, "P": -1.6, "Q": -3.5, "R": -4.5,
    "S": -0.8, "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3,
}

def _residues(sequence: str) -> str:
    """Uppercase sequence with a terminal stop removed. Raises ValueError if nothing remains."""
    residues = sequence.upper().removesuffix(STOP)
    if not residues:
        raise ValueError("Cannot compute features of an empty sequence.")
    return residues

def sequence_length(sequence: str) -> int:
    """Number of residues, excluding a terminal stop '*'."""
    return len(_residues(sequence))

def residue_counts(sequence: str) -> dict[str, int]:
    """Counts for all 20 standard residues (zeros included), case-insensitive,
    in alphabetical order."""
    counts = Counter(_residues(sequence))
    return {aa: counts[aa] for aa in AMINO_ACIDS}

def residue_frequencies(sequence: str) -> dict[str, float]:
    """Count / sequence_length for each of the 20 standard residues.

    Ambiguous and nonstandard residues count toward the length, so frequencies
    sum to less than 1.0 when they are present.
    """
    length = sequence_length(sequence)
    return {aa: count / length for aa, count in residue_counts(sequence).items()}

def unique_residue_types(sequence: str) -> int:
    """Number of distinct standard residue types present (0–20)."""
    return sum(1 for count in residue_counts(sequence).values() if count > 0)

def _group_count(sequence: str, group: frozenset[str]) -> int:
    """Number of residues in the sequence belonging to a residue group."""
    return sum(1 for residue in _residues(sequence) if residue in group)

def hydrophobic_fraction(sequence: str) -> float:
    """Fraction of residues that are hydrophobic (A, V, I, L, M, F, W)."""
    return _group_count(sequence, HYDROPHOBIC_RESIDUES) / sequence_length(sequence)

def charged_fraction(sequence: str) -> float:
    """Fraction of residues charged at physiological pH (D, E, K, R)."""
    return _group_count(sequence, CHARGED_RESIDUES) / sequence_length(sequence)

def aromatic_fraction(sequence: str) -> float:
    """Fraction of residues that are aromatic (F, W, Y)."""
    return _group_count(sequence, AROMATIC_RESIDUES) / sequence_length(sequence)

def acidic_count(sequence: str) -> int:
    """Number of acidic residues (D, E)."""
    return _group_count(sequence, ACIDIC_RESIDUES)

def basic_count(sequence: str) -> int:
    """Number of basic residues (K, R, H)."""
    return _group_count(sequence, BASIC_RESIDUES)

def gravy(sequence: str) -> float:
    """Grand average of hydropathy: mean Kyte-Doolittle value over standard residues.

    Ambiguous and nonstandard residues have no hydropathy value and are excluded
    from both the sum and the denominator. Raises ValueError if no standard
    residues remain.
    """
    values = [KYTE_DOOLITTLE[r] for r in _residues(sequence) if r in KYTE_DOOLITTLE]
    if not values:
        raise ValueError("Cannot compute GRAVY: sequence has no standard residues.")
    return sum(values) / len(values)
