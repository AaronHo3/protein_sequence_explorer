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

# Average residue masses in Da (amino acid minus one water), derived from
# Biopython's IUPACData.protein_weights. Includes U (Sec) and O (Pyl).
RESIDUE_MASSES: dict[str, float] = {
    "A": 71.0779, "C": 103.1429, "D": 115.0874, "E": 129.1140, "F": 147.1738,
    "G": 57.0513, "H": 137.1393, "I": 113.1576, "K": 128.1723, "L": 113.1576,
    "M": 131.1960, "N": 114.1026, "O": 237.2981, "P": 97.1152, "Q": 128.1292,
    "R": 156.1857, "S": 87.0773, "T": 101.1039, "U": 150.0379, "V": 99.1310,
    "W": 186.2099, "Y": 163.1732,
}
WATER_MASS = 18.01528  # average mass of H2O, added once for the chain termini

# Side-chain and terminal pKa values from EMBOSS (iep). U (selenocysteine) is not in
# EMBOSS; 5.2 is a literature estimate. Real pKa values shift with local environment.
POSITIVE_PKA: dict[str, float] = {"Nterm": 8.6, "K": 10.8, "R": 12.5, "H": 6.5}
NEGATIVE_PKA: dict[str, float] = {
    "Cterm": 3.6, "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1, "U": 5.2,
}
PKA_SET = "EMBOSS (U: 5.2, literature estimate)"  # recorded in run metadata
PI_TOLERANCE = 1e-4

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

def molecular_weight(sequence: str) -> float:
    """Average molecular weight in Da: sum of residue masses plus one water.

    Assumes no post-translational modifications and reduced cysteines (each
    disulfide bond would lower the mass by ~2.016 Da). Raises ValueError for
    residues without a defined mass (B, Z, J, X).
    """
    residues = _residues(sequence)
    unknown = sorted(set(residues) - RESIDUE_MASSES.keys())
    if unknown:
        raise ValueError(f"Cannot compute molecular weight: no mass for {', '.join(unknown)}.")
    return sum(RESIDUE_MASSES[r] for r in residues) + WATER_MASS

def _ionisable_group_counts(sequence: str) -> Counter[str]:
    """Counts of each ionisable side chain, plus one N- and one C-terminus."""
    counts = Counter(_residues(sequence))
    counts["Nterm"] = 1
    counts["Cterm"] = 1
    return counts

def net_charge(sequence: str, ph: float = 7.0) -> float:
    """Estimated net charge at a given pH using the Henderson-Hasselbalch equation.

    Treats every ionisable group as independent and fully solvent-exposed, so it
    ignores structural pKa shifts and counts disulfide-bonded cysteines as ionisable.
    """
    counts = _ionisable_group_counts(sequence)
    positive = sum(
        counts[group] / (1 + 10 ** (ph - pka)) for group, pka in POSITIVE_PKA.items()
    )
    negative = sum(
        counts[group] / (1 + 10 ** (pka - ph)) for group, pka in NEGATIVE_PKA.items()
    )
    return positive - negative

def isoelectric_point(sequence: str) -> float:
    """Estimated pI: the pH where net_charge is zero, found by bisection on [0, 14].

    Net charge decreases monotonically with pH, so bisection always converges.
    """
    _residues(sequence)  # validate non-empty before searching
    low, high = 0.0, 14.0
    while high - low > PI_TOLERANCE:
        mid = (low + high) / 2
        if net_charge(sequence, mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2

def hydropathy_profile(sequence: str, window: int = 9) -> list[tuple[int, float]]:
    """Sliding-window Kyte-Doolittle hydropathy.

    Returns (1-based centre position, mean hydropathy) for each full window. Residues
    without a hydropathy value (B/Z/J/X/U/O) are skipped within a window; a window
    with no standard residues gives NaN. Window must be odd and no longer than the
    sequence. Kyte & Doolittle suggest 9 for surface regions, 19 for transmembrane
    helices (mean > 1.6).
    """
    residues = _residues(sequence)
    if window < 1 or window % 2 == 0:
        raise ValueError(f"Window must be a positive odd number, got {window}.")
    if window > len(residues):
        raise ValueError(f"Window {window} is longer than the sequence ({len(residues)}).")

    half = window // 2
    profile: list[tuple[int, float]] = []
    for start in range(len(residues) - window + 1):
        values = [KYTE_DOOLITTLE[r] for r in residues[start:start + window] if r in KYTE_DOOLITTLE]
        mean = sum(values) / len(values) if values else float("nan")
        profile.append((start + half + 1, mean))
    return profile
