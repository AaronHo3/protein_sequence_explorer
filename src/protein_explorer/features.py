# Sequence-level feature calculations. Functions assume the sequence has already
# passed validation (see validation.py) and do not re-validate it.

from collections import Counter

from protein_explorer.validation import STANDARD_RESIDUES, STOP

# Fixed column order for feature vectors: same keys, same order, every run
AMINO_ACIDS: tuple[str, ...] = tuple(sorted(STANDARD_RESIDUES))

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
