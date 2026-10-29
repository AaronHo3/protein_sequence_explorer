# Validates parsed protein sequences: classifies residues, flags invalid characters,
# and detects sequences that are probably nucleic acid. 

from dataclasses import dataclass

STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")
AMBIGUOUS_RESIDUES = frozenset("BZJX")  # Asx, Glx, Xle, unknown
NONSTANDARD_RESIDUES = frozenset("UO")  # selenocysteine, pyrrolysine
NUCLEOTIDE_CHARS = frozenset("ACGTUN")
STOP = "*"

# Below this length, short peptides (e.g. "GAGA") are too easily mistaken for DNA
NUCLEOTIDE_CHECK_MIN_LENGTH = 20
# Real proteins are ~22% A+C+G+T; >= 90% is almost certainly nucleic acid
NUCLEOTIDE_FRACTION_THRESHOLD = 0.9

@dataclass(frozen=True)
class ValidationResult:
    """Outcome of validating one protein sequence."""
    is_valid: bool
    invalid_residues: tuple[tuple[int, str], ...]  # (1-based position, character)
    ambiguous_count: int
    nonstandard_count: int
    has_terminal_stop: bool
    warnings: tuple[str, ...]

def _looks_like_nucleotide(residues: str) -> bool:
    """True if a sequence (uppercase, terminal stop removed) is mostly nucleotide characters."""
    if len(residues) < NUCLEOTIDE_CHECK_MIN_LENGTH:
        return False
    nucleotide_count = sum(1 for char in residues if char in NUCLEOTIDE_CHARS)
    return nucleotide_count / len(residues) >= NUCLEOTIDE_FRACTION_THRESHOLD

def validate_sequence(sequence: str) -> ValidationResult:
    """Validate a protein sequence without modifying it.

    Comparison is case-insensitive; reported characters and positions refer to the
    original sequence.
    """
    if not sequence:
        return ValidationResult(False, (), 0, 0, False, ("Sequence is empty.",))

    upper = sequence.upper()
    has_terminal_stop = upper.endswith(STOP)
    # Residues to classify: the terminal stop is a translation artifact, not a residue
    residues = upper[:-1] if has_terminal_stop else upper

    invalid: list[tuple[int, str]] = []
    ambiguous_count = 0
    nonstandard_count = 0
    for position, char in enumerate(residues, start=1):
        if char in STANDARD_RESIDUES:
            continue
        if char in AMBIGUOUS_RESIDUES:
            ambiguous_count += 1
        elif char in NONSTANDARD_RESIDUES:
            nonstandard_count += 1
        else:
            invalid.append((position, sequence[position - 1]))

    warnings: list[str] = []
    if sequence != upper:
        warnings.append("Sequence contains lowercase residues.")
    if ambiguous_count:
        warnings.append(f"Sequence contains {ambiguous_count} ambiguous residue(s) (B/Z/J/X).")
    if nonstandard_count:
        warnings.append(f"Sequence contains {nonstandard_count} nonstandard residue(s) (U/O).")
    if has_terminal_stop:
        warnings.append("Sequence ends with a stop codon '*'.")
    if invalid:
        warnings.append(f"Sequence contains {len(invalid)} invalid character(s).")

    # A sequence with no determinable residues carries no usable information
    all_ambiguous = len(residues) > 0 and ambiguous_count == len(residues)
    if all_ambiguous:
        warnings.append("Sequence consists only of ambiguous residues.")
    if not residues:
        warnings.append("Sequence contains no residues besides a stop codon.")
    is_nucleotide = _looks_like_nucleotide(residues)
    if is_nucleotide:
        warnings.append("Sequence looks like nucleic acid (DNA/RNA), not protein.")

    return ValidationResult(
        is_valid=not (invalid or all_ambiguous or is_nucleotide or not residues),
        invalid_residues=tuple(invalid),
        ambiguous_count=ambiguous_count,
        nonstandard_count=nonstandard_count,
        has_terminal_stop=has_terminal_stop,
        warnings=tuple(warnings),
    )
