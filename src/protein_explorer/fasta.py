# Manual FASTA parser: reads records and raises FastaFormatError on malforemd input.

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

class FastaFormatError(ValueError):
    """Raised when input is not valid FASTA format. The message must include the line number."""

@dataclass(frozen=True)
class FastaRecord:
    """A single FASTA record."""
    identifier: str # first token after '>'
    description: str #rest of the header
    sequence: str # concatenated sequence lines

def _finish_record(
    identifier: str,
    description: str,
    sequence_lines: list[str],
    header_line_number: int,
    seen_identifiers: set[str],
) -> FastaRecord:
    """Validate a completed record and build it. Single place for record-level checks."""
    if not sequence_lines:
        raise FastaFormatError(
            f"Line {header_line_number}: Record '{identifier}' has no sequence."
        )
    if identifier in seen_identifiers:
        raise FastaFormatError(
            f"Line {header_line_number}: Duplicate identifier '{identifier}'."
        )
    seen_identifiers.add(identifier)
    return FastaRecord(identifier, description, ''.join(sequence_lines))

def parse_fasta_lines(lines: Iterable[str]) -> list[FastaRecord]:
    """Parse FASTA lines into a list of FastaRecord objects."""
    records: list[FastaRecord] = []
    seen_identifiers: set[str] = set()
    identifier: str | None = None
    description = ''
    header_line_number = 0
    sequence_lines: list[str] = []

    for line_number, line in enumerate(lines, start=1):
        line = line.strip()
        if not line:
            continue  # skip empty lines
        if line.startswith('>'):
            if identifier is not None:
                # Save the previous record
                records.append(_finish_record(
                    identifier, description, sequence_lines, header_line_number, seen_identifiers
                ))
            # Start a new record
            parts = line[1:].split(maxsplit=1)
            if not parts:
                raise FastaFormatError(f"Line {line_number}: Header has no identifier.")
            identifier = parts[0]
            description = parts[1] if len(parts) > 1 else ''
            header_line_number = line_number
            sequence_lines = []
        else:
            if identifier is None:
                raise FastaFormatError(f"Line {line_number}: Sequence data found before any header.")
            # Whitespace is formatting, not sequence data; residues themselves are left untouched
            sequence_lines.append(''.join(line.split()))

    # Add the last record if it exists
    if identifier is None:
        raise FastaFormatError("Line 0: No FASTA records found (input is empty).")
    records.append(_finish_record(
        identifier, description, sequence_lines, header_line_number, seen_identifiers
    ))

    return records

def read_fasta(path: str | Path) -> list[FastaRecord]:
    """Read a FASTA file and return a list of FastaRecord objects.

    Raises FileNotFoundError if the path does not exist (from open()).
    """
    with Path(path).open('r', encoding='utf-8') as f:
        return parse_fasta_lines(f)
