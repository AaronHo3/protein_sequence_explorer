# Command-line entry point: `python -m protein_explorer summarize <file.fasta>`.

import argparse
import sys

from protein_explorer.fasta import FastaFormatError, FastaRecord, read_fasta
from protein_explorer.validation import ValidationResult, validate_sequence

ID_WIDTH = 24

def format_summary(records: list[FastaRecord], results: list[ValidationResult]) -> str:
    """Build a plain-text table with one row per record, followed by totals."""
    header = f"{'Identifier':<{ID_WIDTH}} {'Length':>6}  {'Valid':<5}  Warnings"
    lines = [header, "-" * len(header)]
    for record, result in zip(records, results, strict=True):
        warnings = "; ".join(result.warnings) or "-"
        valid = "yes" if result.is_valid else "NO"
        lines.append(
            f"{record.identifier:<{ID_WIDTH}} {len(record.sequence):>6}  {valid:<5}  {warnings}"
        )

    lengths = [len(record.sequence) for record in records]
    n_valid = sum(result.is_valid for result in results)
    lines.append("")
    lines.append(f"Records: {len(records)}  Valid: {n_valid}  Invalid: {len(records) - n_valid}")
    lines.append(f"Length range: {min(lengths)}-{max(lengths)} residues")
    return "\n".join(lines)

def summarize(path: str) -> int:
    """Parse and validate a FASTA file, print the summary, and return an exit code."""
    try:
        records = read_fasta(path)
    except FileNotFoundError:
        print(f"Error: file not found: {path}", file=sys.stderr)
        return 1
    except FastaFormatError as error:
        print(f"Error: {path} is not valid FASTA. {error}", file=sys.stderr)
        return 1

    results = [validate_sequence(record.sequence) for record in records]
    print(format_summary(records, results))
    return 0

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="protein_explorer")
    subcommands = parser.add_subparsers(dest="command", required=True)
    summarize_parser = subcommands.add_parser(
        "summarize", help="Print a validation summary for a FASTA file."
    )
    summarize_parser.add_argument("fasta", help="Path to a protein FASTA file.")

    args = parser.parse_args(argv)
    if args.command == "summarize":
        return summarize(args.fasta)
    return 1

if __name__ == "__main__":
    sys.exit(main())
