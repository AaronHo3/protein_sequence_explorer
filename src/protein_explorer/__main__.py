# Command-line entry point:
#   python -m protein_explorer summarize <file.fasta>
#   python -m protein_explorer analyze <file.fasta> --output results/
#   python -m protein_explorer plot <file.fasta> --output results/figures/

import argparse
import sys
from pathlib import Path

from protein_explorer.dataset import (
    build_composition_table,
    build_feature_table,
    summarize_features,
)
from protein_explorer.export import build_run_info, write_results
from protein_explorer.fasta import FastaFormatError, FastaRecord, read_fasta
from protein_explorer.validation import ValidationResult, validate_sequence

ID_WIDTH = 24
DEFAULT_OUTPUT_DIR = "results"
DEFAULT_PH = 7.0
DEFAULT_FIGURE_DIR = "results/figures"
DEFAULT_WINDOW = 9

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

def _load_records(path: str) -> list[FastaRecord] | None:
    """Read a FASTA file, printing a user-facing error and returning None on failure."""
    try:
        return read_fasta(path)
    except FileNotFoundError:
        print(f"Error: file not found: {path}", file=sys.stderr)
    except FastaFormatError as error:
        print(f"Error: {path} is not valid FASTA. {error}", file=sys.stderr)
    return None

def summarize(path: str) -> int:
    """Parse and validate a FASTA file, print the summary, and return an exit code."""
    records = _load_records(path)
    if records is None:
        return 1
    results = [validate_sequence(record.sequence) for record in records]
    print(format_summary(records, results))
    return 0

def analyze(path: str, output_dir: str, ph: float) -> int:
    """Compute feature tables for a FASTA file and write them to output_dir."""
    records = _load_records(path)
    if records is None:
        return 1
    feature_table = build_feature_table(records, ph=ph)
    composition_table = build_composition_table(records)
    summary_table = summarize_features(feature_table)
    run_info = build_run_info(Path(path), feature_table, ph)

    written = write_results(
        Path(output_dir), feature_table, composition_table, summary_table, run_info
    )
    n_valid = run_info["valid_records"]
    print(f"Analyzed {len(records)} records ({n_valid} valid, {len(records) - n_valid} invalid).")
    for file_path in written:
        print(f"  wrote {file_path}")
    return 0

def plot(path: str, output_dir: str, window: int, fmt: str) -> int:
    """Generate all figures for a FASTA file and write them to output_dir."""
    # Imported here so summarize/analyze don't pay matplotlib's import cost
    from protein_explorer.plots import write_all_figures

    if window < 1 or window % 2 == 0:
        print(f"Error: --window must be a positive odd number, got {window}.", file=sys.stderr)
        return 1
    records = _load_records(path)
    if records is None:
        return 1
    written = write_all_figures(records, Path(output_dir), window=window, fmt=fmt)
    print(f"Wrote {len(written)} figures to {output_dir}")
    return 0

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="protein_explorer")
    subcommands = parser.add_subparsers(dest="command", required=True)

    summarize_parser = subcommands.add_parser(
        "summarize", help="Print a validation summary for a FASTA file."
    )
    summarize_parser.add_argument("fasta", help="Path to a protein FASTA file.")

    analyze_parser = subcommands.add_parser(
        "analyze", help="Compute features and write CSV tables."
    )
    analyze_parser.add_argument("fasta", help="Path to a protein FASTA file.")
    analyze_parser.add_argument(
        "--output", default=DEFAULT_OUTPUT_DIR, help="Output directory (default: results/)."
    )
    analyze_parser.add_argument(
        "--ph", type=float, default=DEFAULT_PH, help="pH for net charge (default: 7.0)."
    )

    plot_parser = subcommands.add_parser("plot", help="Generate figures.")
    plot_parser.add_argument("fasta", help="Path to a protein FASTA file.")
    plot_parser.add_argument(
        "--output", default=DEFAULT_FIGURE_DIR, help="Output directory (default: results/figures/)."
    )
    plot_parser.add_argument(
        "--window", type=int, default=DEFAULT_WINDOW,
        help="Hydropathy window, odd (default: 9; use 19 for transmembrane helices).",
    )
    plot_parser.add_argument(
        "--format", dest="fmt", choices=["png", "svg", "pdf"], default="png",
        help="Image format (default: png).",
    )

    args = parser.parse_args(argv)
    if args.command == "summarize":
        return summarize(args.fasta)
    if args.command == "analyze":
        return analyze(args.fasta, args.output, args.ph)
    if args.command == "plot":
        return plot(args.fasta, args.output, args.window, args.fmt)
    return 1

if __name__ == "__main__":
    sys.exit(main())
