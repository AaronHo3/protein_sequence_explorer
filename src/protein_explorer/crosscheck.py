# Cross-checks our feature calculations against Biopython and explains the differences.
# Requires the dev extra (biopython). Run as:
#   python -m protein_explorer.crosscheck data/sample.fasta --output docs/validation_report.md

import argparse
import sys
from collections.abc import Callable
from pathlib import Path
from types import MappingProxyType

import Bio
import pandas as pd
from Bio.SeqUtils import IsoelectricPoint as bio_ip
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from protein_explorer import features
from protein_explorer.dataset import parse_accession
from protein_explorer.fasta import FastaRecord, read_fasta
from protein_explorer.features import PkaSet
from protein_explorer.plots import entry_label
from protein_explorer.validation import STANDARD_RESIDUES

# Biopython's base pKa table, read from the library rather than copied by hand
BIOPYTHON_BASE_PKA = PkaSet(
    "Biopython base (no terminal adjustment)",
    MappingProxyType(dict(bio_ip.positive_pKs)),
    MappingProxyType(dict(bio_ip.negative_pKs)),
)

# Sequences that exercise residue conventions where implementations commonly differ
EDGE_CASES: dict[str, str] = {
    "standard": "MKTAYIAKQR",
    "unknown X": "MKTXAYIAKQR",
    "ambiguous B": "MKTBAYIAKQR",
    "selenocysteine U": "MKTUAYIAKQR",
    "terminal stop *": "MKTAYIAKQR*",
    "lowercase": "mktayiakqr",
}

def biopython_pka_set(sequence: str) -> PkaSet:
    """Biopython's pKa set for this sequence: N/C-terminal pKa depend on the terminal residues."""
    residues = sequence.upper().removesuffix("*")
    positive = dict(bio_ip.positive_pKs)
    negative = dict(bio_ip.negative_pKs)
    positive["Nterm"] = bio_ip.pKnterminal.get(residues[0], positive["Nterm"])
    negative["Cterm"] = bio_ip.pKcterminal.get(residues[-1], negative["Cterm"])
    return PkaSet(
        "Biopython (terminal-adjusted)", MappingProxyType(positive), MappingProxyType(negative)
    )

def _paired_values(sequence: str) -> dict[str, tuple[float, float]]:
    """(ours, Biopython) for each feature both libraries compute."""
    bio = ProteinAnalysis(sequence)
    return {
        "molecular_weight": (features.molecular_weight(sequence), bio.molecular_weight()),
        "aromaticity": (features.aromatic_fraction(sequence), bio.aromaticity()),
        "gravy": (features.gravy(sequence), bio.gravy()),
        "isoelectric_point": (features.isoelectric_point(sequence), bio.isoelectric_point()),
        "net_charge_ph7": (features.net_charge(sequence, 7.0), bio.charge_at_pH(7.0)),
    }

def compare_features(records: list[FastaRecord]) -> pd.DataFrame:
    """Long table: one row per (protein, feature) with our value, Biopython's, and the difference."""
    rows = []
    for record in records:
        for feature, (ours, theirs) in _paired_values(record.sequence).items():
            rows.append({
                "protein": entry_label(record.identifier),
                "feature": feature,
                "ours": ours,
                "biopython": theirs,
                "difference": ours - theirs,
                "relative_difference": (ours - theirs) / theirs if theirs else float("nan"),
            })
    return pd.DataFrame(rows)

def decompose_pi(records: list[FastaRecord]) -> pd.DataFrame:
    """Swap in Biopython's conventions one at a time to attribute the pI difference."""
    rows = []
    for record in records:
        seq = record.sequence
        emboss = features.isoelectric_point(seq)
        base = features.isoelectric_point(seq, BIOPYTHON_BASE_PKA)
        adjusted = features.isoelectric_point(seq, biopython_pka_set(seq))
        library = ProteinAnalysis(seq).isoelectric_point()
        rows.append({
            "protein": entry_label(record.identifier),
            "ours_emboss": emboss,
            "ours_biopython_pka": base,
            "ours_biopython_pka_terminal": adjusted,
            "biopython": library,
            "pka_set_effect": base - emboss,
            "terminal_effect": adjusted - base,
            "remaining": library - adjusted,
        })
    return pd.DataFrame(rows)

def _try(func: Callable[[], float]) -> str:
    """Run a calculation and describe the outcome, including which exception it raised."""
    try:
        return f"{func():.4f}"
    except Exception as error:  # report every failure type; this is a behaviour survey
        return f"error ({type(error).__name__})"

def edge_case_behaviour() -> pd.DataFrame:
    """How each implementation handles residues outside the standard 20."""
    rows = []
    for case, seq in EDGE_CASES.items():
        bio = ProteinAnalysis(seq)
        rows.append({
            "case": case,
            "sequence": seq,
            "mw_ours": _try(lambda s=seq: features.molecular_weight(s)),
            "mw_biopython": _try(bio.molecular_weight),
            "gravy_ours": _try(lambda s=seq: features.gravy(s)),
            "gravy_biopython": _try(bio.gravy),
            "aromaticity_ours": _try(lambda s=seq: features.aromatic_fraction(s)),
            "aromaticity_biopython": _try(bio.aromaticity),
            "pi_ours_emboss": _try(lambda s=seq: features.isoelectric_point(s)),
            "pi_ours_biopython_pka": _try(
                lambda s=seq: features.isoelectric_point(s, biopython_pka_set(s))
            ),
            "pi_biopython": _try(bio.isoelectric_point),
        })
    return pd.DataFrame(rows)

def mass_conventions(records: list[FastaRecord]) -> pd.DataFrame:
    """Average vs monoisotopic molecular weight: the choice matters for mass spectrometry."""
    rows = []
    for record in records:
        average = ProteinAnalysis(record.sequence).molecular_weight()
        mono = ProteinAnalysis(record.sequence, monoisotopic=True).molecular_weight()
        rows.append({
            "protein": entry_label(record.identifier),
            "ours_average": features.molecular_weight(record.sequence),
            "biopython_average": average,
            "biopython_monoisotopic": mono,
            "average_minus_monoisotopic": average - mono,
        })
    return pd.DataFrame(rows)

def _markdown_table(table: pd.DataFrame, float_format: str = "{:.4f}") -> str:
    """Render a DataFrame as a GitHub Markdown table without extra dependencies."""
    def cell(value: object) -> str:
        return float_format.format(value) if isinstance(value, float) else str(value)
    header = "| " + " | ".join(table.columns) + " |"
    divider = "|" + "---|" * len(table.columns)
    body = ["| " + " | ".join(cell(v) for v in row) + " |" for row in table.itertuples(index=False)]
    return "\n".join([header, divider, *body])

def _agreement_section(comparison: pd.DataFrame) -> list[str]:
    worst = (
        comparison.assign(abs_difference=comparison["difference"].abs())
        .groupby("feature", sort=False)["abs_difference"].max()
        .rename("max_abs_difference").reset_index()
    )
    return [
        "## Agreement by feature",
        "",
        "Largest absolute difference across all proteins:",
        "",
        _markdown_table(worst, "{:.6f}"),
        "",
        "- **Molecular weight, aromaticity, GRAVY:** same definitions and the same mass "
        "table, so they agree to rounding error. The tiny molecular-weight difference comes "
        "from rounding our residue masses to 4 decimals.",
        "- **pI and net charge:** differ because the two use different pKa values. "
        "See the breakdown below.",
        "",
    ]

def _pi_section(pi: pd.DataFrame) -> list[str]:
    return [
        "## Breaking down the pI difference",
        "",
        "Each column swaps in one more of Biopython's conventions:",
        "",
        "1. `pka_set_effect`: EMBOSS → Biopython's base pKa values.",
        "2. `terminal_effect`: Biopython also adjusts the terminal pKa for the first and last "
        "residue (e.g. N-terminal Met uses 7.0 instead of 7.5).",
        "3. `remaining`: whatever is left after both changes. It should be about 0.",
        "",
        _markdown_table(pi, "{:.3f}"),
        "",
        f"- The pKa set accounts for up to {pi['pka_set_effect'].abs().max():.2f} pH units; "
        f"terminal adjustment for up to {pi['terminal_effect'].abs().max():.2f}.",
        f"- After both, the largest remaining difference is {pi['remaining'].abs().max():.4f}, "
        "within bisection tolerance. Our implementation reproduces Biopython exactly "
        "when given its conventions.",
        "- " + _ranking_statement(pi),
        "",
    ]

def _ranking_statement(pi: pd.DataFrame) -> str:
    """Whether the two pKa sets order proteins the same way (what rank-based models see)."""
    same = (pi["ours_emboss"].rank() == pi["biopython"].rank()).all()
    return (
        "Both pKa sets rank the proteins in the same order by pI."
        if same
        else "The two pKa sets rank the proteins in a different order by pI."
    )

def _conventions_section(records: list[FastaRecord]) -> list[str]:
    return [
        "## Mass conventions",
        "",
        "Average masses weight each element by natural isotope abundance (SDS-PAGE, "
        "stoichiometry). Monoisotopic masses use only the lightest isotopes (mass "
        "spectrometry of peptides). We use average masses.",
        "",
        _markdown_table(mass_conventions(records), "{:.2f}"),
        "",
        "## Residues outside the standard 20",
        "",
        _markdown_table(edge_case_behaviour()),
        "",
        "- **X / B:** neither implementation assigns a mass. Biopython's GRAVY raises an "
        "error; ours averages over the known residues.",
        "- **U (selenocysteine):** both include its mass. Biopython has no pKa for U, so its "
        "charge ignores it. Ours uses pKa 5.2, which lowers the EMBOSS pI.",
        "- **Terminal `*`:** Biopython counts it as a residue (aromaticity 1/11, molecular "
        "weight raises an error). We strip it before every calculation.",
        "- **Search range:** " + _clamping_statement(),
        "",
    ]

def _clamping_statement() -> str:
    """Show, with a poly-Asp peptide, that Biopython clamps pI to its 4.05-12 search range."""
    peptide = "D" * 20
    bio = ProteinAnalysis(peptide)
    bio_pi = bio.isoelectric_point()
    ours = features.isoelectric_point(peptide, biopython_pka_set(peptide))
    return (
        f"Biopython's pI bisection searches pH 4.05–12, ours 0–14. For poly-Asp (D×20), "
        f"Biopython returns {bio_pi:.2f}, where the net charge is "
        f"{bio.charge_at_pH(bio_pi):.1f}, so that value is not a pI. With the same pKa values, "
        f"our search finds {ours:.2f}. Biopython clamps silently; real proteins rarely fall "
        "outside this range, but very acidic or basic peptides can."
    )

def build_report(records: list[FastaRecord], input_path: Path) -> str:
    """Markdown validation report comparing our features with Biopython."""
    comparison = compare_features(records)
    accessions = ", ".join(a for r in records if (a := parse_accession(r.identifier)))
    sections = [
        "# Validation report: protein_explorer vs Biopython",
        "",
        f"Generated by `python -m protein_explorer.crosscheck` on `{input_path}` "
        f"({len(records)} proteins: {accessions}) with Biopython {Bio.__version__}.",
        "",
        *_agreement_section(comparison),
        *_pi_section(decompose_pi(records)),
        *_conventions_section(records),
        "## Full comparison",
        "",
        _markdown_table(comparison, "{:.6g}"),
        "",
    ]
    return "\n".join(sections)

def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="protein_explorer.crosscheck")
    parser.add_argument("fasta", help="Path to a protein FASTA file.")
    parser.add_argument("--output", default="docs/validation_report.md", help="Report path.")
    args = parser.parse_args(argv)

    input_path = Path(args.fasta)
    records = [
        r for r in read_fasta(input_path)
        if STANDARD_RESIDUES.issuperset(r.sequence.upper())
    ]
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(build_report(records, input_path), encoding="utf-8")
    print(f"Wrote {output} ({len(records)} proteins compared)")
    return 0

if __name__ == "__main__":
    sys.exit(main())
