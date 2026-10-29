# Plain-language explanations of each feature, shown in the app's expandable panels.

from dataclasses import dataclass

@dataclass(frozen=True)
class FeatureNote:
    name: str
    meaning: str
    calculation: str
    caveat: str

FEATURE_NOTES: tuple[FeatureNote, ...] = (
    FeatureNote(
        "Molecular weight",
        "The mass of the polypeptide chain, roughly 110 Da per residue.",
        "Sum of average residue masses plus one water for the chain ends.",
        "Calculated from the sequence as written. Signal peptides, propeptides, "
        "glycans, phosphates, and disulfide bonds all change the real mass.",
    ),
    FeatureNote(
        "Isoelectric point (pI)",
        "The pH at which the protein carries no net charge. Proteins are least soluble "
        "near their pI.",
        "Henderson–Hasselbalch charge of every ionisable side chain and both termini, "
        "with the zero crossing found by bisection on pH 0–14. EMBOSS pKa values.",
        "pKa values shift with the local environment in the folded protein. Different "
        "pKa sets give pI values up to about 0.5 units apart.",
    ),
    FeatureNote(
        "Net charge",
        "The protein's overall charge at the chosen pH. It affects solubility, binding, "
        "and how it runs in electrophoresis.",
        "Same Henderson–Hasselbalch sum as the pI, evaluated at one pH.",
        "Counts disulfide-bonded cysteines as ionisable, and ignores modifications such "
        "as phosphorylation.",
    ),
    FeatureNote(
        "GRAVY",
        "Average hydropathy. Positive values are hydrophobic, negative are hydrophilic. "
        "Most soluble proteins fall between about −0.5 and 0.",
        "Mean Kyte–Doolittle value over standard residues.",
        "One number for the whole chain, so it hides hydrophobic segments. It predicts "
        "solubility poorly; insulin is positive but soluble.",
    ),
    FeatureNote(
        "Aromaticity",
        "The share of Phe, Trp, and Tyr. Trp and Tyr absorb UV at 280 nm, which is how "
        "protein concentration is often measured.",
        "(F + W + Y) / length, the same definition as Biopython.",
        "His is not counted, though its ring is aromatic.",
    ),
    FeatureNote(
        "Hydropathy profile",
        "Hydrophobicity along the chain. Peaks are buried or membrane-facing segments; "
        "troughs are surface loops or disordered regions.",
        "Mean Kyte–Doolittle value in a sliding window, plotted at the window's centre. "
        "Window 9 shows surface patches; window 19 matches a membrane-spanning helix.",
        "Signal peptides also cross the 1.6 transmembrane threshold, so a peak alone "
        "does not prove a membrane protein.",
    ),
)
