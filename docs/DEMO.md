# Demo script (about 2½ minutes)

A recorded walkthrough of Protein Sequence Explorer. Each step lists what to do on screen and what to say. Timings are approximate.

## Before recording

```bash
uv pip install -e ".[dev,app]"
uv run streamlit run app/streamlit_app.py
```

- Use a browser window about 1440 × 900, light mode, at 100% zoom.
- Start on the Import page with the example dataset loaded. If a previous session left filters or selections, click **Reset workspace**.
- Have `tests/fixtures/mixed_validity.fasta` ready to drag in.

## 1. Load a real dataset (0:00–0:20)

**Do:** On **Import**, point at the eight UniProt proteins and the all-valid status. Drag `mixed_validity.fasta` onto the uploader, show the warning and the two failed records with their reasons, then click **Load example dataset**.

**Say:** "This loads protein FASTA files and validates every record. A sequence with an alignment gap, or a DNA sequence labelled as protein, is reported with the reason, never silently fixed or dropped. I'll go back to eight reviewed UniProt proteins."

## 2. Dataset summary (0:20–0:40)

**Do:** Open **Overview**. Hover a tick in the rug under the pI histogram, then the linear fit on the molecular weight chart.

**Say:** "Eight proteins, about nineteen hundred residues. Each tick under the histogram is one protein. Molecular weight tracks length at about 114 daltons per residue, close to the textbook 110."

## 3. Select a protein (0:40–0:50)

**Do:** In the sidebar protein picker, type "lys" and choose **LYSC_CHICK, Lysozyme C**. Open **Protein**.

**Say:** "Lysozyme: an estimated pI of 9.25 and a charge of about plus nine at pH 7. It's a cationic enzyme that binds negatively charged bacterial cell walls."

## 4. Inspect residues (0:50–1:15)

**Do:** Hover a few residues. Click residue 1, then shift-click residue 18.

**Say:** "Every residue shows its name, codes, class, and hydropathy. Selecting residues 1 to 18 shows a stretch that's 72% hydrophobic, with a mean hydropathy above +2. That's lysozyme's signal peptide, which is cleaved off before the protein is secreted."

## 5. Switch colouring modes (1:15–1:30)

**Do:** Change **Colour residues by** to **Positive vs negative charge**, then to **Hydrophobic vs hydrophilic**.

**Say:** "Charge colouring shows at a glance why lysozyme is basic: lysines and arginines everywhere, few acidic residues. Each mode uses the same residue definitions as the calculated features."

## 6. Hydropathy profile (1:30–1:50)

**Do:** Scroll to **Hydropathy profile**. Point at the shaded selection over residues 1–18. Move the window slider from 9 to 19. Drag across another peak, around residues 95–110, and scroll up to show it selected in the viewer.

**Say:** "Bars are single residues; the line is a sliding-window mean. The selection is linked both ways: dragging on the plot selects those residues in the sequence. The transmembrane reference is off by default, because a signal peptide crosses it too. A peak alone is not evidence of a membrane helix."

## 7. Compare proteins (1:50–2:15)

**Do:** Open **Compare**. Set the axes to isoelectric point and GRAVY. Click the CALM1_HUMAN point. Scroll to **Normalized properties**, then to **Composition heatmap** and choose **Clustered by composition**.

**Say:** "Calmodulin sits at the acidic end, with a pI below 4. Its calcium-binding loops are full of aspartate and glutamate. Clicking it selects it everywhere. The normalized heatmap compares properties with different units on one z-score scale, and the composition heatmap clusters proteins by amino acid makeup."

## 8. Export (2:15–2:30)

**Do:** Open **Export**. Download **Features** and **Analysis report (Markdown)**, then click **Prepare figure bundle**.

**Say:** "Everything is exportable: feature tables ready to use as a machine-learning feature matrix, FASTA that round-trips the input exactly, and a report recording the input checksum, pH, and pKa set, so every number can be reproduced."

## Optional extras

- **Dark mode:** switch the system theme. Chart colours change to a separately validated dark palette.
- **Long sequences:** upload titin (UniProt Q8WZ42, 34,350 residues) to show the viewer handling a very long protein.
- **Validation:** open `docs/validation_report.md` to show agreement with Biopython and the breakdown of the pI difference.
