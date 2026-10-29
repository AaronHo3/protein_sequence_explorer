# Data provenance

`sample.fasta` — 8 reviewed (Swiss-Prot) protein sequences downloaded from UniProtKB on 2026-10-28.

| Accession | Entry name   | Protein                    | Organism          |
|-----------|--------------|----------------------------|-------------------|
| P00698    | LYSC_CHICK   | Lysozyme C                 | *Gallus gallus*   |
| P01308    | INS_HUMAN    | Insulin (preproinsulin)    | *Homo sapiens*    |
| P02769    | ALBU_BOVIN   | Serum albumin              | *Bos taurus*      |
| P04637    | P53_HUMAN    | Cellular tumor antigen p53 | *Homo sapiens*    |
| P0DP23    | CALM1_HUMAN  | Calmodulin-1               | *Homo sapiens*    |
| P42212    | GFP_AEQVI    | Green fluorescent protein  | *Aequorea victoria* |
| P61769    | B2MG_HUMAN   | Beta-2-microglobulin       | *Homo sapiens*    |
| P68871    | HBB_HUMAN    | Hemoglobin subunit beta    | *Homo sapiens*    |

Source: The UniProt Consortium, *UniProt: the Universal Protein Knowledgebase*, https://www.uniprot.org (CC BY 4.0).

Reproduce with:

```bash
curl -s "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=accession:P01308+OR+accession:P68871+OR+accession:P00698+OR+accession:P42212+OR+accession:P04637+OR+accession:P0DP23+OR+accession:P02769+OR+accession:P61769" -o data/sample.fasta
```

Note: these are full UniProt precursor sequences (they include signal peptides and propeptides), not the mature proteins.
