##Example files for testing VRMOD

The repository includes example files that can be used to verify software installation and reproduce a test analysis.

| File Type | Description |
|------------|-------------|
| `Mus_musculus.GRCm38.dna.chromosomeTEST.19.fa` | Example mouse genomic DNA sequence from chr19 in FASTA format used as input for CRM prediction. |
| `*.matrix` | Position Weight Matrix (PWM) files representing conserved transcription factor binding motifs used by VRMOD to identify potential binding sites. |
| `*.desc` | Metadata files associated with each PWM, containing motif statistics, consensus sequence information, and the list of genes that are potentially regulated by this motif. |

### Position Weight Matrix Files (`.matrix`)

The `.matrix` files contain Position Weight Matrices (PWMs) that describe the nucleotide preferences of a conserved DNA motif at each position. These matrices are used by VRMOD and Patser to scan genomic sequences for potential transcription factor binding sites.

### Motif Description Files (`.desc`)

Each `.desc` file contains information describing the corresponding PWM, including:

- Motif identifier
- Motif width (number of nucleotide positions)
- Score threshold
- Information content (IC)
- Statistical significance (P-value)
- Consensus sequence
- Position-specific nucleotide frequencies
- Gene cluster size
- Gene identifiers associated with the motif

Example:

```text
ENSMUSG00000001558_ENSMUST00000001599_1.13

Width     = 22
IC        = 14.846
P-value   = 3.04087E-57
Consensus = tAtRAAAGATGtagtAGGgAAa

Gene Cluster size = 37
```

This example represents a highly conserved 22-bp motif derived from a cluster of 37 homologous genes. The consensus sequence summarizes the most frequent nucleotide observed at each position, while the PWM provides the quantitative nucleotide frequencies used during motif scanning.

### Running the Example Analysis

Example command:

```bash
perl VRMOD_V0.1.pl example \
Mus_musculus.GRCm38.dna.chromosomeTEST.19.fa
```

Expected output:

```text
Mus_musculus.GRCm38.dna.chromosomeTEST.19.fa.pred_mod.txt
```

Successful execution is indicated by the message:

```text
# program finished
```

at the end of the output file.
