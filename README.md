# VRMOD

## Overview

Vertebrate Regulatory Module Detector (VRMOD) is a computational framework for predicting cis-regulatory modules (CRMs) in vertebrate genomes.

The software identifies clusters of transcription factor binding sites and predicts putative regulatory modules from genomic DNA sequences.

---

## Repository Contents

| File | Description |
|--------|-------------|
| VRMOD_V0.1.pl | Main VRMOD prediction script |
| VRMOD_V0.1_RunOnHTCF.pl | Script for running VRMOD on HTCF |
| alphabet | Alphabet file required by Patser |
| patser-v3e | Patser executable |
| Mouse_phylonet_L30_Short_combined_Strategy2_final.tgz | Position weight matrix collection |
| example | Test files with description |

---

## Requirements

### Software

- Perl (version 5.20 or newer recommended)
- Linux/Unix operating system
- Patser v3e

### Included Files

The repository contains:

- Patser executable (`patser-v3e`)
- Alphabet file (`alphabet`)

---

## Installation

### Step 1: Clone Repository

```bash
git clone https://github.com/GuoyanZhao-Lab/VRMOD.git
cd VRMOD
```

### Step 2: Download PWM Database

Download:

```
Mouse_phylonet_L30_Short_combined_Strategy2_final.tgz
```

Extract:

```bash
tar -xzvf Mouse_phylonet_L30_Short_combined_Strategy2_final.tgz
```

### Step 3: Configure Paths

Open:

```bash
VRMOD_V0.1.pl
```

Modify the following variables:

```perl
my $path_patser = "/path/to/patser-v3e";
my $path_alphabet = "/path/to/alphabet";

my $dir = "/path/to/Mouse_phylonet_L30_Short_combined_Strategy2_final";
```

Replace with the appropriate paths on your system.

---

## Input

VRMOD accepts genomic DNA sequences in FASTA format.

Example:

```fasta
>sequence1
ATGCGTACGATCGATCGATCGATCGATCG
```

---

## Running VRMOD

Example command:

```bash
perl VRMOD_V0.1.pl PWM_directory input.fa
```

Example:

```bash
perl VRMOD_V0.1.pl Mouse_phylonet_L30_Short_combined_Strategy2_final input.fa
```

---

## Output

For each input sequence file, VRMOD generates:

```
input.fa.pred_mod.txt
```

The output contains predicted cis-regulatory modules identified within the input sequence.

Successful completion is indicated by:

```
# program finished
```

at the end of the output file.

---

## Example Analysis (more details in the example folder)

Example input:

```bash
Mus_musculus.GRCm38.dna.chromosomeTEST.19.fa
```

Run:

```bash
perl VRMOD_V0.1.pl example Mus_musculus.GRCm38.dna.chromosomeTEST.19.fa
```

Expected output:

```bash
Mus_musculus.GRCm38.dna.chromosomeTEST.19.fa.pred_mod.txt
```

---

## Troubleshooting

### Patser not found

Check:

```perl
$path_patser
```

points to the correct executable.

### PWM directory not found

Verify:

```perl
$dir
```

matches the location of the extracted PWM files.

### No output generated

Confirm:

- Input FASTA file is correctly formatted.
- PWM directory contains matrix files.
- Patser executable has execute permissions.

---

## Citation

If you use VRMOD in your research, please cite:

Gonçalves T, Stewart C, Baxley S, Xu J, Boyer K, George B, Li D, Yang C, Gabel H, Piao X, Cruchaga C, Li Y, Wang T, Avraham O, Zhao G. Unlocking cis-regulatory landscapes across 500 million years of evolution and disease mechanisms. Manuscript under review at NAR Genomics and Bioinformatics. 2026 (In Revision).

---

## Contact

Guoyan Zhao Lab

GitHub:
https://github.com/GuoyanZhao-Lab/VRMOD

   
