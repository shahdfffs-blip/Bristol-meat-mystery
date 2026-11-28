# Bristol-meat-mystery

MSc Bioinformatics assignment. 100 DNA sequences (one protein-coding region per bacterial isolate)
from 8 food outlets were analysed to find how many bacterial sequence groups there are, which protein
the sequences encode, which taxa the groups belong to, and which outlets may be linked to the outbreak.

## Repository structure

```
data/
  sequences_by_outlet.zip        original data (8 outlet folders, 100 FASTA files)
working_sequences/               Step 1 output (100 FASTA files)
scripts/
  1-PrepareData.sh               Step 1 (Bash)
  2-SequenceGroups.Rmd           Steps 2 and 3 (R)
  4-ORFprediction.py             Step 4 (Python)
  6-CleanHeaders.py              Step 6, combine sequences and clean headers (Python)
  7-PhylogeneticTree.Rmd         Step 7 (R)
  pcoa_by_outlet.png             Step 2 output (Figure 1)
  pcoa_by_group.png              Step 2 output (Figure 2)
  representatives.fasta          Step 3 output (DNA of the 5 representatives)
  representative_proteins.fasta  Step 4 output (translated ORFs)
  blast/nr/                      Step 5, BLASTP results against nr (CSV)
  blast/outgroup/                Step 5, BLASTP results against Swiss-Prot, limited to E. coli (CSV)
  references.fasta               Step 5, 20 references + outgroup downloaded from NCBI
  all_sequences.fasta            Step 6 output (26 sequences with clean headers)
  alignment.aln                  Step 6, MSA before trimming
  alignment_clean.aln            Step 6, MSA after Gblocks
  tree.png                       Step 7 output, rooted tree (Figure 3)
  tree.nwk                       Step 7 output, rooted tree in Newick format
```

All scripts and their outputs are kept in the `scripts` folder, because each step reads the
output of the step before it from the same folder.

## Software

- Git Bash (Windows): bash, unzip, sed
- R and RStudio: ape, ggplot2, cluster, DECIPHER (installed with BiocManager), phangorn
- Python 3: biopython (Step 4 only)
- Web tools: NCBI BLASTP, NCBI Protein, Phylogeny.fr

Versions used: R 4.6.1, Python 3.14.9, biopython 1.88
(R package versions are printed at the end of each Rmd file.)

## How to run the analysis

### Step 1 – Prepare the data (Bash)
Run in Git Bash from the main folder of the repository:

```
unzip data/sequences_by_outlet.zip -d data/
bash scripts/1-PrepareData.sh
```

Copies the FASTA files from the 8 outlet folders into `working_sequences` and adds the outlet
name to each header, e.g. `>seq013` becomes `>outlet_01|seq013`. The DNA is not changed.

Check (should print 100):

```
cat working_sequences/*.fasta | grep -c ">"
```

### Steps 2 and 3 – Sequence groups and representatives (R)
Open `scripts/2-SequenceGroups.Rmd` in RStudio and run all chunks.

- Sequences aligned with `AlignSeqs()` (DECIPHER), because they have different lengths.
- Pairwise distances: `dist.dna()`, model JC69.
- PCoA: `cmdscale()`, 10 axes.
- k-means on the 10 PCoA axes, `nstart = 25`, k = 2 to 10, `set.seed(123)`.
- Number of groups = k with the highest mean silhouette width (k = 5, 0.811).
- Representative of each group = sequence with the smallest mean JC69 distance to the other
  sequences in its group.

### Step 4 – ORF prediction and translation (Python)
Run from the `scripts` folder:

```
python 4-ORFprediction.py
```

Based on the longest ORF code from the Python practical, changed to:
- check all 6 reading frames (3 forward and 3 reverse)
- use the bacterial genetic code (`table=11`)
- start the ORF at the first start codon (ATG, GTG or TTG)
- select the longest ORF and translate it

### Step 5 – BLASTP (manual, NCBI website)
See "Manual steps" below.

### Step 6 – Multiple sequence alignment
1. Combine the reference files downloaded from NCBI into `scripts/references.fasta`, then run
   from the `scripts` folder:

```
python 6-CleanHeaders.py
```

   This writes `all_sequences.fasta` (5 representatives + 20 references + outgroup) with no spaces
   or special characters in the headers.
2. Align on Phylogeny.fr (see "Manual steps").

### Step 7 – Phylogenetic tree (R)
Open `scripts/7-PhylogeneticTree.Rmd` in RStudio and run all chunks.

- Input: `scripts/alignment_clean.aln`
- Distances: `dist.ml()`, protein model JTT.
- Tree: neighbour joining, `nj()`; negative branch lengths set to 0.
- Rooted with the outgroup (*E. coli*) using `root()`.
- 100 bootstrap replicates (`bootstrap.phyDat`, `set.seed(123)`), shown on the nodes.
- Representative sequences shown in red.

## Manual steps

### Step 5 – BLASTP (NCBI, 26–27 November 2026)
1. BLASTP of the 5 proteins in `representative_proteins.fasta` against **nr**,
   Max target sequences = 5000. CSV files saved in `scripts/blast/nr/`.
2. Four references were chosen for each representative using these cut-offs:
   E-value ≤ 1e-5, query coverage ≥ 97%, percent identity ≥ 77%.
   Partial sequences, MAGs, unnamed species ("sp.") and hits with doubtful taxonomy were excluded.
   For breadth: closest named species, another species of the same or closest genus,
   another genus of the same family, and a more distant genus.
3. Outgroup: BLASTP of the 5 proteins against **UniProtKB/Swiss-Prot**, Organism =
   *Escherichia coli* (taxid:562). Chosen: Q0TA78.1 (*E. coli* 536). CSV files saved in `scripts/blast/outgroup/`.
4. The 21 accessions were downloaded in FASTA format from NCBI Protein (Send to > File > FASTA).
   HHA6352772.1 and MGP5079620.1 (WGS records) were downloaded separately.
   MGO1758743.1 was replaced by WP_224183614.1 (*Macrococcus armenti*) because it was a MAG.

Accessions: Q88XZ3.2, WP_260210480.1, WP_065937617.1, KRN81423.1, HHA6352772.1, WP_031787941.1,
WP_224183614.1, WP_340403060.1, WP_012751038.1, WP_345360685.1, WP_278856784.1, WP_481737786.1,
MGP5079620.1, WP_096798595.1, WP_058235310.1, WP_152230074.1, WP_013582405.1, WP_251868742.1,
WP_415694640.1, WP_485090712.1, Q0TA78.1 (outgroup).

### Step 6 – Phylogeny.fr ("A la Carte")
`scripts/all_sequences.fasta` was uploaded with MUSCLE (alignment), Gblocks (curation) and default
settings. `alignment.aln` and `alignment_clean.aln` were downloaded to `scripts/`.
The alignment had 1579 positions before and 933 after Gblocks. Only the trimmed alignment was used;
the tree was made in R (Step 7).

## Reflection on the methods

- **Bash for Step 1**, as asked in the brief. A simple `for` loop and `sed` change only the header
  lines, so the DNA stays the same.
- **Alignment before distances (Step 2).** `dist.dna()` compares sequences position by position,
  so the sequences of different lengths had to be aligned first. I used DECIPHER in R because
  Phylogeny.fr (used in the Week 3 practical) accepts a maximum of 50 sequences.
- **JC69 distances, PCoA, k-means and silhouette** follow the Phylogenetics practical.
  k was chosen by silhouette and not fixed in advance. Groups 3, 4 and 5 overlapped on PCoA axes
  1 and 2, so I also plotted axes 3 and 4, where they were clearly separated.
- **Representatives** were the most central sequence of each group (smallest mean distance),
  so they best represent the group. Each was also the longest in its group, so not truncated.
- **ORF prediction** used the practical's longest-ORF code. I added the reverse strand and the
  bacterial genetic code because the samples are bacteria, and many sequences start with GTG or TTG.
- **BLASTP against nr** was used to identify the protein and choose references, and the cut-offs
  are reported as recommended in the Week 2 practical summary.
- **The outgroup** is outside both phyla in the dataset (Bacillota and Actinomycetota) but is
  still a homolog of all five proteins, following the Week 3 practical.
- **MUSCLE + Gblocks** are tools named in the Week 3 lecture. Gblocks removed poorly aligned and
  gapped regions, such as extra regions found only in *E. coli*.
- **Neighbour joining in R** was used for the tree because it is the method used in the practicals.
  R was chosen over Phylogeny.fr + iTOL so the tree, rooting and bootstrap are reproducible from
  code. JTT is a protein model, because the tree is built from protein sequences.
- **Limitations:** the tree is based on a single gene (rpoB), NJ is a simple method, and rpoB
  cannot show whether a strain carries toxin genes.
