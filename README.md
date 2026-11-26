# Bristol-meat-mystery



## Software

- Git Bash (Windows): bash, unzip, sed

- R and RStudio: ape, ggplot2, cluster, DECIPHER (installed with BiocManager)

- Python 3 (VS Code): biopython



## Step 1 – Prepare the data (Bash)

Run in Git Bash:



```

unzip data/sequences\_by\_outlet.zip -d data/

bash scripts/1-PrepareData.sh

```



The script copies the FASTA files from the 8 outlet folders into one folder

(`working_sequences`) and adds the outlet name to each header,

e.g. `>seq013` becomes `>outlet_01|seq013`. The DNA sequences are not changed.


Check (should print 100):


```

cat working\_sequences/\*.fasta | grep -c ">"

```


## Step 2 – Sequence groups (R)

Open `scripts 2-SequenceGroups.Rmd` in RStudio and run all chunks (or Knit).



- The sequences have different lengths, so they are aligned first with `AlignSeqs()` (DECIPHER).

- Pairwise distances: `dist.dna()` with the JC69 model.

- PCoA: `cmdscale()` with 10 axes.

- k-means on the 10 PCoA axes (`nstart = 25`) for k = 2 to 10.

- The number of groups is the k with the highest mean silhouette width.

- `set.seed(123)` is used so k-means gives the same result each time.



Output: `pcoa\_by\_outlet.png`, `pcoa\_by\_group.png`



## Step 3 – Representative sequences (R)

This is the last section of `2-SequenceGroups.Rmd`.



For each group, the representative is the sequence with the smallest mean JC69

distance to the other sequences in its group.



Output: `representatives.fasta`



## Step 4 – ORF prediction and translation (Python)

Run from the `scripts` folder:



```

pip install biopython

python 4-ORFprediction.py

```



Based on the longest ORF code from the practical, with these changes:

- all 6 reading frames are checked (3 forward and 3 reverse)

- the bacterial genetic code is used (`table=11`)

- the ORF starts at the first start codon (ATG, GTG or TTG)

- the longest ORF is selected and translated



Output: `representative\_proteins.fasta`

