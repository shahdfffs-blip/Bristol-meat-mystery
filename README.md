# Bristol-meat-mystery

Software: Git Bash (Windows) – uses bash, unzip and sed.

## Step 1 – Prepare the data

Run in Git Bash:

unzip data/sequences_by_outlet.zip -d data/

bash scripts/01_prepare_data.sh



The script copies the FASTA files from the 8 outlet folders into one folder

(working_sequences) and adds the outlet name to each header,

e.g. >seq013 becomes >outlet_01|seq013. The DNA sequences are not changed.



Check: cat working_sequences/*.fasta | grep -c ">"   (should print 100)

