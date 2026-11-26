# Bristol-meat-mystery

Software: Git Bash (Windows) – uses bash, unzip and sed.

## Step 1 – Prepare the data

Run in Git Bash:

unzip data/sequences\_by\_outlet.zip -d data/

bash scripts/1-PrepareData.sh



The script copies the FASTA files from the 8 outlet folders into one folder

(working\_sequences) and adds the outlet name to each header,

e.g. >seq013 becomes >outlet\_01|seq013. The DNA sequences are not changed.



Check: cat working\_sequences/\*.fasta | grep -c ">"   (should print 100)

