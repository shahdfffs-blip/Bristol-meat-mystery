# Step 6 - Combine sequences and clean the FASTA headers
# Puts the 5 representative proteins, the 20 references and the outgroup in one file.
# Headers have no spaces or special characters (as in the Week 3 practical),
# because they become the tip labels of the tree.
# Written with plain Python (files, loops, strings, dictionaries), without Biopython.


def read_fasta(filename):
    # returns a dictionary: full header line -> sequence
    sequences = {}
    header = ""
    for line in open(filename):
        line = line.strip()
        if line.startswith(">"):
            header = line[1:]
            sequences[header] = ""
        elif line != "":
            sequences[header] = sequences[header] + line
    return sequences


# sequence group of each representative (from Step 2)
groups = {"seq092": "Group1", "seq077": "Group2", "seq033": "Group3",
          "seq074": "Group4", "seq053": "Group5"}

# Swiss-Prot headers have no species name, so I add it here
swissprot_species = {"Q88XZ3.2": "Lactiplantibacillus plantarum WCFS1",
                     "Q0TA78.1": "Escherichia coli 536"}

outgroup = "Q0TA78.1"

all_sequences = {}   # new header -> sequence

# representatives, e.g. "outlet_08|seq092 predicted protein" -> "Group1_outlet_08_seq092"
representatives = read_fasta("representative_proteins.fasta")
for header in representatives:
    first_word = header.split()[0]
    outlet, seq_id = first_word.split("|")
    new_id = groups[seq_id] + "_" + outlet + "_" + seq_id
    all_sequences[new_id] = representatives[header]

# references, e.g. "WP_260210480.1 ... [Lactiplantibacillus pentosus]"
#              -> "Lactiplantibacillus_pentosus_WP_260210480.1"
references = read_fasta("references.fasta")
for header in references:
    first_word = header.split()[0]
    if first_word.startswith("sp|"):
        accession = first_word.split("|")[1]
        species = swissprot_species[accession]
    else:
        accession = first_word
        species = header.split("[")[1].replace("]", "")

    new_id = species.replace(" ", "_") + "_" + accession
    if accession == outgroup:
        new_id = "OUTGROUP_" + new_id
    all_sequences[new_id] = references[header]

# write the new file
output = open("all_sequences.fasta", "w")
for new_id in all_sequences:
    print(new_id)
    output.write(">" + new_id + "\n" + all_sequences[new_id] + "\n")
output.close()

print("Total sequences:", len(all_sequences))
