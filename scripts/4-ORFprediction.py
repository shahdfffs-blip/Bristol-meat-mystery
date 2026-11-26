# Step 4 - ORF prediction and translation
# Based on the longest ORF code from the Python for Bioinformaticians practical (8.3).
# Changes I made:
#   - read all 6 reading frames (3 forward + 3 on the reverse strand), not only 3
#   - use the bacterial genetic code (table 11) because the samples are bacteria
#   - start the ORF at the first start codon (ATG, GTG or TTG in table 11)

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

start_codons = ["ATG", "GTG", "TTG"]

records = list(SeqIO.parse("representatives.fasta", "fasta"))
proteins = []

for record in records:
    print("\n===", record.id, "===")

    best_orf = ""      # longest stretch of amino acids without a stop codon
    best_dna = ""      # the DNA of that stretch
    best_strand = None
    best_frame = None

    # forward strand and reverse strand
    for strand, dna_seq in [("+", record.seq), ("-", record.seq.reverse_complement())]:

        # three reading frames on each strand
        for frame in range(3):
            frame_seq = dna_seq[frame:]
            # ignore incomplete codon at the end (as in the bonus solution)
            frame_seq = frame_seq[:int(len(frame_seq) / 3) * 3]
            protein = frame_seq.translate(table=11)
            fragments = str(protein).split("*")

            position = 0   # where each fragment starts (in amino acids)
            for fragment in fragments:
                if len(fragment) > len(best_orf):
                    best_orf = fragment
                    best_strand = strand
                    best_frame = frame
                    start = frame + position * 3
                    # + 3 to keep the stop codon
                    end = start + len(fragment) * 3 + 3
                    best_dna = dna_seq[start:end]
                position = position + len(fragment) + 1

            longest = max(fragments, key=len)
            print(
                f"strand {strand} frame {frame}: longest ORF = {len(longest)} aa")

    # move to the first start codon in the longest ORF
    orf_dna = best_dna
    for i in range(0, len(best_dna), 3):
        if str(best_dna[i:i + 3]) in start_codons:
            orf_dna = best_dna[i:]
            break

    start_codon = str(orf_dna[0:3])
    has_stop = str(orf_dna.translate(table=11)).endswith("*")

    protein = orf_dna.translate(table=11, to_stop=True)
    protein = "M" + protein[1:]   # the first codon is read as methionine

    print(f"Selected ORF: strand {best_strand}, frame {best_frame}, "
          f"start codon {start_codon}, stop codon {has_stop}, length {len(protein)} aa")

    proteins.append(SeqRecord(protein, id=record.id,
                    description="predicted protein"))

SeqIO.write(proteins, "representative_proteins.fasta", "fasta")
