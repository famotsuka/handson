from Bio.Data import CodonTable
from Bio.Seq import Seq
from Bio import SeqIO
import ast
import csv

# Load dual mutation positions
with open("mutant_sequences/dual_muts.txt") as f:
    dual_muts = ast.literal_eval(f.read().strip())

# Read FASTA
fasta_file = "mutant_sequences/ddlA_frag.fa"
record = next(SeqIO.parse(fasta_file, "fasta"))
seq = str(record.seq).lower()

# Find start of coding region
start_index = seq.find("tatg")
if start_index == -1:
    raise ValueError("No 'tatg' found in sequence.")
# Get flanking sequences
prefix = seq[:start_index + 1]  # before 'tatg'
joint_index = seq.find("gagacc")

# Extract coding region between tatg and (maybe) gagacc
coding_seq = seq[start_index + 1:]
if joint_index != -1 and joint_index >= len(seq) - 30:
    stop_relative_index = joint_index - (start_index + 1)
    if stop_relative_index > 0:
        coding_seq = coding_seq[:stop_relative_index]
        suffix = seq[joint_index:]  # after 'gagacc'

# Trim to full codons
trimmed_seq = coding_seq[:(len(coding_seq) // 3) * 3]

# Get standard codon table
standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

# Invert codon table: amino acid -> list of codons
syn_dict = {}
for codon, aa in standard_table.forward_table.items():
    syn_dict.setdefault(aa, []).append(codon.upper())

# Include stop codons if needed
# syn_dict['*'] = list(standard_table.stop_codons)

# Generate mutated sequences
mutated_records = []
csv_rows = []
print("\n Protein sequences of your double mutants:")

for pos1, pos2 in dual_muts:
    output = f"[{pos1}, {pos2}] -> "

    for pos in (pos1, pos2):
        codon_index = pos - 1
        codon_start = codon_index * 3
        if codon_start + 3 > len(trimmed_seq):
            print(f"Position {pos} is out of range.")
            continue

    # Get wild-type codons and amino acids
    codon1 = trimmed_seq[(pos1-1)*3:(pos1-1)*3+3].upper()
    codon2 = trimmed_seq[(pos2-1)*3:(pos2-1)*3+3].upper()
    aa1 = str(Seq(codon1).translate())
    aa2 = str(Seq(codon2).translate())

    # Get synonymous codons
    syn_codons1 = [c for c in syn_dict.get(aa1, []) if c != codon1]
    syn_codons2 = [c for c in syn_dict.get(aa2, []) if c != codon2]

    # Skip if no synonymous codons available (nothing to mutate)
    if not syn_codons1 or not syn_codons2:
        print(f"Skipping pair [{pos1}, {pos2}] — no synonymous codons to mutate.")
        continue
    for sc1 in syn_codons1:
        for sc2 in syn_codons2:
            # Create a mutable copy of the coding region
            new_seq = list(trimmed_seq.upper())

            # Replace the codons at the specified positions
            new_seq[(pos1 - 1) * 3: (pos1 - 1) * 3 + 3] = sc1
            new_seq[(pos2 - 1) * 3: (pos2 - 1) * 3 + 3] = sc2
            new_seq_str = ''.join(new_seq)

            # Sanity check: translate the coding region only
            mutated_trimmed = new_seq_str[: (len(new_seq_str) // 3) * 3]
            protein_seq = str(Seq(mutated_trimmed).translate(to_stop=False))

            # Build full nucleotide sequence (with prefix/suffix)
            full_seq = prefix + new_seq_str.lower()
            suffix = seq[len(full_seq):]
            full_seq += suffix

            # Highlight mutated codons in uppercase
            mut_pos1_start = len(prefix) + (pos1 - 1) * 3
            mut_pos2_start = len(prefix) + (pos2 - 1) * 3

            # Replace the codons at the correct positions in full_seq (uppercase)
            full_seq = (
                full_seq[:mut_pos1_start] +
                sc1 +
                full_seq[mut_pos1_start + 3:mut_pos2_start] +
                sc2 +
                full_seq[mut_pos2_start + 3:]
            )

            # Build FASTA record (no translation in header)
            fasta_id = f"ddlA_{pos1}{aa1}_{sc1}_{pos2}{aa2}_{sc2}"
            fasta_record = f">{fasta_id}\n{full_seq}"
            mutated_records.append(fasta_record)

            # CSV record
            csv_rows.append([fasta_id, full_seq])

            # Optional: print sanity check (remove or comment if too verbose)
            print(f"{fasta_id} | Protein: {protein_seq}")

# Output to screen (or write to a file if you want)
print("\n Nucleotide sequences of our mutants:")
print("\n".join(mutated_records))

# Save FASTA file
with open("outputs/dual_mutants.fa", "w") as f:
    f.write("\n".join(mutated_records) + "\n")

# Save CSV file
with open("outputs/dual_mutants.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Name", "Sequence"])
    writer.writerows(csv_rows)
