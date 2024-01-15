from Bio import SeqIO

def calculate_average_seqlength(input_file):
    sequences = list(SeqIO.parse(input_file, "fasta"))
    total_length = sum(len(seq_record.seq) for seq_record in sequences)
    return total_length / len(sequences)

def filter_by_length(input_file, minimum_length):
    good_seqs = []
    for seq_record in SeqIO.parse(input_file, "fasta"):
        if len(seq_record.seq) > minimum_length:
            good_seqs.append(seq_record)
    return good_seqs

input_file = "combined_slph_for_alignment.fasta"
output_file = "seq_right_length.fasta"
min_length = 1300

filtered_sequences = filter_by_length(input_file, min_length)

average_length = calculate_average_seqlength(input_file)
print(f"Average sequence length: {average_length:.2f} bases")

SeqIO.write(filtered_sequences, output_file, "fasta")
print(f"{len(filtered_sequences)} sequences above {min_length} bases written to {output_file}")

