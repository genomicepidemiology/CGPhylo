from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner
import json

def align_sequences(seq_a, seq_b):
    aligner = PairwiseAligner()
    aligner.mode = 'global'  # Use global alignment
    # Set the gap penalties
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    # Set the end gap penalties
    aligner.target_end_gap_score = -10.0
    aligner.query_end_gap_score = -10.0

    alignments = aligner.align(seq_a, seq_b)
    gap_positions_a, gap_positions_b = extract_gap_positions(alignments[0])
    # Check if the aligned sequences have the same length
    if len(gap_positions_a) != len(gap_positions_b):
        print("Error: Aligned sequences do not have the same length.")

    return gap_positions_a, gap_positions_b, alignments[0].format()

def extract_gap_positions(alignment):
    seqs_alignment = alignment.format()
    seqs_a_alignment = seqs_alignment.split('\n')[0]
    seqs_b_alignment = seqs_alignment.split('\n')[2]
    return seqs_a_alignment, seqs_b_alignment

def find_gap_positions(aligned_sequence):
    gap_positions = []
    for i, char in enumerate(aligned_sequence):
        if char == '-':
            gap_positions.append(str(i))
    return ','.join(gap_positions)

def find_gap_strings(fasta_file):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    gene_sequences = {}
    gene_alignments = {}


    # Group sequences by gene name
    for header, seq_record in sequences.items():
        gene_sequences[header] = str(seq_record.seq)
        gene_alignments[header] = {}

    t = 0
    for header_a in gene_alignments:
        gene_a = header_a.split('_')[:-2]
        gene_a = '_'.join(gene_a)
        for header_b in gene_alignments:
            gene_b = header_b.split('_')[:-2]
            gene_b = '_'.join(gene_b)
            if gene_a == gene_b and header_a != header_b:
                if header_b not in gene_alignments[header_a]:
                    gap_positions_a, gap_positions_b, alignment = align_sequences(gene_sequences[header_a], gene_sequences[header_b])
                    a_gaps = find_gap_positions(gap_positions_a)
                    b_gaps = find_gap_positions(gap_positions_b)
                    gene_alignments[header_a][header_b] = [a_gaps, b_gaps]
                    gene_alignments[header_b][header_a] = [b_gaps, a_gaps]
                    t += 1
                    if t % 100 == 0:
                        print(t)

    return gene_alignments

gene_alignments = find_gap_strings("b2250.fasta")

# Convert to JSON and save to a file
#print (gene_alignments)
with open('gap_map.json', 'w') as outfile:
    json.dump(gene_alignments, outfile, indent=4)

print("Alignment data saved to gap_map.json")