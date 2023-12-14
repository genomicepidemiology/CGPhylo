from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner
import json

def align_sequences(seq_a, seq_b):
    aligner = PairwiseAligner()
    aligner.mode = 'local'  # Use local alignment
    alignments = aligner.align(seq_a, seq_b)
    gap_positions_a, gap_positions_b = extract_gap_positions(alignments[0])
    return gap_positions_a, gap_positions_b

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


def recreate_alignment(seq, gap_string):
    if not gap_string:
        return seq
    result = []
    gap_positions = set(map(int, gap_string.split(',')))
    for i, char in enumerate(seq):
        if i in gap_positions:
            result.append('-')
        result.append(char)
    return ''.join(result)

def recreate_alignment(seq, gap_string):
    if not gap_string:
        return seq
    result = []
    gap_positions = set(map(int, gap_string.split(',')))
    for i, char in enumerate(seq):
        if i in gap_positions:
            result.append('-')
        result.append(char)
    return ''.join(result)

def find_gap_strings(fasta_file):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    longest_seqs = {}
    gene_alignments = {}

    # Find the longest sequence for each gene
    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        seq = str(seq_record.seq)
        if gene_name not in longest_seqs or len(seq) > len(longest_seqs[gene_name][1]):
            longest_seqs[gene_name] = (header, seq)

    for item in longest_seqs:
        gene_alignments[longest_seqs[item][0]] = {}

    t = 0
    # Align each sequence to the longest one of the same gene
    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        seq = str(seq_record.seq)
        longest_seq_header, longest_seq = longest_seqs[gene_name]
        if header != longest_seq_header:
            gap_positions_a, gap_positions_b = align_sequences(longest_seq, seq)
            a_gaps = find_gap_positions(gap_positions_a)
            b_gaps = find_gap_positions(gap_positions_b)
            gene_alignments[longest_seq_header][header] = [a_gaps, b_gaps]


    return gene_alignments


gene_alignments = find_gap_strings("consensus_genes_2.fasta")

# Convert to JSON and save to a file
with open('gap_map.json', 'w') as outfile:
    json.dump(gene_alignments, outfile, indent=4)

print("Alignment data saved to gene_alignments.json")