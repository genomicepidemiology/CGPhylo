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

def find_gap_strings(fasta_file):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    gene_sequences = {}

    # Group sequences by gene name
    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        if gene_name not in gene_sequences:
            gene_sequences[gene_name] = []
        gene_sequences[gene_name].append((header, str(seq_record.seq)))

    print (len(gene_sequences))
    sys.exit()

    gene_alignments = {}

    # Align each sequence with every other sequence of the same gene
    for gene_name, gene_seqs in gene_sequences.items():
        for header_a, seq_a in gene_seqs:
            for header_b, seq_b in gene_seqs:
                if header_a != header_b:
                    if header_a not in gene_alignments:
                        gene_alignments[header_a] = {}
                    gap_positions_a, gap_positions_b = align_sequences(seq_a, seq_b)
                    a_gaps = find_gap_positions(gap_positions_a)
                    b_gaps = find_gap_positions(gap_positions_b)
                    gene_alignments[header_a][header_b] = [a_gaps, b_gaps]

    return gene_alignments

gene_alignments = find_gap_strings("consensus_genes_2.fasta")

# Convert to JSON and save to a file
with open('gap_map.json', 'w') as outfile:
    json.dump(gene_alignments, outfile, indent=4)

print("Alignment data saved to gene_alignments.json")