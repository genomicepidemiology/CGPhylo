import os
import sys
from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner
import re


def align_sequences(fasta_file):
    # Step 1: Parse the FASTA file
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Step 2: Identify the longest sequence for each gene
    longest_seqs = {}
    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        seq = str(seq_record.seq)  # Convert Seq object to string
        if gene_name not in longest_seqs or len(seq) > len(longest_seqs[gene_name][1]):
            longest_seqs[gene_name] = (header, seq)

    # Step 3: Align each sequence to the longest one of the same gene
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    gene_alignments = {}

    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        seq = str(seq_record.seq)  # Convert Seq object to string
        longest_seq_header, longest_seq = longest_seqs[gene_name]

        if header != longest_seq_header:
            alignment = aligner.align(longest_seq, seq)
            print (alignment)
            # Step 4: Generate CIGAR string
            cigar = alignment[0].cigar
            # Step 5: Store in a nested dictionary
            if gene_name not in gene_alignments:
                gene_alignments[gene_name] = {}
            gene_alignments[gene_name][header] = cigar

    return gene_alignments

# Usage
gene_alignments = align_sequences("consensus_genes.fasta")

for item in gene_alignments:
    print(item)
    for item2 in gene_alignments[item]:
        print(item2, gene_alignments[item][item2])
    print('\n')