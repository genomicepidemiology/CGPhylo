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
    for header, seq in sequences.items():
        gene_name = header.split('_')[0]
        if gene_name not in longest_seqs or len(seq) > len(longest_seqs[gene_name]):
            longest_seqs[gene_name] = seq

    for item in longest_seqs:
        print (item, len(longest_seqs[item]))
    sys.exit()

    # Step 3: Align each sequence to the longest one of the same gene
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    gene_alignments = {}

    for header, seq in sequences.items():
        gene_name = header.split('_')[0]
        if header != gene_name + "_len_" + str(len(longest_seqs[gene_name])):
            alignment = aligner.align(longest_seqs[gene_name], seq)
            # Step 4: Generate CIGAR string
            cigar = alignment[0].cigar
            # Step 5: Store in a nested dictionary
            if gene_name not in gene_alignments:
                gene_alignments[gene_name] = {}
            gene_alignments[gene_name][header] = cigar

    return gene_alignments

# Usage
gene_alignments = align_sequences("consensus_genes.fasta")