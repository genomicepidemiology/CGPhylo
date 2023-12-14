from Bio import SeqIO, Align
from Bio.Align import PairwiseAligner

def align_sequences(fasta_file):
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    longest_seqs = {}

    # Find the longest sequence for each gene
    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        seq = str(seq_record.seq)
        if gene_name not in longest_seqs or len(seq) > len(longest_seqs[gene_name][1]):
            longest_seqs[gene_name] = (header, seq)

    aligner = PairwiseAligner()
    aligner.mode = 'local'  # Set aligner mode to local
    gene_alignments = {}

    # Align each sequence to the longest one of the same gene
    for header, seq_record in sequences.items():
        gene_name = header.split('_')[0]
        seq = str(seq_record.seq)
        longest_seq_header, longest_seq = longest_seqs[gene_name]

        if header != longest_seq_header:
            alignments = aligner.align(longest_seq, seq)
            best_alignment = alignments[0]
            cigar = create_cigar_string(best_alignment)

            # Print the alignment
            #print(f"Alignment between {longest_seq_header} and {header}:\n{best_alignment}")

            if gene_name not in gene_alignments:
                gene_alignments[gene_name] = {}
            gene_alignments[gene_name][header] = cigar

    return gene_alignments

def create_cigar_string(alignment):
    gicar = []  # GICAR string to store only gap information
    ref_pos, query_pos = 0, 0  # Initialize reference and query positions

    for aligned_pairs in zip(*alignment.aligned):
        for start, end in zip(*aligned_pairs):
            ref_gap = start - ref_pos
            query_gap = start - query_pos
            if ref_gap > 0:
                gicar.append(f"{ref_gap}D")  # Deletion in query
            if query_gap > 0:
                gicar.append(f"{query_gap}I")  # Insertion in target
            ref_pos = end
            query_pos = end

    # Check for trailing gaps
    if len(alignment.target) > ref_pos:
        gicar.append(f"{len(alignment.target) - ref_pos}D")
    if len(alignment.query) > query_pos:
        gicar.append(f"{len(alignment.query) - query_pos}I")

    return "".join(gicar)


# Usage example
gene_alignments = align_sequences("consensus_genes_2.fasta")

for gene in gene_alignments:
    print(f"Gene: {gene}")
    for allele in gene_alignments[gene]:
        print(f"  {allele}: {gene_alignments[gene][allele]}")
