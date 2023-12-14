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
            print(f"Alignment between {longest_seq_header} and {header}:\n{best_alignment}")

            if gene_name not in gene_alignments:
                gene_alignments[gene_name] = {}
            gene_alignments[gene_name][header] = cigar

    return gene_alignments

def create_cigar_string(alignment):
    cigar = []
    match_count = 0
    last_end = 0

    for aligned_pairs in zip(*alignment.aligned):
        for start, end in zip(*aligned_pairs):
            if start > last_end:
                cigar.append(f"{start - last_end}I")  # Insertion
            cigar.append(f"{end - start}M")  # Match/Mismatch
            last_end = end

    # Check for trailing insertions
    if len(alignment.target) > last_end:
        cigar.append(f"{len(alignment.target) - last_end}I")

    return "".join(cigar)

# Usage example
gene_alignments = align_sequences("your_fasta_file.fasta")

for gene in gene_alignments:
    print(f"Gene: {gene}")
    for allele in gene_alignments[gene]:
        print(f"  {allele}: {gene_alignments[gene][allele]}")
