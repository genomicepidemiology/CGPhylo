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
    aligner.mode = 'global'
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

            if gene_name not in gene_alignments:
                gene_alignments[gene_name] = {}
            gene_alignments[gene_name][header] = cigar

    return gene_alignments

def create_cigar_string(alignment):
    cigar = []
    match_count = 0
    gap_count = 0
    for start, end in zip(*alignment.aligned):
        if start - match_count > 0:
            if gap_count:
                cigar.append(f"{gap_count}I")
                gap_count = 0
            cigar.append(f"{start - match_count}M")
        gap_count += end - start
        match_count = end

    if gap_count:
        cigar.append(f"{gap_count}I")
    if len(alignment.query) > match_count:
        cigar.append(f"{len(alignment.query) - match_count}M")

    return "".join(cigar)

# Usage
# gene_alignments = align_sequences("path_to_your_fasta_file.fasta")

# Usage
gene_alignments = align_sequences("consensus_genes.fasta")

for item in gene_alignments:
    print(item)
    for item2 in gene_alignments[item]:
        print(item2, gene_alignments[item][item2])
    print('\n')