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