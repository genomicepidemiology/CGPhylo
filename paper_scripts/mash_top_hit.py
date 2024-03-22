def find_highest_overlap(mash_results_file):
    highest_overlap_reference = None
    highest_overlap_count = -1
    highest_identity_score = -1.0

    with open(mash_results_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue  # Skip incomplete lines

            reference, _, identity_score, _, overlap = parts
            identity_score = float(identity_score)
            overlap_count = int(overlap.split('/')[0])

            # Update if this line has a higher overlap or equal overlap with a higher identity score
            if overlap_count > highest_overlap_count or (
                    overlap_count == highest_overlap_count and identity_score > highest_identity_score):
                highest_overlap_reference = reference
                highest_overlap_count = overlap_count
                highest_identity_score = identity_score

    return highest_overlap_reference


# Example usage
mash_results_file = 'mash_results'
highest_overlap_reference = find_highest_overlap(mash_results_file)
print(f'The reference with the highest overlap is: {highest_overlap_reference}')
