from Bio import Phylo
import numpy as np

# Load Newick file
newick_file = '/home/people/malhal/papers/cgphylo/MINTyper/results/mintyper_illumina_no_prune/tree_taxa.newick'

# Function to calculate distance matrix and sorted taxa labels
def calculate_distance_matrix(tree):
    # Sort taxa lexicographically by their names
    taxa = sorted(tree.get_terminals(), key=lambda x: x.name)
    n = len(taxa)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):  # Only fill upper triangle, it's symmetric
            distance = tree.distance(taxa[i], taxa[j])
            matrix[i, j] = matrix[j, i] = distance
    # Return the matrix and the sorted taxa names
    return matrix, [taxon.name for taxon in taxa]

# Function to print the distance matrix in the specified format
def print_formatted_matrix(matrix, taxa):
    print(f"{'':<20}", f"{len(taxa)}")
    for i, taxon in enumerate(taxa):
        print(f"{taxon:<20}", end='')
        for j in range(i):
            print(f"{int(matrix[i, j]):<20}", end='')
        print()  # New line at the end of each row

# Read the Newick tree
tree = Phylo.read(newick_file, 'newick')

# Calculate the distance matrix with taxa sorted lexicographically
distance_matrix, sorted_taxa_labels = calculate_distance_matrix(tree)

# Print the formatted distance matrix
print_formatted_matrix(distance_matrix, sorted_taxa_labels)
