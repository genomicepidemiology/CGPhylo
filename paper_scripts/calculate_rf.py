from ete3 import Tree

# Function to read a tree from a Newick file
def read_tree_from_file(file_path):
    with open(file_path, 'r') as file:
        return Tree(file.read().strip())

# Paths to the Newick files
newick_file_1 = '/home/people/malhal/papers/cgphylo/MINTyper/results/mintyper_illumina_no_prune/tree_taxa.newick'
newick_file_2 = '/home/people/malhal/papers/cgphylo/saffrontree/output/mintyper_illumina/kmer_tree_taxa.newick'


# Load the trees from the Newick files
tree1 = read_tree_from_file(newick_file_1)
tree2 = read_tree_from_file(newick_file_2)

# Use the Robinson-Foulds comparison and capture all returned values in a variable
rf_results = tree1.robinson_foulds(tree2, unrooted_trees=True)

# Extract the RF distance and maximum possible RF distance from the results
rf_distance, rf_max = rf_results[0], rf_results[1]

# Print the Robinson-Foulds distance and the maximum possible RF distance
print(f"Robinson-Foulds Distance: {rf_distance}")
print(f"Maximum Possible RF Distance: {rf_max}")

print(tree1)
print(tree2)