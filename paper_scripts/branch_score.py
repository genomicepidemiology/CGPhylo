from ete3 import Tree


def read_tree_from_file(file_path):
    with open(file_path, 'r') as file:
        return Tree(file.read().strip())


def calculate_branch_score_distance(tree1, tree2):
    # Assuming both trees have a compatible structure for direct branch length comparison

    # Flatten the list of distances for each tree
    distances_tree1 = [node.dist for node in tree1.traverse() if not node.is_root()]
    distances_tree2 = [node.dist for node in tree2.traverse() if not node.is_root()]

    # Calculate the sum of absolute differences in distances between the trees
    # Ensuring both trees have the same number of nodes for a valid comparison
    if len(distances_tree1) != len(distances_tree2):
        raise ValueError("The trees have a different number of nodes; cannot calculate Branch Score Distance.")

    branch_score_distance = sum(abs(d1 - d2) for d1, d2 in zip(distances_tree1, distances_tree2))
    return branch_score_distance


# Paths to the Newick files
newick_file_1 = '/home/people/malhal/papers/cgphylo/MINTyper/results/mintyper_illumina_no_prune/tree_taxa.newick'
newick_file_2 = '/home/people/malhal/papers/cgphylo/saffrontree/output/mintyper_illumina/kmer_tree_taxa.newick'

# Load the trees from the Newick files
tree1 = read_tree_from_file(newick_file_1)
tree2 = read_tree_from_file(newick_file_2)

# Calculate the Branch Score Distance
try:
    bsd = calculate_branch_score_distance(tree1, tree2)
    print(f"Branch Score Distance: {bsd}")
except ValueError as e:
    print(e)
