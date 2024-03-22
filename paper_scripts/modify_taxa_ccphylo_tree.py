from ete3 import Tree

# Function to modify the taxa names in the tree
def modify_taxa_names(tree):
    for leaf in tree.iter_leaves():
        leaf.name = leaf.name.replace('.fsa', '')
    return tree

# Path to your input Newick file
input_newick_path = '/home/people/malhal/papers/cgphylo/MINTyper/results/mintyper_illumina_no_prune/tree.newick'
# Path to your output Newick file
output_newick_path = '/home/people/malhal/papers/cgphylo/MINTyper/results/mintyper_illumina_no_prune/tree_taxa.newick'

# Load the tree from the Newick file
with open(input_newick_path, 'r') as file:
    newick_string = file.read().strip()
    tree = Tree(newick_string)

# Modify the taxa names
modified_tree = modify_taxa_names(tree)

# Write the modified tree to a new Newick file
with open(output_newick_path, 'w') as output_file:
    output_file.write(modified_tree.write(format=1))

print("Modified Newick file has been written to", output_newick_path)
