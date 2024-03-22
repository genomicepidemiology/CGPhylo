from ete3 import Tree

# Function to modify the taxa names in the tree to only include the ID
def modify_taxa_names(tree):
    for leaf in tree.iter_leaves():
        # Split the name to isolate the ID part (assumes a consistent naming format)
        leaf.name = leaf.name.split('/')[-1].split('.fastq')[0]
    return tree

# Path to your input Newick file
input_newick_path = '/home/people/malhal/papers/cgphylo/saffrontree/output/mintyper_illumina/kmer_tree.newick'
# Path to your output Newick file
output_newick_path = '/home/people/malhal/papers/cgphylo/saffrontree/output/mintyper_illumina/kmer_tree_taxa.newick'

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
