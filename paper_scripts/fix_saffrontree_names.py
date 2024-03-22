import re

def shorten_paths_in_newick(newick_file_path, output_file_path):
    # Read the Newick file
    with open(newick_file_path, 'r') as file:
        newick_string = file.read()

    # Use regular expression to find all paths and replace them with just the file names
    # This regex matches the path and captures the file name for replacement
    #shortened_newick_string = re.sub(r'/home/people/malhal/data/cgphylo/mintyper_data/illumina/', r'', newick_string)
    shortened_newick_string = re.sub(r'_q10', r'', newick_string)

    # Write the modified Newick string to a new file
    with open(output_file_path, 'w') as file:
        file.write(shortened_newick_string)

# Example usage
newick_file_path = 'tree'
output_file_path = 'tree'

# Call the function with your file paths
shorten_paths_in_newick(newick_file_path, output_file_path)
