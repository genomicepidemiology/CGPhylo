import os
import sys

def mintyper2_pipeline():
    #KMA ALIGnment
    gene_list = find_common_genes('/home/people/malhal/mintyper2/test/output_cpo_test')
    #Right now we ONLY use perfect length matches
    file_names, sequences_list = extract_sequences('/home/people/malhal/mintyper2/test/output_cpo_test', gene_list)
    distance_matrix = calculate_distance_matrix(sequences_list)
    print_distance_matrix(file_names, distance_matrix)
def print_distance_matrix(file_names, distance_matrix):
    # Find the maximum length of file names for formatting
    max_name_length = max(len(name) for name in file_names)

    # Print header row with file names
    header = " " * (max_name_length + 3)  # Initial spacing for header row
    for name in file_names:
        header += f"{name:>{max_name_length}}  "
    print(header)

    # Print each row of the distance matrix with the corresponding file name
    for i, row in enumerate(distance_matrix):
        row_str = f"{file_names[i]:<{max_name_length}}: "
        for dist in row:
            row_str += f"{dist:>{max_name_length}}  "
        print(row_str)
def calculate_distance_matrix(sequences):
    num_sequences = len(sequences)
    distance_matrix = [[0 for _ in range(num_sequences)] for _ in range(num_sequences)]

    # Compare each sequence with every other sequence
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            # Count differences in nucleotides at each position
            differences = sum(1 for a, b in zip(sequences[i], sequences[j]) if a != b)
            distance_matrix[i][j] = differences
            distance_matrix[j][i] = differences  # Symmetric matrix

    return distance_matrix

def extract_sequences(directory, headers):
    file_names = []
    sequences_list = []

    # Loop through each file in the directory
    for file in os.listdir(directory):
        if file.endswith('.fsa'):
            file_names.append(file.split('.')[0])
            sequences = []
            with open(os.path.join(directory, file), 'r') as fsa_file:
                lines = fsa_file.readlines()
                current_sequence = ''
                add_sequence = False

                # Loop through each line in the file
                for line in lines:
                    if line.startswith('>'):
                        if current_sequence and add_sequence:
                            sequences.append(current_sequence)
                        current_sequence = ''
                        add_sequence = line.split()[0][1:] in headers
                    elif add_sequence:
                        current_sequence += line.strip().upper()

                # Add the last found sequence
                if current_sequence and add_sequence:
                    sequences.append(current_sequence)

            sequences_list.append(sequences)

    return file_names, sequences_list

def load_fsa_gene_files(fsa_file):
    gene_dict = {}
    with open (fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                header = line.strip()[1:]
                gene_name = header.split('_')[0]
                gene_dict[gene_name] = ''
            else:
                gene_dict[gene_name] += line.strip()
    return gene_dict

def find_common_genes(directory_path):
    files = os.listdir(directory_path)

    gene_lists = []

    # Collect genes from each file
    for file in files:
        if file.endswith('.res'):
            genes = []
            with open(os.path.join(directory_path, file), 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        line = line.strip().split('\t')
                        genes.append(line[0])
            gene_lists.append(genes)

    # Find common genes
    if gene_lists:
        common = set(gene_lists[0])
        for gene_list in gene_lists[1:]:
            common.intersection_update(gene_list)
        return common
    else:
        return set()  # Return an empty set if no gene lists were found
