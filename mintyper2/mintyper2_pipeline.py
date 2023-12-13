import os
import sys

def mintyper2_pipeline(args):
    """Main function"""
    os.system('mkdir {}'.format(args.output))
    #Find species
    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]
            cmd = 'kma -i {} -o {}/{} -t_db /home/people/malhal/mintyper2/consensus_genes_db -ID 50 -md 5 -mct 0.5 -t 8 -mem_mode -ref_fsa -ont'.format(file, args.output, name)
            os.system(cmd)
    if args.illumina != []:
        #Look into if we can parallelize this.
        #Could we run multiple KMA alignments at once which 2 threads.
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            cmd = 'kma -i {} {} -o {}/{} -t_db /home/people/malhal/mintyper2/consensus_genes_db -ID 50 -mct 0.5 -md 5 -mem_mode -dense -ref_fsa -t 8'.format(args.illumina[i], args.illumina[i+1], args.output, name)
            os.system(cmd)
    #KMA ALIGnment
    gene_list = find_common_genes(args.output)
    sequences_dict = extract_sequences(args.output, gene_list)
    for key in sequences_dict:
        print(key, len(sequences_dict[key]))
    #Right now we ONLY use perfect length matches
    distance_matrix, file_names = calculate_pairwise_distances(sequences_dict)
    print_distance_matrix_phylip(distance_matrix, file_names)


def print_distance_matrix_phylip(distance_matrix, file_names):
    num_files = len(file_names)

    # Print the number of files first
    print(f"{num_files}")

    # Print each row of the distance matrix in the specified format
    for i, row in enumerate(distance_matrix):
        # Start with the file name
        print(file_names[i], end='')

        # Print the distances for the lower triangular matrix
        for j in range(i):
            print(f"\t{row[j]}", end='')
        print()

def calculate_pairwise_distances(sequences_dict):
    file_names = list(sequences_dict.keys())
    num_files = len(file_names)
    distance_matrix = [[0 for _ in range(num_files)] for _ in range(num_files)]

    # Iterate over each pair of files
    for i in range(num_files):
        for j in range(i + 1, num_files):
            count = 0  # Count of differences

            # Compare each gene's sequences nucleotide by nucleotide
            for gene in sequences_dict[file_names[i]].keys():
                seq1 = sequences_dict[file_names[i]][gene]
                seq2 = sequences_dict[file_names[j]][gene]
                if len(seq1) != len(seq2):
                    print (f"Warning: {gene} has different lengths between {file_names[i]} and {file_names[j]}")

                # Modified comparison to skip lowercase nucleotides
                diff = sum(1 for a, b in zip(seq1, seq2) if a != b and not (a.islower() or b.islower()))

                # Count differences
                count += diff

            # Store the count in the matrix
            distance_matrix[i][j] = count
            distance_matrix[j][i] = count  # Symmetric matrix

    return distance_matrix, file_names

def extract_sequences(directory, headers):
    sequences_dict = {}

    # Loop through each file in the directory
    for file in os.listdir(directory):
        if file.endswith('.fsa'):
            file_name = file.split('.')[0]
            sequences_dict[file_name] = {header: '' for header in headers}

            with open(os.path.join(directory, file), 'r') as fsa_file:
                lines = fsa_file.readlines()
                current_sequence = ''
                current_header = None

                # Loop through each line in the file
                for line in lines:
                    if line.startswith('>'):
                        if current_sequence and current_header:
                            sequences_dict[file_name][current_header] = current_sequence
                        current_header = line.split()[0][1:]
                        current_sequence = ''
                    elif current_header in headers:
                        current_sequence += line.strip()

                # Add the last found sequence
                if current_sequence and current_header:
                    sequences_dict[file_name][current_header] = current_sequence

    return sequences_dict

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
