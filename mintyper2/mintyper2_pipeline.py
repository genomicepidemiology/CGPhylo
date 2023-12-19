import os
import sys

def mintyper2_pipeline(args):
    """Main function"""
    os.system('mkdir {}'.format(args.output))
    #Find species
    """
    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]
            cmd = 'kma -i {} -o {}/{} -t_db /home/people/malhal/mintyper2/consensus_genes_db -ID 90 -md 5 -mct 0.5 -t 8 -mem_mode -dense -ref_fsa -ont'.format(file, args.output, name)
            os.system(cmd)
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            cmd = 'kma -i {} {} -o {}/{} -t_db /home/people/malhal/mintyper2/consensus_genes_db -ID 90 -mct 0.5 -md 5 -mem_mode -dense -ref_fsa -t 8'.format(args.illumina[i], args.illumina[i+1], args.output, name)
            os.system(cmd)
    """
    gene_list, non_shared_genes = find_common_genes(args.output)

    #Distance non_shared_genes
    #Figure how many of the shared genes have same size.
    print (len(gene_list), 'genes shared between all samples')
    print (len(non_shared_genes), 'genes not shared between all sam ples')
    same_length_genes, genes_to_readjust = find_common_genes_with_same_length(args.output, gene_list)
    print (len(same_length_genes))
    print ('genes to fix', len(genes_to_readjust))
    find_lengths_of_genes_to_readjust(args.output, genes_to_readjust)
    #genes_to_readjust holds the identifier for the genes that need to be readjusted. Look up the top scorer and realign.
    sys.exit()
    for item in non_shared_genes[0:20]:
        print (item)
    sequences_dict = extract_sequences(args.output, gene_list)
    for key in sequences_dict:
        print(key, len(sequences_dict[key]))
    distance_matrix, file_names = calculate_pairwise_distances(sequences_dict)
    print_distance_matrix_phylip(distance_matrix, file_names, args.output)

def find_lengths_of_genes_to_readjust(output, genes_to_readjust):
    gene_dict = dict()
    files = os.listdir(output)
    for file in files:
        if file.endswith('.res'):
            with open(os.path.join(output, file), 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        line = line.strip().split('\t')
                        allele = line[0].split('_')[:-2]
                        gene = '_'.join(allele)
                        if gene in genes_to_readjust:
                            if gene not in gene_dict:
                                gene_dict[gene] = [[], []]
                                gene_dict[gene][0].append(line[0])
                                gene_dict[gene][1].append(file)
                            else:
                                gene_dict[gene][0].append(line[0])
                                gene_dict[gene][1].append(file)

    for gene in gene_dict:
        print (gene, gene_dict[gene])
    sys.exit()
def recreate_alignment(seq, gap_string):
    if not gap_string:
        return seq
    result = []
    gap_positions = set(map(int, gap_string.split(',')))
    for i, char in enumerate(seq):
        if i in gap_positions:
            result.append('-')
        result.append(char)
    return ''.join(result)


def find_common_genes_with_same_length(output, gene_list):
    same_length_list = list()
    files = os.listdir(output)


    # Collect genes from each file
    for file in files:
        if file.endswith('.res'):
            genes = set()
            top_score_dict = dict()
            with open(os.path.join(output, file), 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        line = line.strip().split('\t')
                        allele = line[0].split('_')[:-2]
                        gene = '_'.join(allele)
                        if gene in gene_list:
                            if gene not in top_score_dict:
                                top_score_dict[gene] = [line[0], line[1]]
                            elif float(line[1]) > float(top_score_dict[gene][1]):
                                top_score_dict[gene] = [line[0], line[1]]

            for gene in top_score_dict:
                genes.add(top_score_dict[gene][0])
            same_length_list.append(genes)

    # Find common genes
    common = set(same_length_list[0])
    for gene_list in same_length_list[1:]:
        common.intersection_update(gene_list)

    genes_to_reajust = set()
    for gene_list in same_length_list:
        for allele in gene_list:
            if allele not in common:
                gene = '_'.join(allele.split('_')[:-2])
                genes_to_reajust.add(gene)

    return list(common), list(genes_to_reajust)


def print_distance_matrix_phylip(distance_matrix, file_names, output):
    num_files = len(file_names)
    with open(output + '/distance_matrix.phylip', 'w') as w:


        # Print the number of files first
        print(f"{num_files}", file=w)

        # Print each row of the distance matrix in the specified format
        for i, row in enumerate(distance_matrix):
            # Start with the file name
            print(file_names[i], end='', file=w)

            # Print the distances for the lower triangular matrix
            for j in range(i):
                print(f"\t{row[j]}", end='', file=w)
            print(file=w)

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
            genes = set()
            with open(os.path.join(directory_path, file), 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        line = line.strip().split('\t')
                        allele = line[0].split('_')[:-2]
                        gene = '_'.join(allele)
                        genes.add(gene)

            gene_lists.append(genes)
    for item in gene_lists:
        print (len(item))

    # Find common genes
    common = set(gene_lists[0])
    for gene_list in gene_lists[1:]:
        common.intersection_update(gene_list)

    non_shared = set()
    for gene_list in gene_lists:
        for gene in gene_list:
            if gene not in common:
                non_shared.add(gene)

    return common, non_shared


def find_duplicates(strings):
    counts = {}
    duplicates = []

    for string in strings:
        if string in counts:
            counts[string] += 1
        else:
            counts[string] = 1

    for string, count in counts.items():
        if count > 1:
            duplicates.append(string)

    return duplicates