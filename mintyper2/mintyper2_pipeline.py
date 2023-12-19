import os
import sys
import json
def mintyper2_pipeline(args):
    """Main function"""
    os.system('mkdir {}'.format(args.output))
    #Find species and load database
    """
    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]
            cmd = 'kma -i {} -o {}/{} -t_db /home/people/malhal/mintyper2/consensus_genes_db_2 -ID 90 -md 5 -mct 0.5 -t 8 -mem_mode -dense -ref_fsa -ont'.format(file, args.output, name)
            os.system(cmd)
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            cmd = 'kma -i {} {} -o {}/{} -t_db /home/people/malhal/mintyper2/consensus_genes_db_2 -ID 90 -mct 0.5 -md 5 -mem_mode -dense -ref_fsa -t 8'.format(args.illumina[i], args.illumina[i+1], args.output, name)
            os.system(cmd)
    """
    gene_list, non_shared_genes = find_common_genes(args.output)
    file_sequences_dict = load_sequences_from_file(args.output, gene_list)
    file_path = '/home/people/malhal/mintyper2/gap_map.json'
    gap_map = load_json(file_path)
    distance_matrix, file_names, total_length = calculate_pairwise_distances(file_sequences_dict, gap_map)
    print_distance_matrix_phylip(distance_matrix, file_names, args.output, total_length)

def load_json(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

def load_sequences_from_file(output, gene_list):
    file_sequences_dict = dict()
    files = os.listdir(output)
    for file in files:
        if file.endswith('.fsa'):
            name = file.split('.')[0]
            file_sequences_dict[name] = dict()
            with open(os.path.join(output, file), 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        line = line.strip()
                        allele = line[1:]
                        gene = line[1:].split('_')[0]
                        if gene in gene_list:
                            file_sequences_dict[name][gene] = [allele, '']
                    if gene != None and gene in gene_list and not line.startswith('>'):
                        file_sequences_dict[name][gene][1] += line.strip()
    return file_sequences_dict



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


def print_distance_matrix_phylip(distance_matrix, file_names, output, total_length):
    normalization_factor = 1000000 / total_length
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
                print(f"\t{int(row[j])*normalization_factor}", end='', file=w)
            print(file=w)

def calculate_pairwise_distances(sequences_dict, gap_map):
    file_names = list(sequences_dict.keys())
    num_files = len(file_names)
    distance_matrix = [[0 for _ in range(num_files)] for _ in range(num_files)]
    total_length = 0
    for file in sequences_dict:
        for gene in sequences_dict[file]:
            total_length += len(sequences_dict[file][gene][1])
        break # Only need to do this once
    # Iterate over each pair of files
    for i in range(num_files):
        for j in range(i + 1, num_files):
            count = 0  # Count of differences

            # Compare each gene's sequences nucleotide by nucleotide
            for gene in sequences_dict[file_names[i]].keys():
                seq1 = sequences_dict[file_names[i]][gene][1]  # Get the sequence for the first file
                seq2 = sequences_dict[file_names[j]][gene][1]  # Get the sequence for the second file

                # Retrieve gap strings from the gap_map
                if len(seq1) != len(seq2):
                    allele_1 = sequences_dict[file_names[i]][gene][0]
                    allele_2 = sequences_dict[file_names[j]][gene][0]
                    gap_string1, gap_string2 = gap_map[allele_1][allele_2]
                    # Realign sequences
                    realigned_seq1 = recreate_alignment(seq1, gap_string1)
                    realigned_seq2 = recreate_alignment(seq2, gap_string2)
                else:
                    realigned_seq1 = seq1
                    realigned_seq2 = seq2

                # Ensure the sequences are of the same length after realignment
                if len(realigned_seq1) != len(realigned_seq2):
                    print(f"Warning: {gene} has different lengths between {file_names[i]} and {file_names[j]} after realignment")

                # Modified comparison to skip lowercase nucleotides
                diff = sum(1 for a, b in zip(realigned_seq1, realigned_seq2) if a != b and not (a.islower() or b.islower()))
                # Count differences
                count += diff

            # Store the count in the matrix
            distance_matrix[i][j] = count
            distance_matrix[j][i] = count  # Symmetric matrix

    return distance_matrix, file_names, total_length

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