import os
import sys
import json
import re
import hashlib
import time


def mintyper2_pipeline(args):
    """Main function"""

    os.system('mkdir {}'.format(args.output))

    # Check all species

    #exclude_list, top_specie = check_all_species(args)
    top_specie = 'Salmonella enterica'
    species_db_string = get_species_db_string(top_specie, args.db_dir)
    genome_size = find_highest_length_in_spa_files(args.output, top_specie)
    gap_map_path = species_db_string[:-5] + 'gap_map.json'

    """
    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]
            if not name in exclude_list:
                cmd = 'kma -i {} -o {}/{} -t_db {} -ID 90 -md 5 -mct 0.5 -t 8 -mem_mode -dense -ref_fsa -ont'.format(file, args.output, name, species_db_string)
                os.system(cmd)
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            if not name in exclude_list:
                cmd = 'kma -i {} {} -o {}/{} -t_db {} -ID 90 -mct 0.5 -md 5 -mem_mode -dense -ref_fsa -t 8'.format(args.illumina[i], args.illumina[i+1], args.output, name, species_db_string)
                os.system(cmd)
    """
    #gap_map_path = '/home/people/malhal/databases/cgmlst_dbs/cgmlst_db/Escherichia_coli_cgMLST_alleles/Escherichia_coli_cgMLST_alleles_consensus_gap_map.json'
    gene_list, non_shared_genes = find_common_genes(args.output)
    print (len(gene_list), 'genes found in all samples (core genes)')
    print (len(non_shared_genes), 'genes not found in all samples (non-shared genes)')
    file_sequences_dict = load_sequences_from_file(args.output, gene_list)
    gap_map = load_json(gap_map_path)
    t1 = time.time()
    distance_matrix, file_names = calculate_pairwise_distances(file_sequences_dict, gap_map)
    t2 = time.time()
    print ("Time to calculate distance matrix: ", t2-t1)
    normalization_factor = 1000000 / genome_size
    distance_matrix_output_name = 'distance_matrix_1M.txt'
    print_distance_matrix_phylip(distance_matrix, file_names, args.output, distance_matrix_output_name, normalization_factor)
    print("A distance matrix normalized to a genome size of 1.000.000 has been outputted. The identified core genes spanned {} bases.".format(genome_size), file=sys.stderr)
    #TBD should we give an option to give input for normalization factor? Genome size?
    #distance_matrix_output_name = 'distance_matrix_GS.txt'
    #print_distance_matrix_phylip(distance_matrix, file_names, args.output, distance_matrix_output_name, 1)
    #print ("A distance matrix normalized to a genome size of {} has been outputted. The identified core genes spanned {} bases.".format(genome_size, genome_size), file=sys.stderr)


def check_all_species(args):
    top_template_count = dict()
    reference_results = dict()

    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]
            os.system('kma -i {} -o {}{} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                      .format(file, args.output + '/species_mapping_', name, args.db_dir + '/bac_db/bac_db',
                              args.threads))
            top_template = highest_scoring_hit_spa_file(args.output + '/species_mapping_' + name + '.spa')
            specie = top_template.split(' ')[1] + ' ' + top_template.split(' ')[2]
            if specie in top_template_count:
                top_template_count[specie] += 1
            else:
                top_template_count[specie] = 1
            reference_results[name] = specie

    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            os.system('kma -i {} {} -o {}{} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                      .format(args.illumina[i], args.illumina[i+1], args.output + '/species_mapping_', name, args.db_dir + '/bac_db/bac_db',
                              args.threads))
            top_template = highest_scoring_hit_spa_file(args.output + '/species_mapping_' + name + '.spa')
            specie = top_template.split(' ')[1] + ' ' + top_template.split(' ')[2]
            if specie in top_template_count:
                top_template_count[specie] += 1
            else:
                top_template_count[specie] = 1
            reference_results[name] = specie

    top_specie = max(top_template_count, key=top_template_count.get)
    print('The most common specie is {} with {} hits.'.format(top_specie, top_template_count[top_specie]))

    exclude_list = []

    for file in reference_results:
        if reference_results[file] != top_specie:
            exclude_list.append(file)

    return exclude_list, top_specie


def find_highest_length_in_spa_files(directory, species):
    highest_length = 0

    # Check if the directory exists
    if not os.path.exists(directory):
        print(f"Directory '{directory}' does not exist.")
        return None

    # Loop through all files in the directory
    for filename in os.listdir(directory):
        if filename.endswith(".spa"):
            filepath = os.path.join(directory, filename)

            # Open the .spa file
            with open(filepath, 'r') as spa_file:
                # Read each line in the file
                for line in spa_file:
                    # Split the line by tab to access the columns
                    columns = line.strip().split('\t')

                    # Check if the species is a substring in the line
                    if species in columns[0]:
                        # Extract the length from the appropriate column
                        length = int(columns[4])

                        # Update the highest length if necessary
                        if length > highest_length:
                            highest_length = length

    return highest_length




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
                        gene = extract_gene_name(allele)
                        if gene in gene_list:
                            file_sequences_dict[name][gene] = [allele, '']
                    if gene != None and gene in gene_list and not line.startswith('>'):
                        file_sequences_dict[name][gene][1] += line.strip()
    return file_sequences_dict



def recreate_alignment(seq, gap_string):
    if not gap_string:
        return seq

    gap_positions = list(map(int, gap_string.split(',')))
    # Sort gap positions to ensure we insert them in the correct order
    gap_positions.sort()

    result = []
    added_gaps = 0

    for i, char in enumerate(seq):
        # Adjust the position for gaps already added
        while added_gaps < len(gap_positions) and i + added_gaps == gap_positions[added_gaps]:
            result.append('-')
            added_gaps += 1
        result.append(char)

    # Append remaining gaps if any (in case they are at the end of the sequence)
    result.extend('-' * (len(gap_positions) - added_gaps))

    return ''.join(result)



def print_distance_matrix_phylip(distance_matrix, file_names, output, distance_matrix_output_name, normalization_factor):
    num_files = len(file_names)
    with open(output + '/' + distance_matrix_output_name, 'w') as w:


        # Print the number of files first
        print(f"{num_files}", file=w)

        # Print each row of the distance matrix in the specified format
        for i, row in enumerate(distance_matrix):
            # Start with the file name
            print(file_names[i], end='', file=w)

            # Print the distances for the lower triangular matrix
            for j in range(i):
                print(f"\t{round(int(row[j]) * normalization_factor, 2)}", end='', file=w)
            print(file=w)

def calculate_pairwise_distances(sequences_dict, gap_map):
    file_names = list(sequences_dict.keys())
    num_files = len(file_names)
    distance_matrix = [[0 for _ in range(num_files)] for _ in range(num_files)]
    # Iterate over each pair of files
    # TBD consider if this can be done faster
    # TBD consider if we can use itertools.combinations to avoid comparing the same files twice
    for i in range(num_files):
        for j in range(i + 1, num_files):
            #TBD write a speed compare test. Basically, it is quicker to hash both sequences if they have equal length and then only check the SNPs if hashes don't match
            print (f"Comparing {file_names[i]} and {file_names[j]}")
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


                # Dont count gaps test
                #diff = sum(1 for a, b in zip(realigned_seq1, realigned_seq2) if
                #           a != b and a != '-' and b != '-' and not (a.islower() or b.islower()))

                if hash(realigned_seq1) != hash(realigned_seq2):
                    diff = sum(1 for a, b in zip(realigned_seq1, realigned_seq2) if
                               a != b and a != '-' and b != '-' and not (a.islower() or b.islower()))
                else:
                    diff = 0
                #print (diff)

                # Counts gaps. Gaps should not be included in SNPs distances, but consider using this for a another metric in the future.
                #diff = sum(1 for a, b in zip(realigned_seq1, realigned_seq2) if
                #           a != b and ((a == '-' or b == '-') or not (a.islower() or b.islower())))

                #if file_names[i] == 'SRR1188432_1':
                #    if file_names[j] == 'SRR1188445_1':
                #        gene_hash_1 = hashlib.md5(realigned_seq1.encode()).hexdigest()
                #        gene_hash_2 = hashlib.md5(realigned_seq2.encode()).hexdigest()
                #        print (gene_hash_1, gene_hash_2)
                #if diff > 0:
                #    if file_names[i] == 'SRR1188432_1':
                #        if file_names[j] == 'SRR1188445_1':
                #            print(f"{gene} has {diff} differences between {file_names[i]} and {file_names[j]}")
                #            print (len(seq1), len(seq2))
                #            print (realigned_seq1)
                #            print (realigned_seq2)
                #            print ('seq1: ', seq1)
                #            print ('seq2: ', seq2)

                # Count differences
                count += diff

            # Store the count in the matrix
            distance_matrix[i][j] = count
            distance_matrix[j][i] = count  # Symmetric matrix

    return distance_matrix, file_names

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
                        gene = extract_gene_name(line[0])
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


def highest_scoring_hit_spa_file(file_path):
    """
    Identifies and returns the highest scoring template from a tab-separated file.

    This function reads through a specified file, assuming it to be tab-separated,
    and identifies the row with the highest score in a designated score column.
    It returns the template corresponding to this highest score.

    Args:
        file_path (str): The path to the file to be read. Assumes a specific format where
                         the score is in the third column and the template in the first.

    Returns:
        str: The identifier of the highest scoring template.
    """

    highest_score = 0
    highest_scoring_template = ""

    with open(file_path, 'r') as file:
        next(file)  # Skip the header line
        for line in file:
            columns = line.split('\t')
            try:
                # Extract score and compare to find the highest
                score = int(columns[2])  # Score is expected in the 3rd column
                if score > highest_score:
                    highest_score = score
                    highest_scoring_template = columns[0]  # Template is expected in the 1st column
            except ValueError:
                # Skip line if score is not an integer or line is malformed
                continue

    return highest_scoring_template

def get_species_db_string(top_species, db_dir):
    #Update these lists in cgMLST changes are made
    Mycobacterium_list = ['Mycobacterium tuberculosis', 'Mycobacterium bovis', 'Mycobacterium aafricanum', 'Mycobacterium canettii']
    Klebsiella_list = ['Klebsiella pneumoniae', 'Klebsiella variicola', 'Klebsiella quasipneumoniae']
    Cronobacter_list = ['Cronobacter sakazakii', 'Cronobacter malonaticus']
    Campylobacter_list = ['Campylobacter jejuni', 'Campylobacter coli']
    if top_species in Mycobacterium_list:
        db_string = 'Mycobacterium_tuberculosis_bovis_africanum_canettii_cgMLST_alleles'
    elif top_species in Klebsiella_list:
        db_string = 'Klebsiella_pneumoniae_variicola_quasipneumoniae_cgMLST_alleles'
    elif top_species in Cronobacter_list:
        db_string = 'Cronobacter_sakazakii_malonaticus_cgMLST_alleles'
    elif top_species in Campylobacter_list:
        db_string = 'Campylobacter_jejuni_coli_cgMLST_alleles'
    else:
        db_string = "{}_{}_cgMLST_alleles".format(top_species.split(' ')[0], top_species.split(' ')[1])

    db_string = db_dir + '/' + db_string + '/' + db_string + '_consensus_genes'

    if os.path.exists(db_string + '.name'):
        return db_string
    else:
        #TBD Do any acutally exist here, instead we want to exclude this sample and give a warning in the log
        sys.exit('No cgMLST database found for species: ' + top_species)

def extract_gene_name(gene_string):
    match = re.match(r'(.+?)_len_\d+', gene_string)
    if match:
        return match.group(1)
    else:
        return None  # or some error handling