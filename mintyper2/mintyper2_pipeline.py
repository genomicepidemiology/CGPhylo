import os
import sys
import json
import re


def mintyper2_pipeline(args):
    """Main function"""

    os.system('mkdir {}'.format(args.output))
    # Run KMA alignment for bacteria mapping
    # Run KMA alignment for cgMLST mapping
    # TBD Build function which excludes samples with the species.
    """
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
            species_db_string = get_species_db_string(top_template, args.db_dir)
            gap_map_path = species_db_string[:-5] + 'gap_map.json'

            cmd = 'kma -i {} -o {}/{} -t_db {} -ID 90 -md 5 -mct 0.5 -t 8 -mem_mode -dense -ref_fsa -ont'.format(file, args.output, name, species_db_string)
            os.system(cmd)
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            os.system('kma -i {} {} -o {}{} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                      .format(args.illumina[0], args.illumina[1], args.output + '/species_mapping_', name, args.db_dir + '/bac_db/bac_db',
                              args.threads))
            top_template = highest_scoring_hit_spa_file(args.output + '/species_mapping_' + name + '.spa')
            species_db_string = get_species_db_string(top_template, args.db_dir)
            gap_map_path = species_db_string[:-5] + 'gap_map.json'

            cmd = 'kma -i {} {} -o {}/{} -t_db {} -ID 90 -mct 0.5 -md 5 -mem_mode -dense -ref_fsa -t 8'.format(args.illumina[i], args.illumina[i+1], args.output, name, species_db_string)
            os.system(cmd)
    """
    gap_map_path = '/home/people/malhal/databases/cgmlst_dbs/cgmlst_db/Escherichia_coli_cgMLST_alleles/Escherichia_coli_cgMLST_alleles_consensus_gap_map.json'
    gene_list, non_shared_genes = find_common_genes(args.output)
    print (gene_list)
    sys.exit()
    print (len(gene_list))
    print (len(non_shared_genes))
    file_sequences_dict = load_sequences_from_file(args.output, gene_list)
    print (file_sequences_dict)
    sys.exit()
    gap_map = load_json(gap_map_path)
    print ('Loaded gap map')
    distance_matrix, file_names, total_length = calculate_pairwise_distances(file_sequences_dict, gap_map)
    print_distance_matrix_phylip(distance_matrix, file_names, args.output, total_length)
    print("The output distance matrix has been normalized to a genome size of 1.000.000. The identified core genes spanned {} bases.".format(total_length), file=sys.stderr)

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
                print(f"\t{round(int(row[j]) * normalization_factor, 2)}", end='', file=w)
            print(file=w)

def calculate_pairwise_distances(sequences_dict, gap_map):
    file_names = list(sequences_dict.keys())
    num_files = len(file_names)
    distance_matrix = [[0 for _ in range(num_files)] for _ in range(num_files)]
    total_length = 0
    print (len(sequences_dict))
    for file in sequences_dict:
        print (file)
        print (sequences_dict[file])
        for gene in sequences_dict[file]:
            print (gene)
            print (sequences_dict[file][gene][1])
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

def get_species_db_string(top_hit, db_dir):
    #Update these lists in cgMLST changes are made
    top_species = top_hit.split(' ')[1] + ' ' + top_hit.split(' ')[2]
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
        db_string = "{}_{}_cgMLST_alleles".format(top_hit.split(' ')[1], top_hit.split(' ')[2])

    print (db_string)
    db_string = db_dir + '/' + db_string + '/' + db_string + '_consensus_genes'

    if os.path.exists(db_string + '.name'):
        return db_string
    else:
        print (db_string)
        #TBD Do any acutally exist here, instead we want to exclude this sample and give a warning in the log
        sys.exit('No cgMLST database found for species: ' + top_species)

def extract_gene_name(gene_string):
    match = re.match(r'(.+?)_len_\d+', gene_string)
    if match:
        return match.group(1)
    else:
        return None  # or some error handling