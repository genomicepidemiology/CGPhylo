import os
import sys
import json
import re
import hashlib
import time
import logging



def mintyper2_pipeline(args):
    """Main function"""

    os.system('mkdir {}'.format(args.output))

    logging.basicConfig(
        format='%(asctime)s %(message)s',
        filename=args.output + '/cgphylo.log',
        level=logging.INFO)

    if args.fsf:
        exclude_list, top_specie = fast_species_finder(args)
    else:
        exclude_list, top_specie = check_all_species(args)

    species_db_string = get_species_db_string(top_specie, args.db_dir)

    gap_map_path = species_db_string[:-5] + 'gap_map.json'
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
    outliers, non_outliers = find_gene_count_outliers(args.output)
    if len(outliers) > 0:
        logging.info('Outliers: {}. These samples failed to identify enough genes to be included in the analysis.'.format(outliers))
    gene_list, non_shared_genes = find_common_genes(args.output, outliers)
    logging.info('{} genes found in all samples (core genes)'.format(len(gene_list)))
    logging.info('{} genes not found in all samples (non-shared genes)'.format(len(non_shared_genes)))
    file_sequences_dict, cg_nucleotide_count = load_sequences_from_file(args.output, gene_list, outliers)
    logging.info('The core genes spanned {} bases.'.format(cg_nucleotide_count))
    gap_map = load_json(gap_map_path)
    distance_matrix, file_names = calculate_pairwise_distances(file_sequences_dict, gap_map)
    logging.info('Distance matrix: ')
    for item in distance_matrix:
        logging.info(item)
    normalization_factor = 1000000 / cg_nucleotide_count
    distance_matrix_output_name = 'distance_matrix_1M.txt'
    print_distance_matrix_phylip(distance_matrix, file_names, args.output, distance_matrix_output_name, normalization_factor)
    print("A distance matrix normalized to a genome size of 1.000.000 has been outputted. The identified core genes spanned {} bases.".format(cg_nucleotide_count), file=sys.stderr)
    run_ccphylo(args.output + '/' + distance_matrix_output_name, args.output + '/tree.newick')


def fast_species_finder(args):
    excluded_list = []
    top_specie = ''
    if args.nanopore != []:
        subset = " ".join(args.nanopore[0:5])
        os.system('kma -i {} -o {} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                  .format(subset, args.output + '/specie_mapping', args.db_dir + '/bac_db/bac_db',
                          args.threads))
        spa_file = args.output + '/specie_mapping.spa'
        top_template = highest_scoring_hit_spa_file(spa_file)
        if top_template != '':
            top_specie = top_template.split(' ')[1] + ' ' + top_template.split(' ')[2]

    if args.illumina != []:
        subset = " ".join(args.illumina[0:10])
        os.system('kma -i {} -o {} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                  .format(subset, args.output + '/specie_mapping', args.db_dir + '/bac_db/bac_db',
                          args.threads))
        spa_file = args.output + '/specie_mapping.spa'
        top_template = highest_scoring_hit_spa_file(spa_file)
        if top_template != '':
            top_specie = top_template.split(' ')[1] + ' ' + top_template.split(' ')[2]

    return excluded_list, top_specie

#TBD CONTINUE HERE
def find_gene_count_outliers(directory):
    # List all .res files
    files = [f for f in os.listdir(directory) if f.endswith('.res')]

    # Initialize a dictionary to store gene counts for each file
    gene_counts = {}

    # Parse each file and count genes
    for file in files:
        with open(os.path.join(directory, file), 'r') as f:
            gene_count = sum(1 for line in f if not line.startswith('#'))
            gene_counts[file] = gene_count

    # Calculate median and threshold
    median_gene_count = sorted(gene_counts.values())[len(gene_counts) // 2]
    threshold = median_gene_count * 0.25

    # Identify outliers and non-outliers
    outliers = [file.split('.')[0] for file, count in gene_counts.items() if count <= threshold]
    non_outliers = [file.split('.')[0] for file, count in gene_counts.items() if count > threshold]

    return outliers, non_outliers


def check_all_species_with_mash(args):
    top_template_count = dict()
    reference_results = dict()
    exclude_list = []

    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]

            os.system('mash sketch -o {}/{}.msh -p 8 -s 25000 {}'.format(args.output, name, file))
            os.system('mash dist {} {}/{}.msh > {}/{}_mash_results.txt'.format(args.db_dir + '/bac_db/25k_bac_db.msh', args.output, name, args.output, name))
            top_template = find_highest_overlap_mash(args.output + '/' + name + '_mash_results.txt')

            with open(args.db_dir + '/bac_db/bac_db.name', 'r') as f:
                for line in f:
                    if line.startswith(top_template):
                        line = line.split()
                        specie = line[1] + ' ' + line[2]
                        if specie in top_template_count:
                            top_template_count[specie] += 1
                        else:
                            top_template_count[specie] = 1
                        print(top_template, name, specie)
                        reference_results[name] = specie
                        if specie == '' or specie == None:
                            exclude_list.append(name)
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            os.system('mash sketch -o {}/{}.msh -p 8 -s 25000 {} {}'.format(args.output, name, args.illumina[i], args.illumina[i+1]))
            os.system('mash dist {} {}/{}.msh > {}/{}_mash_results.txt'.format(args.db_dir + '/bac_db/25k_bac_db.msh', args.output, name, args.output, name))
            top_template = find_highest_overlap_mash(args.output + '/' + name + '_mash_results.txt')
            with open(args.db_dir + '/bac_db/bac_db.name', 'r') as f:
                for line in f:
                    if line.startswith(top_template):
                        line = line.split()
                        specie = line[1] + ' ' + line[2]
                        reference_results[name] = specie
                        print(top_template, name, specie)
                        if specie in top_template_count:
                            top_template_count[specie] += 1
                        else:
                            top_template_count[specie] = 1
                        reference_results[name] = specie
                        if specie == '' or specie == None:
                            exclude_list.append(name)

    for item in top_template_count:
        print (item, top_template_count[item])

    top_specie = max(top_template_count, key=top_template_count.get)
    print('The most common specie is {} with {} hits.'.format(top_specie, top_template_count[top_specie]))



    for file in reference_results:
        print(file, reference_results[file])
        if reference_results[file] != top_specie:
            exclude_list.append(file)

    return exclude_list, top_specie

def find_highest_overlap_mash(mash_results_file):
    highest_overlap_reference = None
    highest_overlap_count = -1
    highest_identity_score = -1.0

    with open(mash_results_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 5:
                continue  # Skip incomplete lines

            reference, _, identity_score, _, overlap = parts
            identity_score = float(identity_score)
            overlap_count = int(overlap.split('/')[0])

            # Update if this line has a higher overlap or equal overlap with a higher identity score
            if overlap_count > highest_overlap_count or (
                    overlap_count == highest_overlap_count and identity_score > highest_identity_score):
                highest_overlap_reference = reference
                highest_overlap_count = overlap_count
                highest_identity_score = identity_score

    return highest_overlap_reference

def run_ccphylo(distance_matrix_file, output_file):
    cmd = 'ccphylo tree --input {} --output {}'.format(distance_matrix_file, output_file)
    os.system(cmd)


def get_genome_size(args, top_specie):
    top_dict = {}
    file_list = os.listdir(args.output)
    for file in file_list:
        if file.startswith('species_mapping_'):
            with open(args.output + '/' + file, 'r') as f:
                for line in f:
                    num = line.split('\t')[1]
                    if top_specie in line:
                        if num in top_dict:
                            top_dict[num] += 1
                        else:
                            top_dict[num] = 1
    top_num = max(top_dict, key=top_dict.get)
    os.system('kma seq2fasta -t_db {} -seqs {} > {}/top_species.fasta'.format(args.db_dir + '/bac_db/bac_db', top_num, args.output))
    genome_size = count_nucleotides(args.output + '/top_species.fasta')
    return genome_size

def count_nucleotides(fasta_file):
    nucleotide_count = 0

    with open(fasta_file, 'r') as file:
        sequence_started = False

        for line in file:
            line = line.strip()

            if line.startswith('>'):
                # Skip sequence identifier lines
                if sequence_started:
                    break
            else:
                # Count nucleotides in the sequence
                nucleotide_count += len(line)
                sequence_started = True

    return nucleotide_count

def check_all_species(args):
    top_template_count = dict()
    reference_results = dict()
    exclude_list = []

    if args.nanopore != []:
        for file in args.nanopore:
            if len(file.split(' ')) == 1:
                name = file.split('/')[-1].split('.')[0]
            else:
                name = file.split(' ')[0].split('/')[-1].split('.')[0]
            os.system('kma -i {} -o {}{} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                      .format(file, args.output + '/species_mapping_', name, args.db_dir + '/bac_db/bac_db',
                              args.threads))
            spa_file = args.output + '/species_mapping_' + name + '.spa'
            top_template = highest_scoring_hit_spa_file(spa_file)
            if top_template != '':
                specie = top_template.split(' ')[1] + ' ' + top_template.split(' ')[2]
                if specie in top_template_count:
                    top_template_count[specie] += 1
                else:
                    top_template_count[specie] = 1
                reference_results[name] = specie
            else:
                reference_results[name] = 'No hits found'
                exclude_list.append(name)
                logging.info('No hits found for nanopore file: {}'.format(file))
            print (name, specie)
    if args.illumina != []:
        for i in range(0, len(args.illumina), 2):
            name = args.illumina[i].split('/')[-1].split('.')[0]
            os.system('kma -i {} {} -o {}{} -t_db {} -mem_mode -t {} -Sparse -ss c' \
                      .format(args.illumina[i], args.illumina[i+1], args.output + '/species_mapping_', name, args.db_dir + '/bac_db/bac_db',
                              args.threads))
            spa_file = args.output + '/species_mapping_' + name + '.spa'
            top_template = highest_scoring_hit_spa_file(spa_file)
            if top_template != '':
                specie = top_template.split(' ')[1] + ' ' + top_template.split(' ')[2]
                if specie in top_template_count:
                    top_template_count[specie] += 1
                else:
                    top_template_count[specie] = 1
                reference_results[name] = specie
            else:
                reference_results[name] = 'No hits found'
                exclude_list.append(name)
                logging.info('No hits found for nanopore file: {}'.format(file))

    for item in top_template_count:
        print(item, top_template_count[item])

    top_specie = max(top_template_count, key=top_template_count.get)
    print('The most common specie is {} with {} hits.'.format(top_specie, top_template_count[top_specie]))


    for file in reference_results:
        print (file, reference_results[file])
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

def load_sequences_from_file(output, gene_list, outliers):
    file_sequences_dict = dict()
    files = os.listdir(output)
    highest_nucleotide_count = 0
    for file in files:
        if file.endswith('.fsa'):
            file_id = file.split('.')[0]
            if file_id not in outliers:
                current_nucleotide_count = 0
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
                            current_nucleotide_count += len(line.strip())
                if current_nucleotide_count > highest_nucleotide_count:
                    highest_nucleotide_count = current_nucleotide_count
    return file_sequences_dict, highest_nucleotide_count



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


                # Dont count gaps test
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

def find_common_genes(directory_path, outliers):
    files = os.listdir(directory_path)

    gene_lists = []

    for file in files:
        if file.endswith('.res'):
            file_id = file.split('.')[0]
            if file_id not in outliers:
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

    if highest_scoring_template == "":
        print(f"Warning: No highest scoring template found in {file_path}")


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