#!/usr/bin/env python3

import argparse
import sys
import os
import logging
import subprocess

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

#NOTE Make sure to remove Burkholderia mallei (RKI) and only use Burkholderia mallei (FLI). Make sure to rename the folder to Burkholderia_mallei_cgMLST_alleles.

def main(args):
    """Main function"""
    check_kma_installed()
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    decompress_tar_archive(args.zip, args.output)
    allele_folders = os.listdir(args.output)
    for allele_folder in allele_folders:
        if not allele_folder.startswith('.'):
            path = args.output + '/' + allele_folder
            rename_alleles(path)
            concat_renamed_fasta_files(path)
            os.system('kma index -i {}/{}_complete.fasta -o {}/{}_complete'.format(path, allele_folder, path, allele_folder))
            #os.system('rm {}/*.fasta'.format(path))
            # Build consensus alleles
            gene_names = derive_gene_names('{}/{}_complete.name'.format(path, allele_folder))
            gene_lengths = derive_gene_lengths(gene_names, '{}/{}_complete.fasta'.format(path, allele_folder))
            consensus_genes = derive_consensus_genes(gene_lengths, '{}/{}_complete.fasta'.format(path, allele_folder), gene_names)

            # Writing to file
            with open('{}/{}_consensus_genes.fasta'.format(path, allele_folder), 'w') as f:
                for gene in consensus_genes:
                    f.write('>' + gene + '\n')
                    f.write(consensus_genes[gene] + '\n')

            os.system('kma index -i {}/{}_consensus_genes.fasta -o {}/{}_consensus_genes'.format(path, allele_folder, path, allele_folder))


            # Produce gapmap


def derive_gene_names(name_file):
    gene_names = set()
    with open(name_file, 'r') as f:
        for line in f:
            allele_name = line.strip()
            gene_name = allele_name.split('_')[:-1]
            gene_name = '_'.join(gene_name)
            gene_names.add(gene_name)
    return gene_names

def derive_gene_lengths(gene_names, fsa_file):
    gene_dict = {name: {} for name in gene_names}
    sequence = ''
    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    length = len(sequence)
                    gene_dict[gene_name][length] = gene_dict[gene_name].get(length, 0) + 1
                allele_name = line.strip()[1:]
                gene_name = allele_name.split('_')[:-1]
                gene_name = '_'.join(gene_name)
                sequence = ''
            else:
                sequence += line.strip()
    if sequence != '':
        length = len(sequence)
        gene_dict[gene_name][length] = gene_dict[gene_name].get(length, 0) + 1

    return gene_dict

def derive_consensus_genes(gene_lengths, fsa_file):
    gene_dict = {}
    name_set = set()

    # Filter out alleles with occurrences less than 10
    for name in gene_lengths:
        for length, count in gene_lengths[name].items():
            if count >= 10:
                name_set.add(name + '_len_' + str(length))

    for name in name_set:
        gene_dict[name] = [[0, 0, 0, 0] for _ in range(int(name.split('_')[-1]))]

    t = 0
    sequence = ''
    with open(fsa_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if sequence != '':
                    t += 1
                    if t % 10000 == 0:
                        print(t)
                    search_string = gene_name + '_len_' + str(len(sequence))
                    if search_string in gene_dict:
                        for i in range(len(sequence)):
                            if sequence[i].upper() == 'A':
                                gene_dict[search_string][i][0] += 1
                            elif sequence[i].upper() == 'C':
                                gene_dict[search_string][i][1] += 1
                            elif sequence[i].upper() == 'G':
                                gene_dict[search_string][i][2] += 1
                            elif sequence[i].upper() == 'T':
                                gene_dict[search_string][i][3] += 1
                allele_name = line.strip()[1:]
                gene_name = allele_name.split('_')[:-1]
                gene_name = '_'.join(gene_name)
                sequence = ''
            else:
                sequence += line.strip()

    if sequence != '':
        t += 1
        if t % 10000 == 0:
            print(t)
        search_string = gene_name + '_len_' + str(len(sequence))
        if search_string in gene_dict:
            for i in range(len(sequence)):
                if sequence[i].upper() == 'A':
                    gene_dict[search_string][i][0] += 1
                elif sequence[i].upper() == 'C':
                    gene_dict[search_string][i][1] += 1
                elif sequence[i].upper() == 'G':
                    gene_dict[search_string][i][2] += 1
                elif sequence[i].upper() == 'T':
                    gene_dict[search_string][i][3] += 1

    consensus_genes = {}
    for gene in gene_dict:
        consensus_genes[gene] = ''
        for i in range(len(gene_dict[gene])):
            max_position = gene_dict[gene][i].index(max(gene_dict[gene][i]))
            if max_position == 0:
                consensus_genes[gene] += 'A'
            elif max_position == 1:
                consensus_genes[gene] += 'C'
            elif max_position == 2:
                consensus_genes[gene] += 'G'
            elif max_position == 3:
                consensus_genes[gene] += 'T'

    return consensus_genes

def check_kma_installed():
    try:
        # Attempt to execute the KMA command and get its version
        result = subprocess.run(["kma", "-v"], capture_output=True, text=True, check=True)
        if result.returncode == 0:
            return "KMA is installed and available."
        else:
            return "KMA is installed but there was an error executing it."
    except FileNotFoundError:
        return "KMA is not installed or not available in PATH."


def clean(path):
    # Remove hidden files
    files = os.listdir(path)
    for file in files:
        sub_files = os.listdir(path + '/' + file)
        for sub_file in sub_files:
            if sub_file.startswith('._'):
                os.remove(path + '/' + file + '/' + sub_file)
            if sub_file.startswith('rename'):
                os.remove(path + '/' + file + '/' + sub_file)
            if 'complete' in sub_file:
                os.remove(path + '/' + file + '/' + sub_file)


def concat_renamed_fasta_files(path):
    species_name = path.split('/')[-1]
    os.system('cat {}/renamed_* > {}/{}_complete.fasta'.format(path, path, species_name))

def rename_alleles(path):
    files = os.listdir(path)

    for file in files:
        file_name = file.split('.')[0]
        if not file.startswith('.'):
            with open(path + '/renamed_' + file, 'w') as outfile:
                with open(path + '/' + file, 'r') as f:
                    print ('renamed: ' + path + '/' + file)
                    for line in f:
                        if line.startswith('>'):
                            number = line.strip()[1:]
                            print('>{}_{}'.format(file_name, number), file=outfile)
                        else:
                            print(line.strip(), file=outfile)


def decompress_tar_archive(path, output):
    """Decompress a tar archive"""
    os.system('tar -xvf {} -C {}'.format(path, output))


if __name__ == '__main__':
    # initialize the options parser
    parser = argparse.ArgumentParser('build_cg_db', add_help=False)

    parser.add_argument('--zip', action="store", type=str, dest='zip',
                        help='Zip archive of cgmlst db which have been downloaded.')
    parser.add_argument('--output', action="store", type=str, default='output', dest="output",
                        help="Output directory")
    parser.add_argument('-h', '--help', action='help', help='Show this help message and exit')

    args = parser.parse_args()

    main(args)