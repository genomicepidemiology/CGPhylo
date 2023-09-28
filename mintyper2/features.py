import os
import sys
import numpy as np


def produce_features(args):
    os.mkdir('output')
    illumina_list = []
    all_files = args.nanopore
    for i in range(0, len(args.illumina), 2):
        string = args.illumina[i] + ' ' + args.illumina[i+1]
        illumina_list.append(string)
        all_files.append(string)
    produce_kmers(args)
    matrix_5mers = np.zeros((len(all_files), len(all_files)))
    for i in range(len(matrix_5mers)):
        for j in range(len(matrix_5mers)):
            if i == j:
                matrix_5mers[i][j] = 1
            else:
                file_1 = all_files[i]
                file_2 = all_files[j]
                if len(file_1.split(' ')) == 1:
                    name_1 = file_1.split('/')[-1].split('.')[0]
                else:
                    name_1 = file_1.split(' ')[0].split('/')[-1].split('.')[0]
                if len(file_2.split(' ')) == 1:
                    name_2 = file_2.split('/')[-1].split('.')[0]
                else:
                    name_2 = file_2.split(' ')[0].split('/')[-1].split('.')[0]
                read_set_1 = read_file_to_set('output/{}_5mers.txt'.format(name_1))
                read_set_2 = read_file_to_set('output/{}_5mers.txt'.format(name_2))
                matrix_5mers[i][j] = jaccard_index(read_set_1, read_set_2)

    print ('5mers')
    print (matrix_5mers)

    matrix_13mers = np.zeros((len(all_files), len(all_files)))
    for i in range(len(matrix_13mers)):
        for j in range(len(matrix_13mers)):
            if i == j:
                matrix_13mers[i][j] = 1
            else:
                file_1 = all_files[i]
                file_2 = all_files[j]
                if len(file_1.split(' ')) == 1:
                    name_1 = file_1.split('/')[-1].split('.')[0]
                else:
                    name_1 = file_1.split(' ')[0].split('/')[-1].split('.')[0]
                if len(file_2.split(' ')) == 1:
                    name_2 = file_2.split('/')[-1].split('.')[0]
                else:
                    name_2 = file_2.split(' ')[0].split('/')[-1].split('.')[0]
                read_set_1 = read_file_to_set('output/{}_13mers.txt'.format(name_1))
                read_set_2 = read_file_to_set('output/{}_13mers.txt'.format(name_2))
                matrix_13mers[i][j] = jaccard_index(read_set_1, read_set_2)
    print ('13mers')
    print (matrix_13mers)


    matrix_21mers = np.zeros((len(all_files), len(all_files)))
    for i in range(len(matrix_21mers)):
        for j in range(len(matrix_21mers)):
            if i == j:
                matrix_21mers[i][j] = 1
            else:
                file_1 = all_files[i]
                file_2 = all_files[j]
                if len(file_1.split(' ')) == 1:
                    name_1 = file_1.split('/')[-1].split('.')[0]
                else:
                    name_1 = file_1.split(' ')[0].split('/')[-1].split('.')[0]
                if len(file_2.split(' ')) == 1:
                    name_2 = file_2.split('/')[-1].split('.')[0]
                else:
                    name_2 = file_2.split(' ')[0].split('/')[-1].split('.')[0]
                read_set_1 = read_file_to_set('output/{}_21mers.txt'.format(name_1))
                read_set_2 = read_file_to_set('output/{}_21mers.txt'.format(name_2))
                matrix_21mers[i][j] = jaccard_index(read_set_1, read_set_2)
    print ('21mers')
    print (matrix_21mers)

    return 'test'
def produce_kmers(args):
    """Produces kmer files for the input file."""
    for item in args.illumina:
        name = item.split('/')[-1].split('.')[0]
        os.system('kmc -k1 -cs1000000000000 -b {} output/{}_1mers . > output/{}_1mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_1mers output/{}_1mers.txt'.format(name, name))
        os.system('kmc -k5 -cs1000000000 {} output/{}_5mers . > output/{}_5mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_5mers output/{}_5mers.txt'.format(name, name))
        os.system('kmc -k13 -cs1000000000 {} output/{}_13mers . > output/{}_13mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_13mers output/{}_13mers.txt'.format(name, name))
        os.system('kmc -k21 -cs1000000000 {} output/{}_21mers . > output/{}_21mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_21mers output/{}_21mers.txt'.format(name, name))

    for item in args.nanopore:
        name = item.split('/')[-1].split('.')[0]
        os.system('kmc -k1 -cs1000000000000 -b {} output/{}_1mers output/ > output/{}_1mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_1mers output/{}_1mers.txt'.format(name, name))
        os.system('kmc -k5 -cs1000000000 {} output/{}_5mers output/ > output/{}_5mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_5mers output/{}_5mers.txt'.format(name, name))
        os.system('kmc -k13 -cs1000000000 {} output/{}_13mers output/ > output/{}_13mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_13mers output/{}_13mers.txt'.format(name, name))
        os.system('kmc -k21 -cs1000000000 {} output/{}_21mers output/ > output/{}_21mers_stats.txt'.format(item, name, name))
        os.system('kmc_dump output/{}_21mers output/{}_21mers.txt'.format(name, name))
def read_file_to_set(filename):
    """Reads kmer strings from a file and returns a set of those kmers."""
    kmers = set()
    if len(filename.split(' ')) == 1:
        with open(filename, 'r') as f:
            for line in f:
                kmer, _ = line.strip().split('\t')
                kmers.add(kmer)
    else:
        with open(filename.split(' ')[0], 'r') as f:
            for line in f:
                kmer, _ = line.strip().split('\t')
                kmers.add(kmer)
        with open(filename.split(' ')[1], 'r') as f:
            for line in f:
                kmer, _ = line.strip().split('\t')
                kmers.add(kmer)
    return kmers

def jaccard_index(set1, set2):
    """Calculates the Jaccard index between two sets."""
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

def main():
    file1 = input("Enter the path to the first file: ")
    file2 = input("Enter the path to the second file: ")

    kmers1 = read_file_to_set(file1)
    kmers2 = read_file_to_set(file2)

    jaccard = jaccard_index(kmers1, kmers2)

    print(f"Jaccard Index between the two files: {jaccard:.4f}")
