import os
import sys
import numpy as np


def produce_features(args):
    os.mkdir('output')
    #produce_kmers(args)
    illumina_list = []
    for i in range(0, len(args.illumina), 2):
        string = args.illumina[i] + ' ' + args.illumina[i+1]
        illumina_list.append(string)

    matrix = np.zeros((len(args.nanopore), len(illumina_list)))
    print (matrix)

def produce_kmers(args):
    """Produces kmer files for the input file."""
    for item in args.illumina:
        name = item.split('/')[-1].split('.')[0]
        os.system('kmc -k1 -cs1000000000000 -b {} output/1mers . > output/{}_1mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/1mers output/1mers.txt')
        os.system('kmc -k5 -cs1000000000 {} output/5mers . > output/{}_5mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/5mers output/5mers.txt')
        os.system('kmc -k13 -cs1000000000 {} output/13mers . > output/{}_13mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/13mers output/13mers.txt')
        os.system('kmc -k21 -cs1000000000 {} output/21mers . > output/{}_21mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/21mers output/21mers.txt')
    for item in args.nanopore:
        name = item.split('/')[-1].split('.')[0]
        os.system('kmc -k1 -cs1000000000000 -b {} output/1mers . > output/{}_1mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/1mers output/1mers.txt')
        os.system('kmc -k5 -cs1000000000 {} output/5mers . > output/{}_5mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/5mers output/5mers.txt')
        os.system('kmc -k13 -cs1000000000 {} output/13mers . > output/{}_13mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/13mers output/13mers.txt')
        os.system('kmc -k21 -cs1000000000 {} output/21mers . > output/{}_21mers_stats.txt'.format(item, name))
        os.system('kmc_dump output/21mers output/21mers.txt')
def read_file_to_set(filename):
    """Reads kmer strings from a file and returns a set of those kmers."""
    kmers = set()
    with open(filename, 'r') as f:
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
