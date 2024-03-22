import os
import sys

path = '/home/people/malhal/cgphylo/test/output_cpo_test'
files = os.listdir(path)

gene_lists = []

for file in files:
    if file.endswith('.res'):
        genes = list()
        with open(path + '/' + file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    genes.append(line[0])
        gene_lists.append(genes)

common = set(gene_lists[0])
for gene_list in gene_lists[1:]:
    common.intersection_update(gene_list)

print (len(common))