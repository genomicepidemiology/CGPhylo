import os
import sys

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

def derive_consensus_genes(gene_lengths, fsa_file, gene_names):
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


# Usage
gene_names = derive_gene_names('new_alleles.name')
gene_lengths = derive_gene_lengths(gene_names, 'new_alleles.fasta')
consensus_genes = derive_consensus_genes(gene_lengths, 'new_alleles.fasta', gene_names)

# Writing to file
with open('consensus_genes_2.fasta', 'w') as f:
    for gene in consensus_genes:
        f.write('>' + gene + '\n')
        f.write(consensus_genes[gene] + '\n')

