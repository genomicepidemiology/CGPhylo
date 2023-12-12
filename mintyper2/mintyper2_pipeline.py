import os
import sys

def mintyper2_pipeline():
    #KMA ALIGnment
    gene_list = find_common_genes('/home/people/malhal/mintyper2/test/output_cpo_test')
    #Right now we ONLY use perfect length matches
    print (gene_list)

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
            genes = []
            with open(os.path.join(directory_path, file), 'r') as f:
                for line in f:
                    if not line.startswith('#'):
                        line = line.strip().split('\t')
                        genes.append(line[0])
            gene_lists.append(genes)

    # Find common genes
    if gene_lists:
        common = set(gene_lists[0])
        for gene_list in gene_lists[1:]:
            common.intersection_update(gene_list)
        return common
    else:
        return set()  # Return an empty set if no gene lists were found
