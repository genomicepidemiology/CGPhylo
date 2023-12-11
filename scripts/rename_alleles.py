import os
import sys

path = '/Users/malhal/Downloads/Escherichia_coli_cgMLST_alleles'
files = os.listdir(path)

for file in files:
    file_name = file.split('.')[0]
    with open (path + '/renamed_' + file, 'w') as outfile:
        with open(path + '/' + file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    number = line.strip()[1:]
                    print ('>{}_{}'.format(file_name, number), file=outfile)
                else:
                    print (line.strip(), file=outfile)
