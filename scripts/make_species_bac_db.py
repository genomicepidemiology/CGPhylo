import os
import sys

species_string_list = list()

with open('cgmlst_list.txt', 'r') as f:
    for line in f:
        if line.startswith('Brucella spp.'):
            species_string_list.append('Brucella sp.') #Naming difference between RefSeq and PubMLST
        elif line.startswith('Burkholderia mallei (FLI)'):
            pass
        elif line.startswith('Campylobacter'):
            species_string_list.append('Campylobacter jejuni')
            species_string_list.append('Campylobacter coli')
        elif line.startswith('Cronobacter'):
            species_string_list.append('Cronobacter sakazakii')
            species_string_list.append('Cronobacter malonaticus')
        elif line.startswith('Klebsiella'):
            species_string_list.append('Klebsiella pneumoniae')
            species_string_list.append('Klebsiella variicola')
            species_string_list.append('Klebsiella quasipneumoniae')
        elif line.startswith('Mycobacterium'):
            species_string_list.append('Mycobacterium tuberculosis')
            species_string_list.append('Mycobacterium bovis')
            species_string_list.append('Mycobacterium africanum')
            species_string_list.append('Mycobacterium canettii')
        else:
            species_string_list.append(line.strip().split(' ')[0] + ' ' + line.strip().split(' ')[1])

db_path = '/home/people/malhal/contamErase_db/bac_db.name'

string_list = list()

t = 1
with open(db_path, 'r') as f:
    for line in f:
        line = line.split(' ')
        specie = line[1] + ' ' + line[2]
        if specie in species_string_list:
            string_list.append(t)
            t += 1

seq_string = ",".join(str(x) for x in string_list)

os.system('kma seq2fasta -t_db /home/people/malhal/contamErase_db/bac_db -seqs {} > cgmlst_species.fasta'.format(seq_string))