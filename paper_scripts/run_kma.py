import os
import sys

path = '/home/people/malhal/cgphylo/data/030'
files = os.listdir(path)

output = 'output_cpo_test'
os.system('mkdir {}'.format(output))

for file in files:
    name = file.split('.')[0]
    cmd = 'kma -i {} -o {}/{} -t_db /home/people/malhal/cgphylo/consensus_genes_db -t 8 -ID 30 -ont -md 5 -eq 14 -mct 0.5'.format(path + '/' + file, output, name)
    os.system(cmd)

