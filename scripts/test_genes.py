import os
import sys

path = '/home/people/malhal/mintyper2/data/030'
files = os.listdir(path)

for file in files:
    name = file.split('.')[0]
    cmd = 'kma -i {} -o {} -t_db /home/people/malhal/mintyper2/new_alleles -ID 50 -ont -md 1.5 -mct 0.5 -t 8'.format(path + '/' + file, name)
    os.system(cmd)