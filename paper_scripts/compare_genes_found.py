import os
import sys

def find_common_elements(matrix):
    if not matrix:
        return []

    # Convert the first sublist to a set
    common_elements = set(matrix[0])

    # Iterate over the rest of the sublists and update the set
    for sublist in matrix[1:]:
        common_elements = common_elements.intersection(sublist)

    # Convert the set back to a list before returning
    return list(common_elements)

path = '/home/people/malhal/cgphylo/test/output_cpo_test'
files = os.listdir(path)

data = []
index = 0
for file in files:
    if file.endswith('.res'):
        data.append([])
        with open(file, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    data[index].append(line[0].split('_')[0].strip())
        index += 1

common = find_common_elements(data)
print (len(common))


