def process_matrix_file(input_file_path, output_file_path):
    with open(input_file_path, 'r') as input_file, open(output_file_path, 'w') as output_file:
        for line in input_file:
            parts = line.strip().split('\t')
            if parts[0]:  # Check if there's a sample name to process
                # Process the sample name
                parts[0] = parts[0].split('.')[0]
                sample_name_parts = parts[0].split('_')[:2]
                new_sample_name = '_'.join(sample_name_parts)
                parts[0] = new_sample_name
            # Write the processed line to the output file
            output_file.write('\t'.join(parts) + '\n')

# Define the input and output file paths
#input_files = ['distance_matrix_1M.txt']
#input_files = ['distmatrix.txt']
input_files = ['corrected_matrix']
#output_files = ['corrected_distance_matrix_1M.txt']
output_files = ['corrected_matrix_final']
#input_files = ['phylip_matrix']
#output_files = ['corrected_matrix']


# Process each file
for input_file, output_file in zip(input_files, output_files):
    process_matrix_file(input_file, output_file)
