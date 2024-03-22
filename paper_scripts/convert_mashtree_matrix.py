def convert_to_phylip_corrected(input_filename, output_filename):
    with open(input_filename, 'r') as infile:
        lines = infile.readlines()

        # Extract headers (sequence identifiers)
        headers = lines[0].strip().split('\t')[1:]  # Skip the first empty cell

        # Initialize an empty list to hold the matrix rows
        matrix = []

        # Process each line to build the full matrix
        for line in lines[1:]:
            parts = line.strip().split('\t')
            row = [float(dist) for dist in parts[1:]]  # Convert distances to float
            matrix.append(row)

        # Prepare the PHYLIP output with the corrected format
        with open(output_filename, 'w') as outfile:
            # Write the number of sequences
            outfile.write(f"{len(headers)}\n")

            # Write each row in PHYLIP format, placing zeros in the upper half
            for i, header in enumerate(headers):
                distances = [f"{matrix[i][j]:.6f}" if j < i else "0.000000" for j in range(len(headers))]
                distances_str = '\t'.join(distances)  # Create the tab-separated string of distances
                outfile.write(f"{header}\t{distances_str}\n")  # Use variable outside of curly braces


# Specify your input and output filenames
input_filename = 'matrix'
output_filename = 'phylip_matrix'

# Convert the matrix to the corrected PHYLIP format
convert_to_phylip_corrected(input_filename, output_filename)

print(f"Matrix has been correctly formatted and saved to {output_filename}")