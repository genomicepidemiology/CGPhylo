import pandas as pd
import sys

def count_unique_center_names(file_path):
    try:
        # Load the CSV file
        df = pd.read_csv(file_path)

        # Randomly sample 100 lines
        sample_df = df.sample(n=100, random_state=1)

        # Count unique CenterName values
        unique_center_names_count = sample_df['CenterName'].nunique()

        # Print the number of unique CenterName values to stderr
        print(f"Number of unique CenterName values in the sample: {unique_center_names_count}", file=sys.stderr)

        # Return the sampled DataFrame
        return sample_df
    except Exception as e:
        return f"An error occurred: {e}"

def write_run_ids_to_file(df, filename):
    try:
        # Extract Run numbers and write to a file
        run_ids = df['Run']
        run_ids.to_csv(filename, index=False, header=False)
        print(f"Run IDs written to {filename}")
    except Exception as e:
        print(f"An error occurred while writing Run IDs: {e}")

# File path
file_path = '/Users/malhal/Downloads/SraRunInfo_ont.csv'

# Execute the function to get the sampled data
sampled_df = count_unique_center_names(file_path)

# Print the CSV of the sampled DataFrame to stdout
if isinstance(sampled_df, pd.DataFrame):
    print(sampled_df.to_csv(index=False))
    # Write Run IDs to run_ids.txt
    write_run_ids_to_file(sampled_df, 'run_ids.txt')
else:
    print(sampled_df)
