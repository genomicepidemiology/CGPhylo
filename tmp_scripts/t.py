import pandas as pd

def filter_sequences(file_path):
    try:
        # Load the CSV file
        df = pd.read_csv(file_path)

        # Filter for TaxID = 1280 (as integer)
        df = df[df['TaxID'] == 1280]

        # Convert bases to numeric and filter between 100M and 2bn
        df['bases'] = pd.to_numeric(df['bases'], errors='coerce')
        df = df[df['bases'].between(100000000, 2000000000)]

        # Convert LoadDate to datetime and filter for 2019 or later
        df['LoadDate'] = pd.to_datetime(df['LoadDate'], errors='coerce')
        df = df[df['LoadDate'].dt.year >= 2019]

        return df
    except Exception as e:
        return str(e)

# File path
file_path = '/Users/malhal/Downloads/SraRunInfo_ont.csv'

try:
    # Apply the filters and get the result
    filtered_df = filter_sequences(file_path)

    # Print the filtered dataframe
    print(filtered_df.to_csv(index=False))
except Exception as e:
    print("An error occurred:", e)
