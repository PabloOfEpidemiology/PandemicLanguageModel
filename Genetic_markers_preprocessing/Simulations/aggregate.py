import pandas as pd
import glob
import os
import re

def determine_type(filename):
    """
    Extract the Type from the filename by taking the substring before the numbers start.
    """
    # Remove the trailing '_simulations.csv' from the filename
    base_filename = filename.replace('_simulations.csv', '')

    # Use regex to extract the Type
    match = re.match(r'^([^\d]+?)_(?=\d)', base_filename)
    if match:
        return match.group(1)
    else:
        # If no match, return the base filename
        return base_filename

def main():
    # Define the subdirectories to search
    subdirectories = ['FG', 'Spike']

    # List to hold all matching CSV files
    csv_files = []

    # Iterate over each subdirectory and collect CSV files
    for subdir in subdirectories:
        # Construct the search pattern for each subdirectory
        pattern = os.path.join(subdir, "*_simulations.csv")
        found_files = glob.glob(pattern)
        csv_files.extend(found_files)

    if not csv_files:
        print("No simulation CSV files found in 'FG/' or 'Spike/' directories.")
        return

    # List to hold individual DataFrames
    df_list = []

    for file in csv_files:
        try:
            # Read the CSV file
            df = pd.read_csv(file)

            # Determine the Type based on the filename
            type_value = determine_type(os.path.basename(file))

            # Add the Type column
            df['Type'] = type_value

            # Optionally, you can add a column for the filename or other metadata
            # df['Source_File'] = os.path.basename(file)

            df_list.append(df)
            print(f"Processed file: {file} as Type: {type_value}")
        except Exception as e:
            print(f"Error processing file {file}: {e}")

    if df_list:
        # Concatenate all DataFrames
        aggregate_df = pd.concat(df_list, ignore_index=True)

        # Define the output filename (saved in the WD)
        output_filename = "aggregate_28thNov.csv"

        # Save the aggregated DataFrame to CSV
        aggregate_df.to_csv(output_filename, index=False)
        print(f"\nAggregation complete. Saved to {output_filename}")
    else:
        print("No dataframes to concatenate.")

if __name__ == "__main__":
    main()
