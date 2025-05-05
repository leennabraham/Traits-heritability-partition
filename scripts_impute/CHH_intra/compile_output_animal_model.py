import os
import pandas as pd

# Define the directory containing the CSV files
directory = "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/intra_impute/intra_CHH/output_files"
output_file = "/projects/ag-demeaux/labraha3/EpiDom_A_lyrata/Animal_model/intra_impute/intra_CHH/output_files/Intra_compiled_animal_model_CHH_impute.txt"

# Get a list of all CSV files in the directory
data_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith('.csv')]

# Initialize an empty list to hold the data from each file
data_list = []

# Loop through each file, read it, and append the data to the list
for file in data_files:
    # Read the file with the correct delimiter and quote handling
    df = pd.read_csv(file, sep="\t", quotechar='"', engine='python')
    
    # Append the dataframe to the list
    data_list.append(df)

# Ensure all dataframes have the same columns by reindexing them with the first file's columns
first_columns = data_list[0].columns  # Assuming all files should have the same columns
data_list = [df.reindex(columns=first_columns) for df in data_list]

# Concatenate all the dataframes in the list into a single dataframe
combined_data = pd.concat(data_list, ignore_index=True)

# Save the combined data to a text file
combined_data.to_csv(output_file, sep="\t", index=False)

print(f"Compiled data has been saved to {output_file}")