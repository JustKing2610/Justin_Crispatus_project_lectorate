from Bio import SeqIO
import pandas as pd
import re

# Path to the input FASTA file
fasta_file_path = "/path/to/fasta/all_slph.fasta"

# Path to the Excel file with identifiers and sequence types
excel_file_path = "path/to/excel/identifiers_STs.xlsx"

# Path to the output FASTA file
output_fasta_path = "/Path/to/new/fasta/with/STheaders/all_slph_st.fasta "

# Read Excel file into a pandas DataFrame
df = pd.read_excel(excel_file_path, names=['Identifier', 'Sequence_types'])

# Create a dictionary to map identifiers to sequence types
identifier_map = dict(zip(df['Identifier'], df['Sequence_types']))

# Open the output FASTA file for writing
with open(output_fasta_path, 'w') as output_fasta:
    # Iterate over sequences in the input FASTA file
    for record in SeqIO.parse(fasta_file_path, 'fasta'):
        # Use regular expression to extract everything starting from 'GCA_' until the first period
        match = re.match(r"^.*?(GCA_\d+)\..*$", record.id)
        
        if match:
            base_identifier = match.group(1)
            
            # Check if the base identifier is in the mapping dictionary
            if base_identifier in identifier_map:
                # Concatenate the sequence type to the modified identifier
                new_identifier = f'>{identifier_map[base_identifier]}_{base_identifier}'
                
                # Create a new SeqRecord with the updated identifier and the original sequence
                updated_record = record
                updated_record.id = new_identifier
                updated_record.description = ''
                
                # Write the updated SeqRecord to the output FASTA file
                SeqIO.write(updated_record, output_fasta, 'fasta')
                
                # Print for debugging
                print(f"Updated: {new_identifier}")
        else:
            print(f"Skipping record with invalid identifier format: {record.id}")

print(f"Output written to {output_fasta_path}")

