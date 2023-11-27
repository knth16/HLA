################################################################################
# Script: table_extraction.py
# Author: Kenneth Valerio Aguilar
# Date: 11/27/2023
# Version: 1.0
# Purpose: Create analysis table from HLA genotyping data
# Input Requirements: Requires results.txt files with genotyping data
# Usage: run the script
# Output: summary .csv table
################################################################################

# import the required libraries
import pandas as pd
import os
import re

# set root directory
root_dir = "/molbio/projects/hla_genotyping/hlahd/CD_analysis"

# initialize an empty DataFrame
Final_table = pd.DataFrame(columns=['Region', 'Gtype1', 'Gtype2', 'Basename'])

# iterate over subdirectories in CD_analysis
for subdir in os.listdir(root_dir):
    subdir_path = os.path.join(root_dir, subdir)

    # check if it's a directory
    if os.path.isdir(subdir_path):
        for file_name in os.listdir(os.path.join(subdir_path, 'Result')):
            if file_name.endswith("_final.result.txt"):
                file_path = os.path.join(subdir_path, 'Result', file_name)

                # read file content
                with open(file_path, 'r') as file:
                    file_content = file.read()
                    updated_content1 = file_content.replace('Not typed', 'Not_typed')
                    updated_content = re.sub(r'\s{2,}', ' ', updated_content1)

                # write back updated content
                with open(file_path, 'w') as updated_file:
                    updated_file.write(updated_content)

                #extract information and append to the Final_table dataframe
                data = pd.read_csv(file_name,sep=' ',header=None,names=['Region','Gtype1','Gtype2'])
                data2=data.iloc[4:6]
                data2['Basename']= file_name.split('_')[0]

                #concatenate and update final table
                Final_table = pd.concat([Final_table, data2], ignore_index=True)
            else:
                continue
    else:
        continue
#Export the final table
Final_table.to_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/summary_table.csv', index=False)