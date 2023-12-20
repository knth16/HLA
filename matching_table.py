###############################################################################
# Script: matching_table.py
# Author: Kenneth Valerio Aguilar
# Date: 12/20/2023
# Version: 1.0
# Purpose: Create table with matching information between the genotyping and pedigree
# Input Requirements: serotype_table.tsv and summary_table.csv
# Usage: run the script like: python3 serotype_table.py
# Output: summary serotype_table.tsv 
################################################################################

#Import required packages
import pandas as pd
import os
import re

#read in data
serotypes = pd.read_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/serotype_table.tsv' sep='\t')
data = pd.read_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/summary_table.csv',sep=',')

#create empty lists to fill with for loop
serotypes_dict = {}

for index1, row1 in data.iterrows():
    for index2, row2 in serotypes.iterrows():
 # Check for DQA1_Gtype1 and DQB1_Gtype1
        if (row2['DQA1'] in row1['DQA1_Gtype1'][9:14] or row2['DQA1'] in row1['DQA1_Gtype2'][9:14]) and \
           (row2['DQB1'] in row1['DQB1_Gtype1'][9:14] or row2['DQB1'] in row1['DQB1_Gtype2'][9:14]):
            if row1['Basename'] not in serotypes_dict.keys():
                serotypes_dict[row1['Basename']] = [row2['Subtype']]
            else:
                serotypes_dict[row1['Basename']].append(row2['Subtype'])
        # Check for DQA1_Gtype1 and DQB1_Gtype2
        elif (row2['DQA1'] in row1['DQA1_Gtype1'][9:14] or row2['DQA1'] in row1['DQA1_Gtype2'][9:14]) and \
             (row2['DQB1'] in row1['DQB1_Gtype2'][9:14]):
            if row1['Basename'] not in serotypes_dict.keys():
                serotypes_dict[row1['Basename']] = [row2['Subtype']]
            else:
                serotypes_dict[row1['Basename']].append(row2['Subtype'])
        # Check for DQA1_Gtype2 and DQB1_Gtype1
        elif (row2['DQA1'] in row1['DQA1_Gtype2'][9:14]) and \
             (row2['DQB1'] in row1['DQB1_Gtype1'][9:14] or row2['DQB1'] in row1['DQB1_Gtype2'][9:14]):
            if row1['Basename'] not in serotypes_dict.keys():
                serotypes_dict[row1['Basename']] = [row2['Subtype']]
            else:
                serotypes_dict[row1['Basename']].append(row2['Subtype'])
        # Check for DQA1_Gtype2 and DQB1_Gtype2
        elif (row2['DQA1'] in row1['DQA1_Gtype2'][9:14]) and \
             (row2['DQB1'] in row1['DQB1_Gtype2'][9:14]):
            if row1['Basename'] not in serotypes_dict.keys():
                serotypes_dict[row1['Basename']] = [row2['Subtype']]
            else:
                serotypes_dict[row1['Basename']].append(row2['Subtype'])

for basename, subtypes in serotypes_dict.items():
    while len(subtypes) < 4:
        subtypes.append("NA")



result_df = pd.DataFrame.from_dict(serotypes_dict, orient='index')
# Rename columns
result_df.columns = ['Serotype_1', 'Serotype_2', 'Serotype_3', 'Serotype_4']

# Reset index to move 'Basename' from index to a regular column
result_df.reset_index(inplace=True)
result_df.rename(columns={'index': 'Basename'}, inplace=True)

#lets add column with Korponay coding
Program_obtained_coding = [] 
Korponay_coding = {'1':2.5,'2':2.2,'3':[7.2,7.4,7.5,7.6],'8':8.1}

#create function to classify
def apply_korponay_coding(row):
    coding_result = []

    for key, value in Korponay_coding.items():
        if isinstance(value, list):
            if any(subtype in value for subtype in row[['Serotype_1', 'Serotype_2', 'Serotype_3', 'Serotype_4']]):
                coding_result.append(key)
        else:
            if value in row[['Serotype_1', 'Serotype_2', 'Serotype_3', 'Serotype_4']].values:
                coding_result.append(key)

    return ''.join(coding_result) if coding_result else 'NA'

# Apply the function to create a new column
result_df['program_Coding'] = result_df.apply(apply_korponay_coding, axis=1)
result_df['program_Coding'] = result_df['program_Coding'].apply(lambda x: x + 'X' if len(x) == 1 else x)


# export result
result_df.to_csv('/molbio/projects/hla_genotyping/hlahd/CD_analysis/matching_table.csv', sep='\t', index=False)


