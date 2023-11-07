# This program works with the prokka tsv results file containing the COG names
# The COG files are from https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/
# The COG files were modified to include column headers
# COG character IDs are listed in importance if there are multiple IDs (letters)

# Example
# python prokka-cog.py 3300042090_263_prokka.tsv

import os, sys
import pandas as pd
import numpy as np

## Add a check to see if there is the right number of arguments 
if len(sys.argv) != 2 :
    print("This program takes as input a prokka result tsv file")
    sys.exit(2)
    
prokka_file = sys.argv[1];

# read in files

prokka = pd.read_csv(prokka_file,sep = '\t')
cog_def = pd.read_csv('cog-20.def.tsv',sep = '\t')
cog_fun = pd.read_csv('fun-20.tsv',sep = '\t')

# merge prokka and cog_def dataframes
prokka_cog = pd.merge(prokka, cog_def, on = 'COG', how = 'left')
 
########### Get only the first character of the cog charater ID
# Copy dataframe
prokka_cog1 = prokka_cog.copy(deep=True) 
# iterate over rows to find get just the cog charater ID
for ind in prokka_cog1.index:
    cog_characters = str(prokka_cog1['COG_category_ID'][ind])
    if cog_characters !=  'nan':
        first_letter = cog_characters[0]
        prokka_cog1.loc[ind, "COG_category_ID"] = first_letter

# merge with cog descriptions
prokka_cog1_fun = pd.merge(prokka_cog1, cog_fun, on = 'COG_category_ID', how = 'left')

# write file 
prokka_cog1_fun.to_csv('prokka-cog_first_letter.tsv', sep="\t")

############## Make new rows with each cog character ID/letter
# Copy the dataframe
prokka_cog2 = prokka_cog.copy(deep=True)
# Add new rows for each character ID of a multiple character ID call
for ind in prokka_cog2.index:
    row_list = prokka_cog2.loc[ind, :].values.flatten().tolist()
    cog_characters = str(row_list[7])
    if cog_characters !=  'nan':
        if len(cog_characters) > 1:
            # delete original row
#            prokka_cog = prokka_cog.drop(index=ind)
            for letter in cog_characters:
                row_list[7] = letter
                prokka_cog2.loc[len(prokka_cog2)] = row_list
                
# remmove rows with mulitple character IDS
for ind in prokka_cog2.index:
    row_list = prokka_cog2.loc[ind, :].values.flatten().tolist()
    cog_characters = str(row_list[7])
    if cog_characters !=  'nan':
        if len(cog_characters) > 1:
            prokka_cog2 = prokka_cog2.drop(index=ind)

# sort based on locus tag to clean up
sorted_prokka_cog2 = prokka_cog2.sort_values(by='locus_tag')

# merge with cog descriptions
prokka_cog2_fun = pd.merge(sorted_prokka_cog2, cog_fun, on = 'COG_category_ID', how = 'left')

# write file
prokka_cog2_fun.to_csv('prokka-cog_multi_letters.tsv', sep="\t")