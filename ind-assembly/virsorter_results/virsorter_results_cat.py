# virsorter_resulter_cat.py
# A program to concatenate results from genome bins

# example
# python virsorter_results_cat.py filter_A1

import os, sys, glob

## Add a check to see if there is the right number of arguments 
if len(sys.argv) != 2 :
    print("This program takes as input a directory containing the virsorter results")
    sys.exit(2)
    
virsorter_dir = sys.argv[1];
# move into the directory containing the sequences
os.chdir(virsorter_dir)

# Create file for results
outfile = open('virsorter_results_cat.tsv', 'w')
    
# move into the directory containing the sequences
os.chdir("Bins")

# get a list of the output directories
out_dir_list = glob.glob('*.out')

flag = "no"
# run virsorter on each of these files
for out_dir in out_dir_list :
    os.chdir(out_dir)
    final_score = open('final-viral-score.tsv', 'r')
    col_names = final_score.readline()
    bin_name = out_dir.replace('.fasta.out', '')
    if flag == "no" :
        outfile.write('bin\t%s' % (col_names))
        flag = "yes"
    file_contents = final_score.readlines()
    for line in file_contents :
        outfile.write('%s\t%s' % (bin_name, line))
    os.chdir("..")
    
# Move back to your virsorter directory
os.chdir("..")
os.chdir("..")