# virsorter_gather_all_contigs.py
# A program to get all contigs even if they have not been annotated by virsorter and their bins

# example
# python virsorter_gather_all_contigs.py filter_A1

import os, sys, glob

## Add a check to see if there is the right number of arguments 
if len(sys.argv) != 2 :
    print("This program takes as input a directory containing the virsorter results")
    sys.exit(2)
    
virsorter_dir = sys.argv[1];
# move into the directory containing the sequences
os.chdir(virsorter_dir)

# Create file for results
outfile = open('virsorter_allcontigs.tsv', 'w')

# move into the directory containing the sequences
os.chdir("Bins")

# get a list of the output directories
fasta_list = glob.glob('*.fasta')

# write column headers
outfile.write('bin\tseqname\t%GC\n')

for file in fasta_list :
    current_fasta_ID = ''
    dict_fasta = {}          # Make an empty dictionary of genes
    final_score = open(file, 'r')
    file_contents = final_score.readlines()
    bin_name = file.replace('.fasta', '')
    for line in file_contents :
        line = line.strip()
        if line.startswith(">") :
            line = line.replace('>', '')
            current_fasta_ID = line
            dict_fasta[line] = ""
        else :
            dict_fasta[current_fasta_ID] += line
    for fasta_ID in dict_fasta:
        seq = dict_fasta[fasta_ID]
        GC = (seq.count('G') + seq.count('C')) / len(seq) * 100
        outfile.write('%s\t%s\t%.0f\n' % (bin_name, fasta_ID, GC))
            
            
# Move back to your virsorter directory
os.chdir("..")
os.chdir("..")