# hallmark-genes.GV.py
# A program to run virsorter on multiple fasta files in the same directory

# open the outfile
outfile = open('hallmark.gv.faa', 'w')

# open hallmark genes
infile1 = open('hallmark-gene.list.GV.txt', 'r')
hallmark_list = infile1.readlines()

# open virus gene file
infile2 = open('filter_A1/Bins/bin.049.fasta.out/iter-0/all.pdg.Viruses.hmmtbl', 'r')
virus_list = infile2.readlines()

# open faa file
infile3 = open('filter_A1/Bins/bin.049.fasta.out/iter-0/all.pdg.faa', 'r')
faa_list = infile3.readlines()

# make a dictionary of the virus gene file with ID and contig name

dict_virus = {}       
for line in virus_list :
    if not line.startswith('#') :
        line_list = line.split(' -          ')
        contig_ID = line_list[0]
        contig_ID = contig_ID.strip()
        virus_ID = line_list[1]
        virus_ID = virus_ID.strip()
        dict_virus[contig_ID] = virus_ID

# make a dictionary of contig IDs with hallmark genes and make a dictionary with annotation

dict_hallmark = {} 
dict_annotation = {}
for gene in hallmark_list :
    gene = gene.strip()
    gene_list = gene.split('\t')
    ID = gene_list[0]
    annotation = gene_list[1]
    dict_annotation[ID] = annotation
    for key in dict_virus :
        virus_ID = dict_virus[key]
        if virus_ID == ID :
            dict_hallmark[key] = ID
            
# get fasta sequence for each hallmark gene

dict_faa = {} 
current_faa_ID = ''
for line in faa_list :
    line = line.strip()
    if line.startswith(">") :
        line = line.replace('>', '')
        line_list = line.split(' # ')
        current_faa_ID = line_list[0]
        dict_faa[current_faa_ID] = ""
    else :
        dict_faa[current_faa_ID] += line

#  Write out a file with the file mark gee

for key in dict_hallmark :
    virus_gene = dict_hallmark[key]
    protein = dict_faa[key]
    annotation = dict_annotation[virus_gene]
    outfile.write('>%s|%s|%s\n' % (key, virus_gene, annotation))
    outfile.write('%s\n' % (protein))