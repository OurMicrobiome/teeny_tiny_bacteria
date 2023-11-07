# A program to merge the depth, bins and results files
# It will make suggested viral calls

# Example
# python virsorter_detect_viral_bins.py filter_A1

import os, sys
import pandas as pd
import numpy as np

## Add a check to see if there is the right number of arguments 
if len(sys.argv) != 2 :
    print("This program takes as input a directory containing the virsorter results")
    sys.exit(2)
    
virsorter_dir = sys.argv[1];
# move into the directory containing the sequences
os.chdir(virsorter_dir)

# read in files
depth = pd.read_csv('merged.depth.txt',sep = '\t') 
depth.rename(columns = {'contigName':'seqname'}, inplace = True)

bins = pd.read_csv('virsorter_allcontigs.tsv',sep = '\t')

results = pd.read_csv('virsorter_results_cat.tsv',sep = '\t') 
results[['seqname', 'full-partial']] = results['seqname'].str.split('|', 1, expand=True)

# merge dataframes
depth_bins = pd.merge(depth, bins, on = 'seqname', how = 'left')
depth_bins_results = pd.merge(depth_bins, results, on = 'seqname', how = 'left')

# sort merged data frame
depth_bins_results.rename(columns = {'bin_x':'bin'}, inplace = True)
depth_bins_results = depth_bins_results.sort_values(['bin', 'contigLen'], ascending = (True, False))

# Replace NA in max group score
depth_bins_results.max_score_group = depth_bins_results.max_score_group.fillna('Bacteria?')

# filter data frame to a certain contig length
depth_bins_results_mask=depth_bins_results['contigLen']>=5000
depth_bins_results_5k = depth_bins_results[depth_bins_results_mask]

# print dataframe
depth_bins_results_5k.to_csv(virsorter_dir + '_depth_bin_virsorter.tsv', sep="\t")

# Get length, hallmark gene sum and GC mean of each bin

bin_length = depth_bins_results_5k.groupby(['bin']).agg(BinContigLenSum=('contigLen','sum'), 
                                            hallmarkSum = ('hallmark','sum'),
                                            gcMean = ('%GC','mean'))
bin_length['hallmarkSum'] = bin_length['hallmarkSum'].astype(int)
bin_length['gcMean'] = bin_length['gcMean'].astype(int)
#bin_length = depth_bins_results_5k.groupby(['bin'])['contigLen'].sum().reset_index()
bin_length.to_csv('bin_length.tsv', sep="\t")

# Get the largest contig in each bin
# print(depth_bins_results.groupby(['bin'])['contigLen'].max().reset_index())

# Get the largest contig row with the largest contig for each bin
contig_calls = depth_bins_results_5k.loc[depth_bins_results_5k.groupby(['bin'])['contigLen'].idxmax()]
contig_calls_select = contig_calls.loc[:, ('seqname', 'bin', 'contigLen', 'max_score', 'max_score_group', 'full-partial')]
contig_calls_select.rename(columns = {'max_score_group':'max_contigLenGroup'}, inplace = True)
contig_calls_select.rename(columns = {'contigLen':'max_contigLen'}, inplace = True)
contig_calls_select.to_csv(virsorter_dir + '_contig_calls.tsv', sep="\t")

# Calculate total contig length and percentages of contig length for each type for each bin

max_score_group = depth_bins_results_5k.groupby(['bin', 'max_score_group']).agg(bin_high_ContigLenSum=('contigLen','sum'))
max_score_group['GroupPercent'] = max_score_group.groupby(level=0).apply(lambda x:
                                                 100 * x / float(x.sum()))
max_score_group['GroupPercent'] = max_score_group['GroupPercent'].astype(int)
max_score_group.to_csv(virsorter_dir + '_max_score_group.tsv', sep="\t")

# Select the group with the highest percent total contig length as the main group for the bin

bin_calls_pre = max_score_group.loc[max_score_group.groupby(['bin'])['bin_high_ContigLenSum'].idxmax()]
bin_calls_pre.to_csv('bin_calls_pre.tsv', sep="\t")

# remove the unbinned from the bin calls
# I can't seem to use bin_calls1 dataframe, but I can read it from the file and then remove. Still learning pandas

bin_calls = pd.read_csv('bin_calls_pre.tsv',sep = '\t')
bin_calls.rename(columns = {'max_score_group':'bin_high_group'}, inplace = True)
bin_calls['GroupPercent'] = bin_calls['GroupPercent'].astype(int)
bin_calls['bin_high_ContigLenSum'] = bin_calls['bin_high_ContigLenSum'].astype(int)
bin_calls = bin_calls[bin_calls.bin != 'bin.unbinned']
bin_calls.to_csv(virsorter_dir + '_bin_calls.tsv', sep="\t")

# merge contig calls and percent total calls

final_calls1 = pd.merge(contig_calls_select, bin_calls, on = 'bin', how = 'left')
final_calls = pd.merge(bin_length, final_calls1, on = 'bin', how = 'left')
final_calls.to_csv(virsorter_dir + '_final_calls.tsv', sep="\t")
