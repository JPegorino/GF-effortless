#!/usr/bin/env python3
import sys, os, io, re, argparse, pandas as pd, numpy as np
sys.path.append('os.path(__file__)')
import GFF # GFF parser

def parse_args():
    parser = argparse.ArgumentParser(
        description="Add e.g. orthogroup statistics from a data table to genome GFFs.")
    parser.add_argument("data_table", 
        help="Path to the input data table (CSV or TSV).")
    parser.add_argument("gff_input_file", 
        help="input GFF file.")
    parser.add_argument("-s", "--stat_name",
        type=str, default="COG",
        help="Optional label to prefix stat columns. Default: 'COG'.")
    parser.add_argument("-l", "--last_column_to_add",
        type=int, default=1,
        help="Column number of final column to include, (e.g. before the first genome column in a pangenome table). Default: Add first column only.")
    return parser.parse_args()

# collect user arguments
args = parse_args()
data_table = args.data_table
gff_input_file = args.gff_input_file
stat_name = args.stat_name
last_column_to_add = args.last_column_to_add

# define functions
def drop_special_chars(string):
    for character in list('() []{}\\|/",;:`'+"'"):
        string = string.replace(character,'')
    return string

# load data
assert os.path.isfile(gff_input_file), "requires an input gff file."
assert str(gff_input_file).endswith('gff'), "input gff file must be a .gff file."
gff_id = os.path.splitext(os.path.basename(gff_input_file))[0]
input_gff = GFF.GFF(gff_input_file,update_feature_stats=True)
assert os.path.isfile(data_table), "requires an input file with tabular data per CDS for each gff (columns)."
delim = ',' if data_table.endswith('csv') else '\t'
df = pd.read_csv(data_table, sep=delim)
pangenome_stats_to_add = df.columns.to_list()

# update out list of stats to add to include the stat_name identifier
if last_column_to_add: # if a column number corresponding to the first genome column was provided
    pangenome_stats_to_add = pangenome_stats_to_add[0:last_column_to_add] # use it to filter the table
else: # otherwise, try to infer from the gff files in the directory (which is likely to include all of them)
    pangenome_stats_to_add = [stat for stat in pangenome_stats_to_add if stat != gff_id]
pangenome_stats_to_add[0] = stat_name
pangenome_stats_to_add[1:] = ['_'.join([stat_name,drop_special_chars(stat)]) for stat in pangenome_stats_to_add[1:]]
old_columns = df.columns.to_list()[0:len(pangenome_stats_to_add)]
df.rename(columns=dict(zip(old_columns,pangenome_stats_to_add)),inplace=True)

# process data for output file
pangenome_stats_to_add.append(gff_id) # add the genome coulmnd for the gff file
df_id = df[pangenome_stats_to_add] # drop columns other than those with stats we want to add
df_id = df_id[df_id[gff_id].notna() & (df_id[gff_id] != "")] # filter out rows with accessory genes that are absent in each genome
# gff_to_COG = df_id.set_index(gff_id)[df.columns[0]].to_dict() # for just column 1 (typically the COG names)
nested_dict = df_id.set_index(gff_id).to_dict(orient="index") # for all the COG stats
feature_count = 0
for feature,stats in nested_dict.items(): # iterate through the table values to add stats for
    try:
        input_gff.features.get(feature).update(stats_to_add=stats)
        feature_count +=1
    except:
        for split_feature in feature.replace(':',';').split(';'):
            try:
                input_gff.features.get(drop_special_chars(split_feature)).update(stats_to_add=stats)
                feature_count +=1
            except:
                print('warning! Could not parse stats for {} from {}'.format(feature,data_table))

print('stats added to {} of {} features in {}'.format(feature_count,len(input_gff.features.keys()),gff_input_file))
# Save updated GFF to new file
input_gff.to_newfile(out_file=gff_input_file.replace('.gff','_{}-data.gff'.format(stat_name)))