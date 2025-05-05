#!/usr/bin/env python3
import sys, os, io, re, pandas as pd, numpy as np
sys.path.append('os.path(__file__)')
import GFF # GFF parser

data_table = sys.argv[1]
gff_input_dir = sys.argv[2]
try:
    stat_name = str(sys.argv[3])
except: 
    stat_name = "COG" # assumes COG if no alternate stat name provided

# define functions
def drop_special_chars(string):
    for character in list('() []{}\\|/",;:`'+"'"):
        string = string.replace(character,'')
    return string

# load data
assert os.path.isdir(gff_input_dir), "requires an input directory with gff files."
assert os.path.isfile(data_table), "requires an input file with tabular data per CDS for each gff (columns)."
delim = ',' if data_table.endswith('csv') else '\t'
df = pd.read_csv(data_table, sep=delim)

# define important variables
pangenome_stats_to_add = df.columns.to_list()
genome_list = {}

# extract inputs
for file in os.listdir(gff_input_dir):
    gff_id = os.path.splitext(file)[0]
    file = os.path.join(gff_input_dir,file)
    if not os.path.isfile(file) or not str(file).endswith('gff'):
        print('{} was not detected as a gff file. Skipping...'.format(file))
        continue
    pangenome_stats_to_add.remove(gff_id) # drop the genome from the list of stats to add
    genome_list[gff_id] = file # add the data for the gff to the dictionary

# update out list of stats to add to include the stat_name identifier
pangenome_stats_to_add[0] = stat_name
pangenome_stats_to_add[1:] = ['_'.join([stat_name,drop_special_chars(stat)]) for stat in pangenome_stats_to_add[1:]]
old_columns = df.columns.to_list()[0:len(pangenome_stats_to_add)]
df.rename(columns=dict(zip(old_columns,pangenome_stats_to_add)),inplace=True)

for gff_id,file in genome_list.items():
    input_gff = GFF.GFF(file,update_feature_stats=True)
    df_id = df[df[gff_id].notna() & (df[gff_id] != "")] # filter out rows with accessory genes that are absent in each genome
    subset_columns = pangenome_stats_to_add
    subset_columns.append(gff_id)
    df_id = df_id[subset_columns] # drop columns other than those with stats we want to add
    # gff_to_COG = df_id.set_index(gff_id)[df.columns[0]].to_dict() # for just the COG names
    nested_dict = df_id.set_index(gff_id).to_dict(orient="index") # for all the COG stats
    feature_count = 0
    for feature,stats in nested_dict.items(): # iterate through the table values to add stats for
        try:
            input_gff.features.get(feature).update(stats_to_add=stats)
            feature_count +=1
        except:
            for split_feature in feature.split(';'):
                try:
                    input_gff.features.get(drop_special_chars(split_feature)).update(stats_to_add=stats)
                    feature_count +=1
                except:
                    print('warning! Could not parse stats for {} from {}'.format(feature,data_table))
    
    print('stats added to {} of {} features in {}'.format(feature_count,len(input_gff.features.keys()),file))
    # Save updated GFF to new file
    input_gff.to_newfile(out_file=file.replace('.gff','_{}-data.gff'.format(stat_name)))
    