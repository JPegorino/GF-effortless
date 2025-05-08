#!/usr/bin/env python3
import sys, os, io, re, argparse, pandas as pd, numpy as np
sys.path.append('os.path(__file__)')
import GFF # GFF parser

def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract nucleotide sequences from genome GFFs.")
    parser.add_argument("gff_input", 
        help="Path to your GFF file from which to extract sequences.")
    parser.add_argument("output_file_name", 
        help="Name to use for the file in which to write otuput sequence data.")
    parser.add_argument("-l","--lookup", 
        help="A value or file containing values to extract from the data.")
    parser.add_argument("-s", "--stat_name",
        default="ID",
        help="The name of a GFF file statistic to look up and compare with string. Default: 'ID'.")
    parser.add_argument("-f", "--filtering",
        default=None,
        help="For numeric columns, specifies the arithmetic operation to use for filtering, and the filtering value.\n\
        Input must be three space-delimited strings in the format 'feature stat' 'operation' 'value'\n\
        Operation must be one of '==', '!=', '>=', or '<='.")
    return parser.parse_args()

# collect user arguments
args = parse_args()
gff_input = GFF.GFF(args.gff_input.rstrip('/'))
output_file_name = args.output_file_name
lookup = args.lookup
stat_name = args.stat_name
filtering = args.filtering

# parse data for user arguments and handle poor input values
if stat_name:
    assert stat_name in gff_input.all_recorded_stats, "Error: {} information is not recorded for gff features.".format(stat_name)

if filtering:
    assert filtering[0] in gff_input.all_recorded_stats, "{} information is not recorded for gff features.".format(filtering[0])
    assert filtering[1] in ['==','!=','>=','<='], "Operation must be one of ['==', '!=', '>=', or '<=']".format(filtering[2])
    assert type(filtering[2]) == int, "{} is not a numeric integer".format(filtering[2])

# parse input attritbute list
if os.path.isfile(lookup):
    with open(lookup,'r') as gene_list:
        ids = [l.rstrip('\r').rstrip('\n') for l in lookup_list]
else:
    ids = [lookup] # to look up a single value

# identify features with matching attributes
selected_features = []
for id in ids:
    if stat_name == 'ID':
        selected_feature = gff_input.feature(id,feature_type='CDS')
        selected_features.append(selected_feature)
    else: # elif stat_name: - can this be improved ??
        for progenitor_ID in gff_input.families.keys():
            feature_with_info = gff_input.feature(progenitor_ID).lookup(id)
            if feature_with_info:
                selected_feature = gff_input.feature(progenitor_ID,feature_type=feature_with_info)
                selected_features.append(selected_feature)

# filtering features
if filtering:
    filtered_features = []
    for feature in selected_features:
        filter_data = feature.lookup(filtering[0])
        if filtering[1] == '==':
            if filter_data == filtering[2]:
                filtered_features.append(feature)
        elif filter[1] == '!=':
            if filtering_data != filtering[2]:
                filtered_features.append(feature)
        elif filtering[1] == '>=':
            if filter_data >= filtering[2]:
                filtered_features.append(feature)
        elif filtering[1] == '<=':
            if filter_data <= filtering[2]:
                filtered_features.append(feature)
    selected_features = filtered_features

if len(selected_features) == 0:
    print('no matching {} features for {}- after filtering nothing to write to file. Exiting...'.format(stat_name,lookup))
    sys.exit(0)

# extract nucleotide sequences
with open(output_file_name,'w') as outfile:
    for feature in selected_features:
        sequence = feature.print_sequence(split_every=60)
        outfile.write(sequence + '\n')
