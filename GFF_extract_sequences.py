#!/usr/bin/env python3
import sys, os, io, re, pandas as pd, numpy as np
sys.path.append('os.path(__file__)')
import GFF # GFF parser

gene_list_file = sys.argv[1]
gff_input = GFF.GFF(sys.argv[2].rstrip('/'))
output_file_name = sys.argv[3].rstrip('/')

with open(gene_list_file,'r') as gene_list:
    ids = [l.rstrip('\r').rstrip('\n') for l in gene_list]

with open(output_file_name,'w') as outfile:
    for id in ids:
        sequence = gff_input.feature(id,stat='ID').printing_sequence(split_every=60)
        outfile.write(sequence + '\n')
