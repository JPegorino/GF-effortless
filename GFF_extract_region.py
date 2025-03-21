#!/usr/bin/env python3
import sys, os, io, re, warnings
sys.path.append('os.path(__file__)')
import GFF # GFF parser

start = sys.argv[1]
end = sys.argv[2]
gff_input = GFF.GFF(sys.argv[3].rstrip('/'))
output_file_name = sys.argv[4].rstrip('/')
include_upstream = int(sys.argv[5]) if sys.argv[5] else 0
use_name = sys.argv[6] if sys.argv[6] else False

with open(output_file_name,'w') as outfile:
    start_feature = gff_input.feature(start,stat='ID')
    end_feature = gff_input.feature(end,stat='ID')
    assert start_feature.contig_number == end_feature.contig_number, 'Error: cannot extract region spanning features on different contigs.'
    if start_feature.start < end_feature.start:
        if start_feature.strand != '+' and include_upstream > 0:
            include_upstream = 0
            warnings.warn('Warning: End feature is upstream of start feature, consider reversing.',Warning)
        # correct orientation
        reverse_strand = False
        start = start_feature.start - include_upstream - 1
        end = end_feature.stop
    else:
        reverse_strand = True
        if start_feature.strand != '-' and include_upstream > 0:
            include_upstream = 0
            warnings.warn('Warning: End feature is upstream of start feature, consider reversing.',Warning)
        start = end_feature.start - 1
        end = start_feature.stop + include_upstream
    contig = gff_input.contigs.get(start_feature.contig_name)
    region_sequence = contig.print_sequence(reverse_strand=reverse_strand,split_every=60,region=(start,end),rename=use_name)
    outfile.write(region_sequence + '\n')