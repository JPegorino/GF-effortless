import sys, os, io, re, warnings, argparse
sys.path.append('os.path(__file__)')
import GFF # GFF parser

def parse_args():
    parser = argparse.ArgumentParser(
        description="Search a GFF for features and extract feature data.")
    parser.add_argument("gff_input", 
        help="Path to your GFF file to search.")
    parser.add_argument("search", 
        help="Search criteria.")
    parser.add_argument("-o", "--output_format", 
        help="Output format. Must be one of: TBD")
    parser.add_argument("-d", "--output_delimiter", 
        default="\t",
        help="Output delimiter.")
    parser.add_argument("-fh", "--fasta_header",
        default="ID",
        help="A semi-colon delimited subset of stats to display, or to use for generating the FASTA header. Default: 'ID'.")
    parser.add_argument("-us", "--upstream",
        default=0,
        help="The length (bp) of upstream sequence to include (must be a multiple of 3 for protein). Default: '0'.")
    parser.add_argument("-ds", "--downstream",
        default=0,
        help="The length (bp) of downstream sequence to include (must be a multiple of 3 for protein). Default: '0'.")
    parser.add_argument("-fs", "--fasta_split_every",
        default=60,
        help="The number of characters per line to include in FASTA sequence output (must be an integer). Default: '60'.")
    return parser.parse_args()

# collect user arguments
args = parse_args()
gff_input = GFF.GFF(args.gff_input.rstrip('/'),update_feature_stats=True)
search_feature_info = args.search
out_format = args.output_format
out_delimiter = args.output_delimiter
fasta_header = args.fasta_header 
try:
    us = int(args.upstream)
    ds = int(args.downstream)
    split_every = int(args.fasta_split_every)
except:
    print('upstream and downstream values must be numeric integers')
    sys.exit

out_feature = gff_input.feature(search_feature_info)
if not out_feature:
    out_feature = gff_input.feature(search_feature_info,strictly_first=False,regex=True)
if not out_feature:
    print('Error: Nothing matched.')
    sys.exit(1)

if type(out_feature) != list:
    out_feature = [out_feature]

for my_out in out_feature:
    if not out_format:
        print(my_out.ID)
    elif out_format == 'coords':
        print(my_out.coords)
    elif out_format == 'features':
        print(my_out.feature_info)
    elif out_format == 'number' or out_format == 'index':
        print(my_out.idx)
    elif out_format == 'bed':
        print(my_out.write(bed=True))
    elif out_format == 'gff':
        print(my_out.raw_entry)
    elif out_format == 'stats':
        for key,value in my_out.feature_info.items():
            print(f'{out_delimiter}'.join([key,str(value)]))
    elif out_format == 'subset':
        for stat in fasta_header.split(';'):
            if stat in my_out.feature_info:
                print(f'{out_delimiter}'.join([stat,str(my_out.feature_info.get(value))]))
    elif out_format == 'tab' or out_format == 'table':
        print(f'{out_delimiter}'.join(my_out.feature_info.keys()))
        print(f'{out_delimiter}'.join([str(val) for val in my_out.feature_info.values()]))
    elif out_format == 'subtab' or out_format == 'subset_table':
        stat_subset = {key:my_out.feature_info.get(key) for key in fasta_header.split(';') if key in my_out.feature_info} 
        print(f'{out_delimiter}'.join(key for key in stat_subset.keys()))
        print(f'{out_delimiter}'.join([str(val) for val in stat_subset.values()]))
    elif out_format == 'fasta' or out_format == 'ffn':    
        print(my_out.print_sequence(gff_input.contigs.get(my_out.contig_name),split_every=60,fasta_name_stats=fasta_header.split(';'),us=us,ds=ds))
    elif out_format == 'protein' or out_format == 'faa':    
        print(my_out.print_sequence(gff_input.contigs.get(my_out.contig_name),fasta_header.split(';'),split_every=60,us=us,ds=ds,protein=True))
    elif out_format == 'contig' or out_format == 'region':    
        out_contig = gff_input.contigs.get(my_out.contig_name)
        reverse_strand = True if my_out.strand == '-' else False
        print(out_contig.print_sequence(reverse_strand=reverse_strand,split_every=60,region=(us,ds),rename=fasta_header))
    elif out_format in my_out.feature_info:
        print(my_out.feature_info.get(out_format))
    else:
        print(f'Error: No info {out_format} for {my_out}')
        continue