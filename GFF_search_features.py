import sys, os, io, re, warnings
sys.path.append('os.path(__file__)')
import GFF # GFF parser

gff = GFF.GFF(sys.argv[1].rstrip('/'))
out_feature = gff.feature(sys.argv[2])
if not out_feature:
    out_feature = gff.feature(sys.argv[2],strictly_first=False,regex=True)
if not out_feature:
    print('Error: Nothing matched.')
    sys.exit(1)

if type(out_feature) != list:
    out_feature = [out_feature]

for my_out in out_feature:
    print(my_out.ID)
    print(my_out.coords)
    print(my_out.feature_info)
    print(my_out.print_sequence())