### CLASS DEFINITIONS ###
### SUPPLEMENTARY OBJECT CLASSES

class contig:
    def __init__(self, entry, number):
        self.name = entry.split(' ')[1]
        self.length = entry.split(' ')[3]
        self.sequence = ''
        self.number = number
    def __repr__(self):
        return self.name
    def concatenate_sequence(self,sequence_to_concatenate):
        self.sequence = self.sequence + sequence_to_concatenate

class GFF_entry:
    def __init__(self, entry, contig_number=None):
        self.entry = entry.split('\t')
        self.contig_name = self.entry[0]
        self.seq_type = self.entry[2] ; self.strand = self.entry[6]
        self.frame = int(self.entry[7]) if self.entry[7] in list('012') else 0
        self.start = int(self.entry[3]) ; self.stop = int(self.entry[4])
        self.contig = contig_number if contig_number else self.contig_name
        self.sequence = []
        self.name = '{};{};{}'.format(self.contig_name,self.start,self.stop)
        self.UID = '{};{};{}'.format(self.contig,self.start,self.stop)
        self.entry_data = {stat.split('=')[0]:stat.split('=')[1] for stat in self.entry[len(self.entry)-1].split(';')}
        self.records = self.entry_data.keys()
    def __repr__(self):
        return self.UID
    def stat_info(self,stat):
        return self.entry_data[stat] if stat in self.entry_data else None
    def add_entry_records(self,stats_to_add):
        assert type(stats_to_add) == dict, 'stats must be provided as a dictionary'
        for stat in stats_to_add.keys():
            if stat not in self.entry_data:
                self.entry_data[stat] = stats_to_add.get(stat)
            else: # if a stat is recorded for the locus twice, split data as gene/CDS
                gene_stat = '{}_{}'.format('gene',stat)
                if self.seq_type == 'CDS':
                    self.entry_data[gene_stat] = stats_to_add.get(stat)
                else:
                    self.entry_data[gene_stat] = self.entry_data.get(stat)
                    self.entry_data[stat] = stats_to_add.get(stat)
    def entry_record_by_stat(self,stat,dictionary):
        assert type(dictionary) == dict, 'input dictionary not a dictionary'
        if self.stat_info(stat): # relies on None == False
            dictionary[self.stat_info(stat)] = self.UID
    def add_sequence(self,entry_contig):
        assert type(entry_contig) == contig and len(entry_contig.sequence) > 0, 'contig sequence not parsed'
        parse_sequence = entry_contig.sequence[self.start+self.frame-1:self.stop+self.frame]
        if self.strand == '-':
            parse_sequence = parse_sequence.lower().replace('a','0').replace('t','2').replace('c','1').replace('g','3')
            parse_sequence = parse_sequence.replace('0','T').replace('2','A').replace('1','G').replace('3','C').upper()[::-1]
        self.sequence.append(parse_sequence)
    def printing_sequence(self,fasta_name_stats='ID',split_every=None):
        if type(fasta_name_stats) == str:
            fasta_name_stats = [fasta_name_stats]
        print_info = ['>{}'.format('_'.join([self.entry_data.get(stat) for stat in fasta_name_stats if stat in self.entry_data]))]
        out_sequence = ''.join(self.sequence[0]) ### added 24/02/2023 the [0] to avoid duplicated gene sequences per fasta entry
        split_value = split_every if split_every else len(out_sequence)
        print_info = print_info + [out_sequence[pos:pos+split_value] for pos in range(0, len(out_sequence), split_value)] 
        return '\n'.join(print_info)

### MAIN GFF CLASS

class GFF:
    def __init__(self, file, other_stat=None):
        self.file = file
        self.name  = self.file.split('/')[-1]
        self.file_info = []
        self.contigs = {}
        self.contig_count = len(self.contigs)
        self.contig_sequence = []
        self.entries = {}
        self.entries_uid_only = {}
        self.entries_per_loci = {}
        self.entries_per_stat = {} # another dictionary can be stored for one other stat at a time
        
        with open(self.file,'r') as infile:
            all_recorded_stats = [] # a list to record each unique combination of recorded stats per entry 
            for line in infile:
                line = line.rstrip('\n')
                if line.startswith('##'):
                    if line.startswith('##sequence-region'):
                        self.contig_count += 1
                        current_contig = contig(line,self.contig_count)
                        self.contigs[str(current_contig)] = current_contig
                        if current_contig.name != str(current_contig.number):
                            self.contigs[current_contig.number] = current_contig
                    else:
                        self.file_info.append(line)
                elif line.split('\t')[0] in self.contigs:
                    current_contig = self.contigs[line.split('\t')[0]] # contig object: current contig
                    current_entry = GFF_entry(line,current_contig.number) # entry object: current entry
                    # record the stats provided for each entry to generate a list of total available stats 
                    if list(current_entry.records) not in all_recorded_stats:
                        all_recorded_stats.append(list(current_entry.records))
                    # add new stats from current entry to dictionary, both by unique identifier and per locus tag
                    if current_entry.UID in self.entries:
                        self.entries[current_entry.UID].add_entry_records(current_entry.entry_data)
                    else:    
                        self.entries[current_entry.UID] = current_entry
                        self.entries_uid_only[current_entry.UID] = current_entry
                        if current_entry.name not in self.entries:
                            self.entries[current_entry.name] = current_entry
                    if 'locus_tag' in self.entries[current_entry.UID].entry_data:
                        current_entry.entry_record_by_stat(stat='locus_tag',dictionary=self.entries_per_loci)
                elif line.startswith('>'):
                    if len(self.contig_sequence) > 0: # if the sequence of a previous contig is currently parsed and stored 
                        current_contig.concatenate_sequence(''.join(self.contig_sequence)) # add it to the data for the contig
                        self.contig_sequence = [] # and wipe the storage variable to start parsing the next contig
                    line = line.lstrip('>')
                    for character in range(len(line)):
                        if line in self.contigs or len(line) == 0:
                            if len(line) == 0:
                                print('Warning: No contig name in fasta entry. Contig fasta sequences were incorrectly parsed.')
                            current_contig = self.contigs[line] ; break
                        else:
                            line = line[0:len(line)-1]
                else:
                    self.contig_sequence.append(line)
            current_contig.concatenate_sequence(''.join(self.contig_sequence)) # store the parsed data for the final contig
        
        self.all_recorded_stats = set(sum(all_recorded_stats, []))
         # if contig fasta sequences were provided, add sequences to entries for each locus
        if len(self.contig_sequence) > 0:
            for entry in self.entries.keys():
                current_entry = self.entries[entry]
                current_contig = self.contigs[current_entry.contig_name]
                current_entry.add_sequence(current_contig)
    
    def print_file_info(self):
        for piece_of_information in self.file_info:
            print(piece_of_information.lstrip('#'))

    def fetch_entry_list(self):
        return [entry.UID for entry in self.entries.values()]
    
    def entry(self,entry_lookup,stat=None):
        if entry_lookup in self.entries:
            out_entry = self.entries.get(entry_lookup)
        elif entry_lookup in self.entries_per_loci:
            out_entry = self.entries.get(self.entries_per_loci.get(entry_lookup))
        elif entry_lookup in self.entries_per_stat:
            out_entry = self.entries.get(self.entries_per_stat.get(entry_lookup))
        else:
            assert stat in self.all_recorded_stats, 'file contains no information for "{}"'.format(entry_lookup)
            for entry in self.entries.values():
                entry.entry_record_by_stat(stat=stat,dictionary=self.entries_per_stat)
            out_entry = self.entries.get(self.entries_per_stat.get(entry_lookup))
        return out_entry 

    def fetch_entry_sequences(self,list_of_loci,list_of_fasta_header_categories,split_every=None):
        out_fasta = []
        if type(list_of_loci) == str:
            list_of_loci = [list_of_loci]
        for locus in list_of_loci:
            locus_entry = self.entry(locus)
            out_fasta.append(locus_entry.printing_sequence(list_of_fasta_header_categories,split_every))
        return '\n'.join(out_fasta)
        
### SUPPLEMENTARY CLASS FOR BWA MAPPING DATA ###

class bwa_mapping:
    def __init__(self, mapping):
        self.number_of_mapped_loci = len(mapping.split(','))
        if self.number_of_mapped_loci == 1 and (len(mapping.split(';')) != 4 or len(mapping.split(':')) <= 1):
            print('Warning: either {} is not a bwa-mapping or is incorrectly ";" delimited'.format(mapping))
            mapping = str(mapping) + ';;;'
        self.contig = mapping.split(';')[0] ; self.gene = mapping.split(';')[2]
        self.prev = mapping.split(';')[1] ; self.next = mapping.split(';')[3]
        self.genome = ('_').join(self.contig.split(':')[0].split('|')[-1].split('_')[:-1])
    
    def contig_name_to_UID(self,numbered=False):
        contig_name = self.contig.split(':')[0].split('_')[-1] if numbered else self.contig.split(':')[0]
        hit_position = self.contig.split(':')[1].replace('-',';')
        return str(contig_name) + ';' + str(hit_position)
    
    def location(self):
        if self.number_of_mapped_loci > 1:
            loc = 'multiple loci'
        elif self.gene:
            if len(self.gene.split('|')) > 1:
                loc = 'overlapping adjacent loci'
            elif self.gene == self.prev == self.next:
                loc = 'CDS'
            else:
                loc = 'contig end'
        else:
            loc = 'intergenic' if self.prev and self.next else 'contig end'
        return loc
    
    def cds_entry(self,entry_lookup,gff,stat=None):
        lookups = entry_lookup.split(';') ; lookup_out = None
        for entry in gff.entries.keys():
            if not entry.startswith(lookups[0]):
                continue # skips entries on other contigs to save time
            current_entry = gff.entries[entry]
             # some entries cover slightly more than the CDS, so identify from the start only (not the stop)
            if (current_entry.start <= int(lookups[1])) and (current_entry.stop >= int(lookups[1])):
                lookup_out = current_entry.UID
                if (current_entry.stop <= int(lookups[2])):
                    print('{} goes over {} gene boundary'.format(entry_lookup,current_entry.UID))
                break
        return gff.entry(entry_lookup=lookup_out,stat=stat) if lookup_out else None

### FUNCTIONS FOR PROCESSING DEFINED OBJECTS ###
### CONVERSION FUNCTION FOR LOCUS TAGS OF IDENTICAL GFF FILES

def GFF_convert_locus_tags(f1,f2,locus_tag):
    f1 = GFF(f1); f2 = GFF(f2)
    locus_uid = f1.entry(locus_tag).UID # parse input
    assert locus_uid in f2.entries, 'comparison gff file missing locus.'
    return '\t'.join([locus_tag,f2.entry(locus_uid).stat_info('locus_tag')])

 # define a function that takes a unique id for a genomic locus (its co-ordinates in an assembly) and returns the locus tag from a gff file
def lookup_locus_tag(gff,locus_uid):
    try:
        locus_tag = gff.entry(locus_uid).stat_info('locus_tag')
    except:
        locus_tag = 'no data' # if the lookup fails, the function returns 'no data'
    return locus_tag

