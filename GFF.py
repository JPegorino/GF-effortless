### CLASS DEFINITIONS ###
### SUPPLEMENTARY OBJECT CLASSES

### GFF contig object class
class contig:
    def __init__(self, entry, number):
        self.contig = entry
        self.name = entry.split(' ')[1]
        self.length = entry.split(' ')[3]
        self.sequence = ''
        self.number = number
        
    def __repr__(self):
        return self.name
    
    def rename(self,new_name=None,include_number=False):
        if not new_name:
            new_name = self.name
        if include_number:
            add_number = '_' + str(self.number)
            if not new_name.endswith(add_number):
                new_name = str(new_name) + add_number
        self.contig = self.contig.replace(self.name,new_name)
        self.name = new_name
        
    def concatenate_sequence(self,sequence_to_concatenate):
        self.sequence = self.sequence + sequence_to_concatenate
        
    def print_sequence(self,split_every=None):
        split_value = split_every if split_every else len(self.sequence)
        split_value = split_value
        print_info='>{}'.format(self.name)
        print_info = [print_info] + [self.sequence[pos:pos+split_value] for pos in range(0, int(self.length), split_value)] 
        return '\n'.join(print_info)
        
    def write(self):
        return self.contig + '\n'

### GFF feature object class

class GFF_entry:
    def __init__(self, entry, contig_number, ID_stat='ID', alt_ID_stat=False):
        self.raw_entry = entry
        self.entry = entry.split('\t')
        self.contig_name = self.entry[0]
        self.seq_type = self.entry[2] ; self.strand = self.entry[6]
        self.frame = int(self.entry[7]) if self.entry[7] in list('012') else 0
        self.start = int(self.entry[3]) ; self.stop = int(self.entry[4])
        self.contig_number = contig_number
        self.coords = '{};{};{}'.format(self.contig_number,self.start,self.stop) 
        self.sequence = []
        self.entry_data = {stat.split('=')[0]:stat.split('=')[1] for stat in self.entry[len(self.entry)-1].split(';')}
        if ID_stat in self.entry_data.keys():
            self.ID = self.entry_data[ID_stat]
        elif alt_ID_stat and alt_ID_stat in self.entry_data.keys():
            self.ID = self.entry_data[alt_ID_stat]
        else:
            print('Warning: no {} stat or specified alternative for entry:\n{}\nGenerating ID from co-ordinates...'.format(ID_stat,entry))
            self.ID = self.coords
        self.records = self.entry_data.keys()
        
    def __repr__(self):
        return self.ID
    
    def lookup(self,value):
        if value in self.entry_data:
            return self.entry_data[value]
        else:
            matches = []
            for stat,attribute in self.entry_data.items():
                if attribute == value:
                    matches.append(stat)
            return matches if len(matches) > 1 else matches[0]
    
    def update(self,stats_to_add,retain_old='raw_',overwrite=False):
        assert type(stats_to_add) == dict, 'stats must be provided as a dictionary'
        for stat in stats_to_add.keys():
            if stat in self.entry_data.keys():
                if not overwrite:
                    stat = retain_old + str(stat)
                    assert stat not in self.entry_data, '{}: cannot have duplicate attributes per feature'.format(stat)                    
            self.entry_data[stat] = stats_to_add.get(stat)
        self.records = self.entry_data.keys()
    
    def rename_contig(self,new_contig_name):
        assert type(new_contig_name) == str, 'contig names must be strings.'
        self.raw_entry = self.raw_entry.replace(self.contig_name,new_contig_name)
        self.entry = self.raw_entry.split('\t')
        self.contig_name = new_contig_name
            
    def add_sequence(self,entry_contig,measure_stats=False):
        assert type(entry_contig) == contig and len(entry_contig.sequence) > 0, '{} contig sequence not parsed'.format(entry_contig)
        parse_sequence = entry_contig.sequence[self.start+self.frame-1:self.stop+self.frame].upper()
        if self.strand == '-':
            parse_sequence = parse_sequence.lower().replace('a','0').replace('t','2').replace('c','1').replace('g','3')
            parse_sequence = parse_sequence.replace('0','T').replace('2','A').replace('1','G').replace('3','C').upper()[::-1]
        self.sequence.append(parse_sequence)
        if measure_stats:
            self.sequence_length = len(self.sequence)
            self.GC = 100*len(self.sequence.replace('T','').replace('A',''))/self.sequence_length
            self.contig_boundary_dist = min(min(0,min(self.start,self.stop)),entry_contig.length-max(self.start,self.stop))
        
    def print_sequence(self,fasta_name_stats='ID',split_every=None):
        if type(fasta_name_stats) == str:
            fasta_name_stats = [fasta_name_stats]
        print_info = ['>{}'.format('_'.join([self.entry_data.get(stat) for stat in fasta_name_stats if stat in self.entry_data]))]
        out_sequence = ''.join(self.sequence[0]) ### added 24/02/2023 the [0] to avoid duplicated gene sequences per fasta entry
        split_value = split_every if split_every else len(out_sequence)
        print_info = print_info + [out_sequence[pos:pos+split_value] for pos in range(0, len(out_sequence), split_value)] 
        return '\n'.join(print_info)

    def write(self,update=True):
        if update:
            entry_data_string = ';'.join([str(k) + '=' + str(v) for k,v in self.entry_data.items()])
            return '\t'.join(['\t'.join(self.entry[0:8]), entry_data_string]) + '\n'
        else:
            return self.raw_entry

### Main GFF file object class
class GFF:
    def __init__(self, file, ID_stat='ID', update_entry_stats=False, alt_ID_stat=None):
        self.file = file
        self.name = self.file.split('/')[-1].split('.gff')[0]
        self.file_info = []
        self.contigs = {}
        self.contig_count = len(self.contigs)
        self.contig_sequence = []
        self.entries = {}
        self.family = {} # a dictionary to store Parent/child links   

        with open(self.file,'r') as infile:
            all_recorded_stats = [] # a list to record each unique combination of recorded stats per gff file 
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
                    current_entry = GFF_entry(line,current_contig.number,ID_stat,alt_ID_stat) # entry object: current entry
                    # record the stats provided for each entry to generate a list of total available stats
                    if list(current_entry.records) not in all_recorded_stats:
                        all_recorded_stats.append(list(current_entry.records))
                    # add new stats from current entry to dictionary
                    if current_entry.ID in self.entries:
                        print('Warning: {} file contains multiple entries for {}, keeping first.'.format(self.file,current_entry))
                    else:
                        self.entries[current_entry.ID] = current_entry
                        if 'Parent' in current_entry.records:
                            parent_child = [current_entry.lookup('Parent'),current_entry.ID]
                            for x,y in [parent_child,parent_child[::-1]]:
                                if x in self.family.keys():
                                    if type(self.family.get(x)) == list:
                                        self.family[x] = self.family.get(x) + [y]
                                    else:
                                        self.family[x] = [self.family.get(x) ,y]
                                else:
                                    self.family[x] = y
                    # add sequence data from the end of the file, if it is there
                elif line.startswith('>'):
                    if len(self.contig_sequence) > 0: # if the sequence of a previous contig is not currently parsed and stored 
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
        # summary info
        self.all_recorded_stats = set(sum(all_recorded_stats, []))

        # update info (optional?)
        # if contig fasta sequences were provided, add sequences to entries for each locus
        if len(self.contig_sequence) > 0:
            for entry in self.entries.keys():
                current_entry = self.entries[entry]
                current_contig = self.contigs[current_entry.contig_name]
                current_entry.add_sequence(current_contig, measure_stats=update_entry_stats)
         # and finally, iterate back through the entries to add any holistic stats:
        if update_entry_stats:
            all_coords = []
            locus_index = 0
            for entry_ID,entry in self.entries.items():
                if entry.coords not in all_coords:
                    all_coords.append(entry.coords)
                    locus_index += 1
                more_info_to_add = {'index': locus_index,
                                    'family': entry.family.get(entry_ID),
                                    'sequence_length': entry.sequence_length,
                                    'contig_sequence_length': self.contigs[entry.contig_name].length,
                                    'contig_boundary_distance': entry.contig_boundary_dist,
                                    'GC': entry.GC}
                # Add ID stat to this for entries that were lacking an ID stat 
                if entry_ID == entry.coords:
                    more_info_to_add[ID_stat] = entry.coords
                self.entries[entry_ID].update(more_info_to_add)
             
    def __repr__(self):
        return str(self.name) + '.gff'
    
    def info(self,as_input=False):
        for piece_of_information in self.file_info:
            return piece_of_information if as_input else piece_of_information.lstrip('#')

    def fetch_entry_list(self):
        return [entry.ID for entry in self.entries.values()]
    
    def entry(self,entry_lookup,stat=None):
        if entry_lookup in self.entries:
            out_entry = self.entries.get(entry_lookup)
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
            out_fasta.append(locus_entry.print_sequence(list_of_fasta_header_categories,split_every))
        return '\n'.join(out_fasta)
    
    def measure_FASTA_sequence_split(self):
        sequence_line = False
        sequence_line_lengths = []
        with open(self.file,'r') as infile:
            for line in infile:
                if sequence_line:
                    sequence_line_lengths.append(len(line))
                elif line.startswith('>'):
                    sequence_line = True
        return max(sequence_line_lengths)-1
    
    def rename_contigs(self,rename_contigs=True):      
        rename_contigs = rename_contigs if type(rename_contigs) == str else self.name 
        new_names = {} ; new_contigs = {} # initialise blank dictionaries
        for contig_ID,current_contig in self.contigs.items():
            if rename_contigs in ['n','number','add_number','#']:
                current_contig.rename(new_name=False,include_number=True)
            else:
                current_contig.rename(new_name=rename_contigs,include_number=True)
            if type(contig_ID) == str:
                new_contig_name = current_contig.name
                new_names[contig_ID] = new_contig_name
            else:
                new_contig_name = contig_ID
            # creating a replacement dictionary within this loop preserves the order
            new_contigs[new_contig_name] = current_contig
        self.contigs = new_contigs
        for entry_ID,current_entry in self.entries.items():
            new_contig_name = new_names[current_entry.contig_name]
            current_entry.rename_contig(new_contig_name)            

    def to_newfile(self,out_file=False,
    number_contigs=True,rename_contigs=True,update_stats=True,
                   add_FASTA_sequence=True,FASTA_Split_Every=True,
                   skip_entries=[],skip_contigs=[]):
        out_file = out_file if out_file else self.file + '.out'
        if rename_contigs or number_contigs:
            self.rename_contigs(number_contigs,rename_contigs)
        with open(out_file,'w') as outfile:
            outfile.write(self.info(as_input=True))
            outfile.write('\n')
            for contig in self.contigs.keys():
                if self.contig_count < len(self.contigs.keys()):
                    skip_contigs = skip_contigs + [i for i in range(0,self.contig_count+1)]
                if contig not in skip_contigs:
                    outfile.write(self.contigs[contig].write())
            for entry_ID,current_entry in self.entries.items():
                corresponding_contig = current_entry.contig_name
                if corresponding_contig not in skip_contigs and entry_ID not in skip_entries:
                    outfile.write(self.entries[entry_ID].write(update=True))
            if add_FASTA_sequence:
                outfile.write('##FASTA\n')
                if FASTA_Split_Every and type(FASTA_Split_Every) != int:
                    FASTA_Split_Every = self.measure_FASTA_sequence_split()
                for contig in self.contigs.keys():
                    if contig not in skip_contigs:
                        outfile.write(self.contigs[contig].print_sequence(split_every=FASTA_Split_Every))
                        outfile.write('\n')
