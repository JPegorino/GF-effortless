### DEPENDENCIES ###
import re
### CLASS DEFINITIONS ###

### GFF contig object class
class contig:
    def __init__(self, entry, number):
        self.contig = entry
        self.name = entry.split(' ')[1]
        self.length = int(entry.split(' ')[3])
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
        
    def print_sequence(self,reverse_strand=False,split_every=None,region=False,rename=False):
        out_sequence = self.sequence if not region else self.sequence[region[0]-1:region[1]]
        if reverse_strand:
            out_sequence = out_sequence.lower().replace('a','0').replace('t','2').replace('c','1').replace('g','3')
            out_sequence = out_sequence.replace('0','T').replace('2','A').replace('1','G').replace('3','C').upper()[::-1]
        out_name = rename if rename else self.name
        split_value = split_every if split_every else len(self.sequence)
        print_info='>{}'.format(out_name)
        print_info = [print_info] + [out_sequence[pos:pos+split_value] for pos in range(0, int(len(out_sequence)), split_value)] 
        return '\n'.join(print_info)
        
    def write(self):
        return self.contig + '\n'

### GFF feature object class
class GFF_feature:
    def __init__(self, entry, contig_number, ID_stat='ID', alt_ID_stat=False):
        self.raw_entry = entry
        self.feature = entry.split('\t')
        self.contig_name = self.feature[0]
        self.seq_type = self.feature[2] ; self.strand = self.feature[6]
        self.frame = int(self.feature[7]) if self.feature[7] in list('012') else 0
        self.start = int(self.feature[3]) ; self.stop = int(self.feature[4])
        self.contig_number = contig_number
        self.coords = '{};{};{}'.format(self.contig_number,self.start,self.stop) 
        self.sequence = []
        self.feature_info = {stat.split('=')[0]:stat.split('=')[1] for stat in self.feature[len(self.feature)-1].split(';')}
        if ID_stat in self.feature_info.keys():
            self.ID = self.feature_info[ID_stat]
        elif alt_ID_stat and alt_ID_stat in self.feature_info.keys():
            self.ID = self.feature_info[alt_ID_stat]
        else:
            print('\nWarning: no {} stat or specified alternative for feature:\n{}\nGenerating ID from co-ordinates...'.format(ID_stat,entry))
            self.ID = self.coords
        self.idx = 0
        self.Parent = self.feature_info.get('Parent') if 'Parent' in self.feature_info else self.ID
        self.family = self.Parent
        self.records = self.feature_info.keys()
        
    def __repr__(self):
        return self.ID

    def add_family(self,heirarchy):
        self.family = heirarchy

    def lookup(self,stat):
        if stat in self.feature_info:
            return self.feature_info[stat]
        elif stat in self.family.all_feature_info:
            output = self.family.all_feature_info[stat]
            if len(output.split(';'))>1:
                return {feature_ID:feature.feature_info.get(stat) for feature_ID,feature in self.family.unique_features.items()}
            else:
                return output
        elif stat in self.family.attributes:
            return self.family.attributes[stat]
        else:
            pattern = re.compile(r"{}".format(stat))
            matching_attributes = [attribute for attribute in self.family.attributes.keys() if pattern.match(attribute)]
            return matching_attributes #if len(matching_attributes) > 0 else 'NA'

    def update(self,stats_to_add,retain_old='raw_',overwrite=False):
        assert type(stats_to_add) == dict, 'stats must be provided as a dictionary'
        for stat in stats_to_add.keys():
            if stat in self.feature_info.keys():
                if not overwrite:
                    stat = retain_old + str(stat)
                    assert stat not in self.feature_info, '{}: cannot have duplicate attributes per feature'.format(stat)                    
            self.feature_info[stat] = stats_to_add.get(stat)
        self.records = self.feature_info.keys()

    def rename_contig(self,new_contig_name):
        assert type(new_contig_name) == str, 'contig names must be strings.'
        self.raw_entry = self.raw_entry.replace(self.contig_name,new_contig_name)
        self.feature = self.raw_entry.split('\t')
        self.contig_name = new_contig_name

    def add_sequence(self,feature_contig,measure_stats=False):
        assert type(feature_contig) == contig and len(feature_contig.sequence) > 0, '{} contig sequence not parsed'.format(feature_contig)
        parse_sequence = feature_contig.sequence[self.start+self.frame-1:self.stop+self.frame].upper()
        if self.strand == '-':
            parse_sequence = parse_sequence.lower().replace('a','0').replace('t','2').replace('c','1').replace('g','3')
            parse_sequence = parse_sequence.replace('0','T').replace('2','A').replace('1','G').replace('3','C').upper()[::-1]
        self.sequence.append(parse_sequence)
        if measure_stats:
            measure_sequence = ''.join(self.sequence[0]) ### [0] see below
            self.sequence_length = len(measure_sequence)
            self.GC = 100*len(measure_sequence.replace('T','').replace('A',''))/self.sequence_length
            self.contig_boundary_dist = min(min(self.start,self.stop),feature_contig.length-max(self.start,self.stop))
        
    def print_sequence(self,fasta_name_stats='ID',split_every=None):
        if type(fasta_name_stats) == str:
            fasta_name_stats = [fasta_name_stats]
        print_info = ['>{}'.format('_'.join([self.feature_info.get(stat) for stat in fasta_name_stats if stat in self.feature_info]))]
        out_sequence = ''.join(self.sequence[0]) ### added 24/02/2023 the [0] to avoid duplicated gene sequences per fasta feature
        split_value = split_every if split_every else len(out_sequence)
        print_info = print_info + [out_sequence[pos:pos+split_value] for pos in range(0, len(out_sequence), split_value)] 
        return '\n'.join(print_info)

    def write(self,update=True):
        if update:
            entry_data_string = ';'.join([str(k) + '=' + str(v) for k,v in self.feature_info.items()])
            return '\t'.join(['\t'.join(self.feature[0:8]), entry_data_string]) + '\n'
        else:
            return self.raw_entry

### GFF feature heirarchy class
class GFF_feature_heirarchy:
    def __init__(self, progenitor):
        self.progenitor = progenitor
        self.count = len(self.progenitor)
        self.all_feature_info = {}
        self.attributes = {}
        self.feature_tally = {}
        self.unique_features = {}
        # create featute type tally
        for feature in self.progenitor:
            if feature.seq_type in self.feature_tally:
                self.feature_tally[feature.seq_type] += 1
            else:
                self.feature_tally[feature.seq_type] = 1
        # split between duplicated and unique feature types per progenitor
        duplicates = []
        for feature in self.progenitor:
            if self.feature_tally.get(feature.seq_type) < 2:
                self.unique_features[feature.seq_type] = feature
            else:
                duplicates.append((feature.ID,feature.seq_type))
        # for the duplicated feature types, add a number to the dictionary lookup
        duplicated_feature_types = [ft for ft in self.feature_tally if self.feature_tally[ft] > 1]
        for ft in duplicated_feature_types:
            duplicated_feature_count = 0
            for ID,feature_type in duplicated_features:
                if feature_type == ft:
                    duplicated_feature_count += 1
                    seq_type_numbered = feature_type + '_{}'.format(str(duplicated_feature_count))
                    self.unique_features[seq_type_numbered] = feature
        # next, iterate through and extract the different feature stats/attributes
        for feature_ID,feature in self.unique_features.items():
            for stat,attribute in feature.feature_info.items():
                if stat not in ['Parent','ID']:
                    if stat in self.all_feature_info:
                        self.all_feature_info.get(stat)[feature_ID]=attribute
                    else:
                        self.all_feature_info[stat] = {feature_ID:attribute}
                    if attribute in self.attributes:
                        self.attributes[attribute] = self.attributes.get(attribute) + ',{}'.format(feature_ID)
                    else:
                        self.attributes[attribute] = feature_ID
        # finally, iterate through the extracted stats/attributes one last time to remove duplication
        self.all_feature_info = { stat: '; '.join(sorted(set(attributes.values()))) for stat,attributes in self.all_feature_info.items()}
        self.attributes = {att: ', '.join(sorted(set(stat.split(',')))) for att,stat in self.attributes.items()}        

### Main GFF file object class
class GFF:
    def __init__(self, file, ID_stat='ID', update_feature_stats=False, alt_ID_stat=None):
        self.file = file
        self.name = self.file.split('/')[-1].split('.gff')[0]
        self.file_info = []
        self.contigs = {}
        self.contig_count = len(self.contigs)
        self.contig_sequence = []
        self.feature_count = 1
        self.features = {} # a dictionary for all features (separate parents and children)
        self.progenitors = {} # a dictionary for parents only to match with children 
        self.children = {} # a dictionary for children only  
        self.renamed_progenitors = {} # a dictionary for matching up old and new names for progenitors, if identified in the data
        self.indexed_features = {} # a dictionary for heirarchies (families) by index
        # take carer!: in order for the indexed_features dictionary to be filled correctly, features must be numbered and appear in numbered order             

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
                    current_feature = GFF_feature(line,current_contig.number,ID_stat,alt_ID_stat) # feature object: current feature
                    # record the stats provided for each feature to generate a list of total available stats
                    if list(current_feature.records) not in all_recorded_stats:
                        all_recorded_stats.append(list(current_feature.records))
                    # update current feature index - this will only work if the features are in order                   
                    if re.search(fr"_0*{str(self.feature_count)}\b", current_feature.ID):
                        current_feature.idx = self.feature_count
                    elif re.search(fr"_0*{str(self.feature_count-1)}\b", current_feature.ID):
                        current_feature.idx = self.feature_count - 1
                    # add new stats from current feature to dictionary
                    if current_feature.ID in self.features:
                        print('Warning: {} file contains multiple entries with ID {}, keeping first.'.format(self.file,current_feature))
                    else:
                        self.features[current_feature.ID] = current_feature
                        if 'Parent' in current_feature.records:
                            self.children[current_feature.ID] = current_feature.Parent
                        else:
                            self.progenitors[current_feature.ID] = [current_feature] # enter data as list
                            self.indexed_features[current_feature.idx] = current_feature.ID
                            # some prokaryote pangenome tools rename GFF features by ID, so this should be handled here
                            if 'prev_ID' in current_feature.records: # 'prev_ID' to record the old ID in PIRATE 'modified_gffs'
                                self.renamed_progenitors[current_feature.lookup('prev_ID')] = current_feature.ID
                            # if the file is in order, progenitors should be just before/after children, so this is the best place to increment the feature count 
                            self.feature_count +=1
                    # add sequence data from the end of the file, if it is there
                elif line.startswith('>'):
                    if len(self.contig_sequence) > 0: # if the sequence of a previous contig is not currently parsed and stored 
                        current_contig.concatenate_sequence(''.join(self.contig_sequence)) # add it to the data for the contig
                        self.contig_sequence = [] # and wipe the storage variable to start parsing the next contig
                    line = line.lstrip('>')
                    for character in range(len(line)):
                        if line in self.contigs or len(line) == 0:
                            if len(line) == 0:
                                print('Warning: No contig name in fasta feature. Contig fasta sequences were incorrectly parsed.')
                            current_contig = self.contigs[line] ; break
                        else:
                            line = line[0:len(line)-1]
                else:
                    self.contig_sequence.append(line)
            current_contig.concatenate_sequence(''.join(self.contig_sequence)) # store the parsed data for the final contig
        
        # summary info
        self.all_recorded_stats = set(sum(all_recorded_stats, []))

        # update offspring list in progenitor dictionary now that all entries have been parsed
        for child_ID,parent_ID in self.children.items(): # add family info for children
            if self.progenitors.get(parent_ID): # this can fail if the feature has been renamed but the parent ID has been left the same
                self.progenitors.get(parent_ID).append(self.features.get(child_ID))
            elif self.renamed_progenitors.get(parent_ID): # if the progenitor ID has been renamed, this needs to be used to stop it from 
                guess_parent_ID = self.renamed_progenitors.get(parent_ID)
                self.progenitors.get(guess_parent_ID).append(self.features.get(child_ID))
            else: # if it does still fail, it might still be possible to identify the parent from the index if the file is in order 
                child_index = self.features.get(child_ID).idx
                guess_parent_ID = self.indexed_features.get(child_index)
                try:
                    self.progenitors.get(guess_parent_ID).append(self.features.get(child_ID))
                except:
                    print('\nError! Parent ID {} not found in file and no Prev_ID or equivalent stat was determined\n'.format(parent_ID))
                    print('Unclear file ordering... cannot guess parent ID\n'.format(guess_parent_ID))
        for progenitor,relatives in self.progenitors.items():
            family_info = GFF_feature_heirarchy(relatives)
            self.progenitors[progenitor] = family_info
            for relative in relatives:
                relative.family = family_info

        # update info (optional?)
        # if contig fasta sequences were provided, add sequences to entries for each locus
        if len(self.contig_sequence) > 0:
            for feature in self.features.keys():
                current_feature = self.features[feature]
                current_contig = self.contigs[current_feature.contig_name]
                current_feature.add_sequence(current_contig, measure_stats=update_feature_stats)
        # and finally, iterate back through the entries to add any holistic stats:
        if update_feature_stats:
            for feature_ID,feature in self.features.items():
                if feature.Parent in self.renamed_progenitors:
                    feature.Parent = self.renamed_progenitors.get(feature.Parent)
                more_info_to_add = {'index': feature.idx,
                                    'all_relative_IDs': feature.family.progenitor,
                                    'sequence_length': feature.sequence_length,
                                    'contig_sequence_length': self.contigs[feature.contig_name].length,
                                    'contig_boundary_distance': feature.contig_boundary_dist,
                                    'GC': feature.GC}
                # Add ID stat to this for entries that were lacking an ID stat 
                if feature_ID == feature.coords:
                    more_info_to_add[ID_stat] = feature.coords
                self.features[feature_ID].update(more_info_to_add)
        # including making a record of the data for all the relatives of each top parent

    def __repr__(self):
        return str(self.name) + '.gff'
    
    def info(self,as_input=False):
        for piece_of_information in self.file_info:
            return piece_of_information if as_input else piece_of_information.lstrip('#')

    def fetch_feature_list(self):
        return [feature.ID for feature in self.features.values()]
    
    def search(self,search_term):
        if search_term in self.features:
            return self.features.get(search_term).coords
        else:
            to_print = []
            for ID,feature in self.progenitors.items():
                for branch in feature.family:
                    if branch.lookup(search_term):
                        to_print.append(branch)
            return(to_print)           

    def feature(self,feature_lookup,ftype=None):
        if feature_lookup in self.features:
            out_feature = self.features.get(feature_lookup)
            if ftype in out_feature.family.unique_features.keys():
                replace_feature = feature_lookup.family.unique_features.get(feature_lookup)
                out_feature = features.get(replace_feature)
        else:
            out_feature = feature_lookup.family.unique_features()
        return out_feature

    def fetch_feature_sequences(self,list_of_loci,list_of_fasta_header_categories,split_every=None):
        out_fasta = []
        if type(list_of_loci) == str:
            list_of_loci = [list_of_loci]
        for locus in list_of_loci:
            locus_feature = self.feature(locus)
            out_fasta.append(locus_feature.print_sequence(list_of_fasta_header_categories,split_every))
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
        for feature_ID,current_feature in self.features.items():
            new_contig_name = new_names[current_feature.contig_name]
            current_feature.rename_contig(new_contig_name)            

    def to_newfile(self,out_file=False,
    rename_contigs=True,update_stats=True,
                   add_FASTA_sequence=True,FASTA_Split_Every=True,
                   skip_entries=[],skip_contigs=[]):
        out_file = out_file if out_file else self.file + '.out'
        if rename_contigs: # the default value (of simply 'True') will number the contigs and name them after the genome
            self.rename_contigs(rename_contigs=rename_contigs) # this can be set to false or 'add_number' to simply add the number
        with open(out_file,'w') as outfile:
            outfile.write(self.info(as_input=True))
            outfile.write('\n')
            for contig in self.contigs.keys():
                if self.contig_count < len(self.contigs.keys()):
                    skip_contigs = skip_contigs + [i for i in range(0,self.contig_count+1)]
                if contig not in skip_contigs:
                    outfile.write(self.contigs[contig].write())
            for feature_ID,current_feature in self.features.items():
                corresponding_contig = current_feature.contig_name
                if corresponding_contig not in skip_contigs and feature_ID not in skip_entries:
                    outfile.write(self.features[feature_ID].write(update=True))
            if add_FASTA_sequence:
                outfile.write('##FASTA\n')
                if FASTA_Split_Every and type(FASTA_Split_Every) != int:
                    FASTA_Split_Every = self.measure_FASTA_sequence_split()
                for contig in self.contigs.keys():
                    if contig not in skip_contigs:
                        outfile.write(self.contigs[contig].print_sequence(split_every=FASTA_Split_Every))
                        outfile.write('\n')