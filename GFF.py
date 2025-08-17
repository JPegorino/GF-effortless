### DEPENDENCIES ###
import re, os

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
        self.coords = '{}~{}~{}'.format(self.contig_number,self.start,self.stop) 
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

    def lookup(self,stat,regex=True):
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
        elif regex:
            pattern = re.compile(r"{}".format(stat))
            matching_attributes = [attribute for attribute in self.family.attributes.keys() if pattern.search(attribute)]
            return matching_attributes #if len(matching_attributes) > 0 else 'NA'
        else:
            return None

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

    def sequence(self,feature_contig,us=0,ds=0,protein=False):
        assert type(feature_contig) == contig and len(feature_contig.sequence) > 0, '{} contig sequence not parsed'.format(feature_contig)
        assert type(us) == int and type(ds) == int, 'upstream and downstream values must be numeric integers'
        raw_ds = ds # to provide the unchanged downstream value for the translation function
        if self.strand == '-':
            us,ds = ds,us # swap upstream and downstream values if the sequnce is on the - strand
        parse_sequence = feature_contig.sequence[self.start+self.frame-1-us:self.stop+self.frame+ds].upper()
        if self.strand == '-':
            parse_sequence = parse_sequence.lower().replace('a','0').replace('t','2').replace('c','1').replace('g','3')
            parse_sequence = parse_sequence.replace('0','T').replace('2','A').replace('1','G').replace('3','C').upper()[::-1]
        return translated(parse_sequence,raw_ds) if protein else parse_sequence
        
    def print_sequence(self,feature_contig,fasta_name_stats='ID',split_every=None,us=0,ds=0,protein=False):
        assert type(feature_contig) == contig and len(feature_contig.sequence) > 0, '{} contig sequence not parsed'.format(feature_contig)
        if type(fasta_name_stats) == str:
            fasta_name_stats = [fasta_name_stats]
        print_info = ['>{}'.format('_'.join([str(self.feature_info.get(stat)).replace(' ','-') for stat in fasta_name_stats if stat in self.feature_info]))]
        parse_sequence = self.sequence(feature_contig,us,ds,protein) ### troubleshooting.
        out_sequence = ''.join(parse_sequence) 
        split_value = split_every if split_every else len(out_sequence)
        print_info = print_info + [out_sequence[pos:pos+split_value] for pos in range(0, len(out_sequence), split_value)] 
        return '\n'.join(print_info)

    def write(self,bed=False,raw=False,bed_score='NA_default',bed_color='NA_default'):
        if raw:
            return self.raw_entry
        elif bed:
            if type(self.family) == GFF_feature_heirarchy:
                progenitor = self.family.progenitor_feature
                feature_start = str(progenitor.start-1)
                feature_end = str(progenitor.stop)
            else:
                feature_start = str(self.start-1)
                feature_end = str(self.stop)
            score_value = self.lookup(bed_score,regex=False)
            score = score_value if type(score_value) in [int,float] else self.idx
            color_value = self.lookup(bed_color,regex=False)
            RGB_format = re.compile('(25[0-5]|2[0-4][0-9]|1?[0-9]{1,2}),(25[0-5]|2[0-4][0-9]|1?[0-9]{1,2}),(25[0-5]|2[0-4][0-9]|1?[0-9]{1,2})')
            color = color_value if type(color_value) == 'str' and RGB_format.match(color_value) else '0,0,0'
            bed_list = [self.contig_name,feature_start,feature_end,self.ID,self.strand,str(score),str(self.start-1),str(self.stop),color]
            return '\t'.join(bed_list)
        else:
            entry_data_string = ';'.join([str(k) + '=' + str(v) for k,v in self.feature_info.items()])
            return '\t'.join(['\t'.join(self.feature[0:8]), entry_data_string]) + '\n'

### GFF feature heirarchy class
class GFF_feature_heirarchy:
    def __init__(self, feature_family):
        self.feature_family = feature_family
        self.progenitor = self.feature_family[0].ID
        self.progenitor_feature = self.feature_family[0]
        self.count = len(self.feature_family)
        self.start = 0
        self.stop = 0
        self.all_feature_info = {}
        self.attributes = {}
        self.feature_tally = {}
        self.unique_features = {}
        # create featute type tally - whilst iterating, add start and end values to lists at the same time
        feature_coordinates = []
        for feature in self.feature_family:
            feature_coordinates.append(feature.start)
            feature_coordinates.append(feature.stop)
            if feature.seq_type in self.feature_tally:
                self.feature_tally[feature.seq_type] += 1
            else:
                self.feature_tally[feature.seq_type] = 1
        # split between duplicated and unique feature types per feature_family
        duplicated_features = []
        for feature in self.feature_family:
            if self.feature_tally.get(feature.seq_type) < 2:
                self.unique_features[feature.seq_type] = feature
            else:
                duplicated_features.append((feature.ID,feature.seq_type))
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
        # calculate the max start and stop, so they can be called directly from the GFF_heirarchy instance
        self.start = min(feature_coordinates)
        self.stop = max(feature_coordinates)
        # NB: these do not have directionality ('start' means 'first bp', not 'CDS start' or 'gene start')
        # finally, iterate through the extracted stats/attributes one last time to remove duplication
        self.all_feature_info = { stat: '; '.join(sorted(set(attributes.values()))) for stat,attributes in self.all_feature_info.items()}
        self.attributes = {att: ', '.join(sorted(set(stat.split(',')))) for att,stat in self.attributes.items()}        

    def __str__(self):
        return self.progenitor # return just the progenitor ID if the heirarchy is printed as string

### Main GFF file object class
class GFF:
    def __init__(self, file, ID_stat='ID', update_feature_stats=False, alt_fasta_file=None, alt_ID_stat=None):
        self.file = file
        self.name = self.file.split('/')[-1].split('.gff')[0]
        self.file_info = []
        self.in_order = True
        self.includes_FASTA = False
        self.contigs = {}
        self.contig_count = len(self.contigs)
        self.contig_sequence = []
        self.feature_count = 1
        self.features = {} # a dictionary for all features (separate parents and children) to match feature objects
        self.families = {} # a dictionary for parents only to match with children 
        self.children = {} # a dictionary of children (only) to match with IDs of their parent 
        self.renamed_parents = {} # a dictionary for parents with modified IDs to match with their old IDs 
        self.indexed_features = {} # a dictionary for heirarchies (families) by index
        self.coords = {} # a dictionary for all features families by their co-ordinates

        # read through the file and parse information line by line
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
                        # self.coords[current_contig.number] = {i:"None" for i in range(1,current_contig.length)}
                        # self.coords.get(current_contig.number)[0] = self.coords.get(current_contig.number)[current_contig.length+1] = "contig_break"
                    else:
                        self.file_info.append(line)
                elif line.split('\t')[0] in self.contigs:
                    current_contig = self.contigs[line.split('\t')[0]] # contig object: current contig
                    current_feature = GFF_feature(line,current_contig.number,ID_stat,alt_ID_stat) # feature object: current feature
                    # record the stats provided for each feature to generate a list of total available stats
                    if list(current_feature.records) not in all_recorded_stats:
                        all_recorded_stats.append(list(current_feature.records))
                    # add new stats from current feature to dictionary
                    if current_feature.ID in self.features:
                        print('Warning: {} file contains multiple entries with ID {}, keeping first.'.format(self.file,current_feature))
                    else:
                        self.features[current_feature.ID] = current_feature
                        if 'Parent' in current_feature.records:
                            self.children[current_feature.ID] = current_feature.Parent
                        else: # if there is no parent, we define this feature as a progenitor 
                            self.families[current_feature.ID] = [current_feature] # enter data as list, starting with progenitor
                            # some prokaryote pangenome tools rename GFF features by ID, which should be handled here
                            # e.g. PIRATE modifies parent feature IDs without updating the 'Parent' stat in the child entries
                            if 'prev_ID' in current_feature.records: # 'prev_ID' used by PIRATE 'modified_gffs' to record the old ID
                                self.renamed_parents[current_feature.lookup('prev_ID')] = current_feature.ID
                    # add sequence data from the end of the file, if it is there
                elif line.startswith('>'):
                    self.includes_FASTA = True
                    if len(self.contig_sequence) > 0: # if the sequence of a previous contig is not currently parsed and stored 
                        current_contig.concatenate_sequence(''.join(self.contig_sequence)) # add it to the data for the contig
                        self.contig_sequence = [] # and wipe the storage variable to start parsing the next contig
                    line = line.lstrip('>')
                    for character in range(len(line)):
                        if line in self.contigs or len(line) == 0:
                            if len(line) == 0:
                                print('Warning: Contig names are not matched in FASTA sequence headers.')
                            current_contig = self.contigs[line] ; break
                        else:
                            line = line[0:len(line)-1]
                else:
                    self.contig_sequence.append(line)
            current_contig.concatenate_sequence(''.join(self.contig_sequence)) # store the parsed data for the final contig
        
        # if no sequence data was found in the GFF file, try to find FASTA sequences in a separate file
        if not self.includes_FASTA: # if no contig FASTA sequences were identified in the GFF file
            if not alt_fasta_file: # use any user specified alternative FASTA
                for extension in ['.fna','.fasta','.fa']: # or seatch for an epynomous file (in the same directory) with a FASTA extension
                    alt_fasta_file = self.file.replace('.gff',extension)
                    if os.file.exists(alt_fasta_file):
                        continue
            assert os.file.exists(alt_fasta_file), 'No contigs in GFF file and no alternative FASTA file {} found.'.format(alt_fasta_file)
            with open(alt_fasta_file,'r') as infile:
                for line in infile:
                    if line.startswith('>'):
                        if len(self.contig_sequence) > 0: # if the sequence of a previous contig is not currently parsed and stored 
                            current_contig.concatenate_sequence(''.join(self.contig_sequence)) # add it to the data for the contig
                            self.contig_sequence = [] # and wipe the storage variable to start parsing the next contig
                        line = line.lstrip('>')
                        for character in range(len(line)):
                            if line in self.contigs or len(line) == 0:
                                if len(line) == 0:
                                    print('Warning: Contig names are not matched in FASTA sequence headers.')
                                current_contig = self.contigs[line] ; break
                            else:
                                line = line[0:len(line)-1]
                    else:
                        self.contig_sequence.append(line)
                current_contig.concatenate_sequence(''.join(self.contig_sequence)) # store the parsed data for the final contig
        
        # generate summary info
        self.all_recorded_stats = set(sum(all_recorded_stats, []))

        # re-iterate through the data to add information that requires the entire parsed file.
        # firstly, add offspring to families dictionary (currently just includes progenitors).
        ### CURRENTLY, THIS CODE IS NOT SET UP TO HANDLE INDIRECT PARENTS/CHAINS --- TO FIX ###
        for child_ID,parent_ID in self.children.items():
            if self.families.get(parent_ID):
                self.families.get(parent_ID).append(self.features.get(child_ID))
            # handle children whose parents have a modified ID stat (PIRATE) using a renaming dictionary.
            elif self.renamed_parents.get(parent_ID):
                guess_parent_ID = self.renamed_parents.get(parent_ID)
                self.families.get(guess_parent_ID).append(self.features.get(child_ID))
            else: 
                print('\nError! Parent ID {} not found in file and no "Prev_ID" or equivalent stat for conversion\n'.format(parent_ID))
        # use the completed dictionary to create and add heirarchy objects to the GFF/feature objects and use these to add feature indices
        last_contig,current_index=(0,0) # to keep track of index order (helps check order of features in file)
        for progenitor,family_list in self.families.items():
            family_heirarchy = GFF_feature_heirarchy(family_list)
            self.families[progenitor] = family_heirarchy
            # determine feature indices for progenitors (features in the same family should have matching indices)
            current_contig_number = family_heirarchy.progenitor_feature.contig_number
            if current_contig_number != last_contig:
                # using round prevents floating point errors affecting the dictionary keys
                self.coords[round(current_index+0.000001,6)] = "contig break"
                last_contig,last_feature,current_index = (current_contig_number,0,current_contig_number)
                self.coords[current_index] = "contig break"
            if family_heirarchy.stop > last_feature:
                current_index = round(current_index+0.000001,6)
                last_feature = family_heirarchy.stop
                self.coords[current_index] = family_heirarchy
            elif family_heirarchy.stop < last_feature:
                print('\Warning! Feature order in file does not match ordering on contig. Indices will be invalid.\n'.format(parent_ID))                
            # add heirarchy objects and indices to all features
            for relative in family_list:
                relative.family = family_heirarchy
                relative.idx = current_index

        # and optionally, iterate back through all the features to add any holistic stats:
        if update_feature_stats:
            for feature_ID,feature in self.features.items():
                # If the feature 'Parent' stat had been renamed (i.e. by PIRATE), update/restore this first
                if feature.Parent in self.renamed_parents:
                    feature.Parent = self.renamed_parents.get(feature.Parent)
                # Then, create a dictionary of non-existing features to add 
                more_info_to_add = {'index': feature.idx,
                                    'family': feature.family.progenitor}
                # Include the ID stat for any entries that were missing one 
                if feature_ID == feature.coords:
                    more_info_to_add[ID_stat] = feature.coords
                # If contig fasta sequences were provided, add sequence stats
                if len(self.contig_sequence) > 0:
                    feature_contig = self.contigs[feature.contig_name]
                    feature_sequence = ''.join(feature.sequence(feature_contig))
                    more_info_to_add['sequence_length'] = len(feature_sequence)
                    more_info_to_add['contig_sequence_length'] = feature_contig.length
                    more_info_to_add['contig_boundary_distance'] = min(min(feature.start,feature.stop),feature_contig.length-max(feature.start,feature.stop))
                    more_info_to_add['GC'] = 100*len(feature_sequence.replace('T','').replace('A',''))/len(feature_sequence)
                # Finally, update the feature stats from the dictionary
                self.features[feature_ID].update(more_info_to_add)

    def __repr__(self):
        return str(self.name) + '.gff'
    
    def info(self,as_input=False):
        for piece_of_information in self.file_info:
            return piece_of_information if as_input else piece_of_information.lstrip('#')

    def fetch_feature_list(self):
        return [feature.ID for feature in self.features.values()]
    
    def feature(self,feature_lookup,feature_type=None,strictly_first=True,regex=False):
        if feature_lookup in self.features:
            out_feature = self.features.get(feature_lookup)
             # the code below allows retreival of data for an alternate feature type in the same family (e.g. CDS info for a gene) if there is only one
            if feature_type in out_feature.family.unique_features.keys(): # this only works for 1:1 rep=lationships, e.g. one CDS per gene
                out_feature = out_feature.family.unique_features.get(feature_type)
            elif feature_type:
                print('Warining: no unique {} feature in {} feature heirarchy'.format(feature_type,feature_lookup))
        else:
             # the code below will attempt retrieval of one or more features with data that matches the lookup, if it is not an ID
            families_with_match = []
            for feature_family in self.families.values():
                progenitor_feature = self.features.get(feature_family.progenitor)
                if progenitor_feature.lookup(feature_lookup,regex=regex):
                    families_with_match.append(progenitor_feature)
            if len(families_with_match) < 1:
                out_feature = None
            elif strictly_first or len(families_with_match) == 1:
                out_feature = families_with_match[0]
                if len(families_with_match) > 1:
                    print('Warining: {} feature families with data to match {}, returning first'.format(len(families_with_match),feature_lookup))
                out_feature = families_with_match
            else:
                out_feature = families_with_match
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
                    outfile.write(self.features[feature_ID].write())
            if add_FASTA_sequence:
                outfile.write('##FASTA\n')
                if FASTA_Split_Every and type(FASTA_Split_Every) != int:
                    FASTA_Split_Every = self.measure_FASTA_sequence_split()
                for contig in self.contigs.keys():
                    if contig not in skip_contigs:
                        outfile.write(self.contigs[contig].print_sequence(split_every=FASTA_Split_Every))
                        outfile.write('\n')

### HELPER FUNCTIONS ###

def translated(nucleotide_string,downstream=0):
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
    protein_sequence = []
    # If the length of downstream sequence included is not a multiple of 3, shorten the sequence first until it is
    nucleotide_string = nucleotide_string[0:len(nucleotide_string)-(downstream % 3)]
    print(downstream)
    # Iterate over the sequence in steps of 3 nucleotides (codon)
    frame_shift = len(nucleotide_string) % 3 # removes nucleotides from start of sequence to keep to frame
    for nt_pos in range(0 + frame_shift, len(nucleotide_string) - 2, 3):
        codon = nucleotide_string[nt_pos:nt_pos+3]
        amino_acid = codon_table.get(codon, '?')  # uses '?' to represent unknown codons
        protein_sequence.append(amino_acid)
    return ''.join(protein_sequence)
