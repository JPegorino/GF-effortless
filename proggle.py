### IMPORT DEPENDENCIES ###
import sys, os, io, re, csv, warnings, traceback, argparse

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
    def __init__(self, entry, contig_number, ID_stat='ID', alt_ID_stat='locus_tag', input_index=1, verbose=False):
        self.raw_entry = entry
        self.feature = entry.split('\t')
        self.contig_name = self.feature[0]
        self.seq_type = self.feature[2] ; self.strand = self.feature[6]
        self.frame = int(self.feature[7]) if self.feature[7] in list('012') else 0
        self.start = int(self.feature[3]) ; self.stop = int(self.feature[4])
        self.contig_number = contig_number
        self.coords = '{}~{}~{}'.format(self.contig_number,self.start,self.stop)
        self.feature_info = {stat.split('=')[0]:stat.split('=')[1] for stat in self.feature[len(self.feature)-1].split(';')}
        self.records = self.feature_info.keys()
        if ID_stat in self.records:
            self.ID = self.feature_info[ID_stat]
        elif alt_ID_stat and alt_ID_stat in self.feature_info.keys():
            self.ID = self.feature_info[alt_ID_stat]
        else:
            if verbose:
                warnings.warn(f'No {ID_stat} stat or alternative for {self.seq_type} feature at loci {self.coords}. Adding ID...')
            try:
                self.ID = '_'.join([str(self.contig_name),self.seq_type,str(input_index).zfill(5)])
                ID_added = {'ID': self.ID}
                ID_added.update(self.feature_info)
                self.feature_info = ID_added
                first_stat = list(self.records)[0]
                self.raw_entry = self.raw_entry.replace(first_stat,f'ID={self.ID};{first_stat}')
            except:
                raise ValueError(f'Could not determine or generate {ID_stat} stat or alternative for {self.seq_type} feature at loci {self.coords}.')
        self.idx = 0
        self.Parent = self.feature_info.get('Parent') if 'Parent' in self.feature_info else self.ID
        self.family = self.Parent
    
    def __repr__(self):
        return self.ID

    def add_family(self,heirarchy):
        self.family = heirarchy
    
    def progenitor(self,feature_dictionary,depth=0):
        if depth > 10: # recursion limit - avoid cyclic / infinite recursion 
            raise RuntimeError(f"{self.ID} feature progenitor detection exceeded recursion limit. Cannot handle feature families with > 10 levels.")
        if self.ID == self.Parent: 
            return self.ID
        else:
            try:
                parent_feature = feature_dictionary.get(self.Parent)
                return parent_feature.progenitor(feature_dictionary, depth + 1)
            except:
                raise KeyError(f"{self.ID} feature progenitor detection failed. Incorrect 'Parent' stat {parent_feature.Parent} for {parent_feature.ID} feature.")

    def lookup(self,stat,regex=True,default_value=None):
        if stat in self.feature_info:
            return self.feature_info[stat]
        elif stat in self.family.related_feature_info:
            output = self.family.related_feature_info[stat]
            if len(output.split(';'))>1:
                # handles the issue of features in same family with conflicting data for same stat
                return {feature_ID:feature.feature_info.get(stat) for feature_ID,feature in self.family.unique_features.items()}
            else:
                return output
        elif stat in self.family.attributes:
            return self.family.attributes[stat]
        elif regex:
            pattern = re.compile(r"{}".format(stat))
            matching_attributes = [attribute for attribute in self.family.attributes.keys() if pattern.search(attribute)]
            return matching_attributes #if len(matching_attributes) > 0 else default_value
        else:
            return default_value

    def update(self,stats_to_add,retain_old='raw_',overwrite=False):
        if not type(stats_to_add) == dict:
            raise TypeError('stats must be provided as a dictionary')
        for stat in stats_to_add.keys():
            if stat in self.feature_info.keys():
                if not overwrite:
                    stat = retain_old + str(stat)
                    if not stat not in self.feature_info:
                        raise AssertionError('{}: cannot have duplicate attributes per feature'.format(stat))
            self.feature_info[stat] = stats_to_add.get(stat)
        self.records = self.feature_info.keys()

    def rename_contig(self,new_contig_name):
        if not type(new_contig_name) == str:
            raise TypeError('contig names must be strings.')
        self.raw_entry = self.raw_entry.replace(self.contig_name,new_contig_name)
        self.feature = self.raw_entry.split('\t')
        self.contig_name = new_contig_name

    def sequence(self,feature_contig,us=0,ds=0,protein=False):
        if not type(feature_contig) == contig and len(feature_contig.sequence) > 0:
            raise Exception('{} contig sequence not parsed'.format(feature_contig))
        if not type(us) == int and type(ds) == int:
            raise TypeError('upstream and downstream values must be numeric integers')
        if protein:
            us *= 3 ; ds *= 3
        if self.strand == '-':
            us,ds = ds,us # swap upstream and downstream values if the sequnce is on the - strand
        parse_sequence = feature_contig.sequence[self.start+self.frame-1-us:self.stop+self.frame+ds].upper()
        if self.strand == '-':
            parse_sequence = parse_sequence.lower().replace('a','0').replace('t','2').replace('c','1').replace('g','3')
            parse_sequence = parse_sequence.replace('0','T').replace('2','A').replace('1','G').replace('3','C').upper()[::-1]
        return translated(parse_sequence) if protein else parse_sequence

    def print_sequence(self,feature_contig,fasta_name_stats='ID',split_every=None,us=0,ds=0,protein=False):
        if not type(feature_contig) == contig and len(feature_contig.sequence) > 0:
            raise Exception('{} contig sequence not parsed'.format(feature_contig))
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
                progenitor = self.family.progenitor
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
        self.progenitor = self.feature_family[0]
        self.count = len(self.feature_family)
        self.start = 0
        self.stop = 0
        self.related_feature_info = {}
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
                    if stat in self.related_feature_info:
                        self.related_feature_info.get(stat)[feature_ID]=attribute
                    else:
                        self.related_feature_info[stat] = {feature_ID:attribute}
                    if attribute in self.attributes:
                        self.attributes[attribute] = self.attributes.get(attribute) + ',{}'.format(feature_ID)
                    else:
                        self.attributes[attribute] = feature_ID
        # calculate the max start and stop, so they can be called directly from the GFF_heirarchy instance
        self.start = min(feature_coordinates)
        self.stop = max(feature_coordinates)
        # NB: these do not have directionality ('start' means 'first bp', not 'CDS start' or 'gene start')
        # finally, iterate through the extracted stats/attributes one last time to remove duplication
        self.related_feature_info = { stat: '; '.join(sorted(set(attributes.values()))) for stat,attributes in self.related_feature_info.items()}
        self.attributes = {att: ', '.join(sorted(set(stat.split(',')))) for att,stat in self.attributes.items()}

    def __str__(self):
        return self.progenitor.ID # return just the progenitor ID if the heirarchy is printed as string

    def locus_tag(self,as_feature=True): # determine a locus tag 
        if 'locus_tag' in self.attributes:
            locus_tag = self.attributes.get('locus_tag')
        else:
            candidates = [f.ID for f in self.feature_family]
            shortest = min(candidates, key=len)
            if all(shortest in feature_ID for feature_ID in candidates):
                candidates = [f for f in candidates if shortest not in f]
                locus_tag = shortest
            else:
                locus_tag = self.progenitor.ID
        return self.feature_family.get(locus_tag) if as_feature else locus_tag

### Main GFF file object class
class GFF:
    def __init__(self, file, ID_stat='ID', update_feature_stats=False, alt_fasta_file=None, alt_ID_stat=None, verbose=False):
        self.file = file
        self.name = os.path.basename(os.path.splitext(self.file)[0])
        self.metadata = {'Genome': self.name}
        self.file_info = []
        self.in_order = True
        self.includes_FASTA = False
        self.feature_type_dict = {}
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
            duplicate_count = 0 # a record of how many troublesome features there are
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
                    current_feature = GFF_feature(line,current_contig.number,ID_stat,alt_ID_stat,verbose=verbose) # feature object: current feature
                    # record the stats provided for each feature to generate a list of total available stats
                    if list(current_feature.records) not in all_recorded_stats:
                        all_recorded_stats.append(list(current_feature.records))
                    if current_feature.seq_type not in self.feature_type_dict:
                        self.feature_type_dict[current_feature.seq_type] = 1
                    else:
                        self.feature_type_dict[current_feature.seq_type] += 1
                    # add new stats from current feature to dictionary
                    if current_feature.ID in self.features:
                        duplicate_count +=1
                        # 
                        alternate_feature = GFF_feature(line,current_contig.number,ID_stat,alt_ID_stat,verbose=verbose,input_index=duplicate_count)
                        if current_feature.ID == alternate_feature.ID: # handle exceptions where IDs are genuine duplicates
                            warnings.warn('{} file contains multiple entries with ID {}, keeping first.'.format(self.file,current_feature))
                            duplicate_count -=1
                            continue
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
                                warnings.warn('Contig names are not matched in FASTA sequence headers.')
                            current_contig = self.contigs[line] ; break
                        else:
                            line = line[0:len(line)-1]
                else:
                    self.contig_sequence.append(line)
            current_contig.concatenate_sequence(''.join(self.contig_sequence)) # store the parsed data for the final contig

        # if no sequence data was found in the GFF file, try to find FASTA sequences in a separate file
        if not self.includes_FASTA: # if no contig FASTA sequences were identified in the GFF file
            if not alt_fasta_file: # use any user specified alternative FASTA
                for fasta_extension in ['.fna','.fasta','.fa']: # or search for an epynomous file (in the same directory) with a FASTA extension
                    alt_fasta_file = self.file.replace('.gff',fasta_extension)
                    if os.path.exists(alt_fasta_file):
                        break
            validate_input_file(alt_fasta_file, enforced_extension=False, more_context='No sequences in GFF file and no alternative FASTA path provided.\nInferred alternative ')
            with open(alt_fasta_file,'r') as infile:
                for line in infile:
                    line = line.rstrip('\n')
                    if line.startswith('>'):
                        if len(self.contig_sequence) > 0: # if the sequence of a previous contig is not currently parsed and stored
                            current_contig.concatenate_sequence(''.join(self.contig_sequence)) # add it to the data for the contig
                            self.contig_sequence = [] # and wipe the storage variable to start parsing the next contig
                        line = line.lstrip('>')
                        for character in range(len(line)):
                            if line in self.contigs or len(line) == 0:
                                if len(line) == 0:
                                    warnings.warn('Contig names are not matched in FASTA sequence headers.')
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
        ### CURRENTLY, THIS CODE IS NOT SET UP TO HANDLE PIRATE MODIFIED GFFS WITH INDIRECT PARENTS/CHAINS ###
        for child_ID,parent_ID in self.children.items():
            if self.families.get(parent_ID):
                self.families.get(parent_ID).append(self.features.get(child_ID))
            # handle children whose parents have a modified ID stat (PIRATE) using a renaming dictionary.
            elif self.renamed_parents.get(parent_ID):
                guess_parent_ID = self.renamed_parents.get(parent_ID)
                self.families.get(guess_parent_ID).append(self.features.get(child_ID))
            elif parent_ID in self.children: # if the parent is not the root parent
                progenitor_ID = self.features.get(parent_ID).progenitor(self.features)
                self.families.get(progenitor_ID).append(self.features.get(child_ID))
            else:
                raise KeyError(f'Could not trace feature heirarchy of {child_ID}. Parent ID {parent_ID} may be invalid.\nNo alternative "Prev_ID" or equivalent stat for conversion.')
        # use the completed dictionary to create and add heirarchy objects to the GFF/feature objects and use these to add feature indices
        last_contig,current_index=(0,1000000-1) # to keep track of index order (helps check order of features in file)
        for progenitor,family_list in self.families.items():
            family_heirarchy = GFF_feature_heirarchy(family_list)
            self.families[progenitor] = family_heirarchy
            # determine feature indices for progenitors (features in the same family should have matching indices)
            current_contig_number = family_heirarchy.progenitor.contig_number
            if current_contig_number != last_contig: # start of new contig
                last_contig,last_feature,current_index = (current_contig_number,0,1000000*current_contig_number)
            if family_heirarchy.stop > last_feature:
                current_index +=1
                last_feature = family_heirarchy.stop
                self.indexed_features[current_index] = family_heirarchy
            elif family_heirarchy.stop < last_feature:
                warnings.warn('Feature order in file does not match ordering on contig. Indices will be invalid.')
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
                                    'family': feature.family.progenitor.ID}
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

    def rename_contigs(self,rename_contigs=True):
        rename_contigs = rename_contigs if type(rename_contigs) == str else self.name
        new_names = {} ; new_contigs = {} # initialise blank dictionaries
        for contig_ID,current_contig in self.contigs.items():
            if rename_contigs in ['n','number','#','Number','idx']:
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

    def feature(self,feature_lookup,as_family=False,feature_type=None,default_value=None):
        if feature_lookup in self.features:
            out_feature = self.features.get(feature_lookup)
        elif feature_lookup in self.indexed_features:
                out_feature = self.indexed_features.get(feature_lookup)
                try:
                    out_feature = out_feature.progenitor
                except: 
                    warnings.warn(f'No progenitor feature for {feature_lookup} in feature family\
                     {self.indexed_features.get(feature_lookup)}, using default: {default_value}.')
                    out_feature = default_value
        else:
            out_feature = default_value
        if isinstance(out_feature, GFF_feature):
            # the code below allows retreival of data for an alternate feature type in the same family (e.g. CDS info for a gene) if there is only one
            if feature_type in out_feature.family.unique_features.keys(): # this only works for 1:1 relationships, e.g. one CDS per gene
                out_feature = out_feature.family.unique_features.get(feature_type)
            elif feature_type:
                warnings.warn('No unique {} feature in {} feature heirarchy'.format(feature_type,feature_lookup))
        if as_family:
                out_feature = out_feature.family
        return out_feature # returns exactly one feature object, or 'None'

    def search_features(self,feature_lookup,just_feature_type=None,regex=False):
     # returns a list of features that match a search query (lookup)
        features_with_match = []
        if feature_lookup in self.features:
            out_feature = self.feature(feature_lookup,feature_type=just_feature_type)
            features_with_match.append(out_feature)
        elif feature_lookup in self.indexed_features:
            out_feature_family = self.feature(feature_lookup)
            if just_feature_type:
                for feature_type in out_feature_family.unique_features.keys():
                    if feature_type.startswith(just_feature_type):
                        features_with_match.append(out_feature)
            else:
                for feature in out_feature_family.family:
                    features_with_match.append(out_feature)
        else:
            for feature_family in self.families.values():
                progenitor_feature = self.features.get(feature_family.progenitor.ID)
                if progenitor_feature.lookup(feature_lookup,regex=regex):
                    features_with_match.append(progenitor_feature)
            features_with_match = features_with_match if len(features_with_match) > 0 else None
        return features_with_match # returns a list of feature objects, or 'None'

    def feature_region(self,feature_lookup,feature_type='CDS',bf=0,af=0,in_bp=False,strand_aware=False,include_self=True):
        out_features = []
        current_feature = self.feature(feature_lookup)
        if not current_feature:
            raise ValueError(f'input {feature_lookup} must be ID or proggle-formatted index of feature in file.')
        current_feature_index = current_feature.idx
        if strand_aware and current_feature.strand == '-':
            bf,af = af,bf # swap before and after distances if the sequence is on the - strand
        if in_bp: # bf and af are distances (bp) befopre/after gene (instead of feature counts)
            current_contig_length = self.contigs.get(current_feature.contig_name).length
            bf, af = max(current_feature.start - bf, 1), min(current_feature.stop + af, current_contig_length)
            current_index, current_start = current_feature.idx, current_feature.start
            while current_start > bf:
                current_index -= 1
                bf_feature = self.feature(current_index,feature_type=feature_type)
                if not bf_feature:
                    break
                current_start = bf_feature.start
            bf = current_feature_index - current_index
            current_index, current_stop = current_feature.idx, current_feature.stop
            while current_stop < af:
                current_index += 1
                af_feature = self.feature(current_index,feature_type=feature_type)
                if not af_feature:
                    break
                current_stop = af_feature.stop
            af = current_index - current_feature_index
        bf = current_feature_index - bf
        af = current_feature_index + af
        for i in range(bf,current_feature_index):
            if i in self.indexed_features:
                bf_feature = self.feature(i,feature_type=feature_type)
                out_features.append(bf_feature)
        if include_self:
            out_features.append(current_feature)
        for i in range(current_feature_index+1,af+1):
            if i in self.indexed_features:
                af_feature = self.feature(i,feature_type=feature_type)
                out_features.append(af_feature)
        return out_features # returns a list of feature objects, or 'None'

    def fetch_feature_list(self):
        return [feature.ID for feature in self.features.values()]

    def fetch_feature_sequences(self,list_of_loci,list_of_fasta_header_categories,split_every=None):
        out_fasta = []
        if type(list_of_loci) == str:
            list_of_loci = [list_of_loci]
        for locus in list_of_loci:
            locus_feature = self.feature(locus)
            out_fasta.append(locus_feature.print_sequence(list_of_fasta_header_categories,split_every))
        return '\n'.join(out_fasta)

    def add_feature_data(self,analysis_file,analysis_type='general',pad_missing=True,header_line_expected=True,only_extract_columns=False):
        input_delimiter = ',' if analysis_file.endswith('.csv') else '\t' # always assumes tab-delimited input unless explicit csv file
        analysis_types = { # a dictionary to match each analysis tool with the locus ID column and default keep columns
            'pangenome' : (None, [0], 1), # currently, only pangenome tables are tested and supported
            'general' : (0, None, header_line_expected)
        } # a dictionary of sensible default values for the output tables of different tools
        feature_column, keep_columns, header_line_expected = analysis_types.get(analysis_type)
        header_line = validate_input_file(analysis_file,extract_header=True,delimiter=input_delimiter) 
        # replace None-values with sensible alternatives
        if not header_line_expected: # custom headers (validate_input_file just selects line 1 of file)
            header_line = [ analysis_type + '_' + str(i) for i in range(1,len(header_line)+1) ]
        if not feature_column: # genome (GFF file) name (assumes pangenome table)
            estimate_feature_column_1 = estimate_feature_column_2 = self.name
            for special_char in list('-.?!()[]{}|,'):
                estimate_feature_column_2 = estimate_feature_column_2.replace(special_char,'_')
            if estimate_feature_column_1 in header_line:
                feature_column = header_line.index(estimate_feature_column_1)
            elif estimate_feature_column_2 in header_line:
                feature_column = header_line.index(estimate_feature_column_2)
            else:
                raise ValueError(f"Cannot add {analysis_type} data - feature column {feature_column} not identified in \n{header_line}.")
                if analysis_type == 'pangenome':
                    raise ValueError(f"Estimated alternatives {estimate_feature_column_1} and {estimate_feature_column_2} not identified either.")

        if not keep_columns: # all columns
            keep_columns = [ i for i in range(1,len(header_line)+1) ]
        # if there are user-specified columns to keep, override the defaults with these 
        if only_extract_columns:
            only_extract_columns, parse_column_numbers = (only_extract_columns.split(','), [])
            for i in only_extract_columns:
                if '-' in i:
                    try:
                        i_start, i_end = map(int, i.split("-"))
                        j = list(range(i_start, i_end+1))
                    except:
                        raise ValueError("only_extract_columns must be integers or a numeric range.")
                else:
                    j = [int(i)]
                parse_column_numbers = parse_column_numbers + j
                only_extract_columns = sorted(parse_column_numbers)
            try: # if all the values are numeric integers, convert them to column numbers
                keep_columns = [int(i)-1 for i in only_extract_columns]
            except:
                keep_columns = [header_line.index(header)-1 for header in only_extract_columns]
        # extract the data from the file
        add_data = {}
        with open(analysis_file, 'r', newline='', encoding='utf-8-sig') as add_file:
            analysis_data = csv.reader(add_file, delimiter=input_delimiter) # quoting=csv.QUOTE_NONE
            for row in analysis_data:
                key_data = row[feature_column]
                if ';' in key_data or ':' in key_data:
                    key_list = row[feature_column].replace('(','').replace('(','').split(':')
                    for keys in key_list:
                        for key in keys.split(';'):
                            add_data[key] = {header_line[i].replace(' ','_'):row[i] for i in keep_columns}
                else:
                    add_data[key_data] = {header_line[i].replace(' ','_'):row[i] for i in keep_columns}
        for feature_name,feature in self.features.items():
            if feature_name in add_data.keys():
                feature_analysis_data = add_data.get(feature_name)
                feature.update(feature_analysis_data)
            else:
                if pad_missing:
                    filler_dict = {header_line[i].replace(' ','_'):'None' for i in keep_columns}
                    feature.update(filler_dict)

    def add_metadata(self,metadata_filepath,only_extract_columns=False,input_delimiter='\t'):
        header_line = validate_input_file(metadata_file,extract_header=True,delimiter=input_delimiter)
        input_delimiter = ',' if metadata_file.endswith('.csv') else '\t' # always assumes tab-delimited input unless explicit csv file
        # if there are user-specified columns to keep, extract indices for these 
        if only_extract_columns:  
            only_extract_columns = only_extract_columns.split(';')
            try: # if all the values are numeric integers, convert them to column numbers
                keep_columns = [int(i)-1 for i in only_extract_columns]
            except:
                keep_columns = [header_line.index(header)-1 for header in only_extract_columns]
        else:
            keep_columns = [ i for i in range(1,len(header_line)+1) ] # all columns
        # extract the data from the file
        with open(metadata_filepath, 'r', newline='', encoding='utf-8-sig') as metadata_file:
            metadata = csv.reader(metadata_file, delimiter=input_delimiter)
            for row in metadata:
                if row[0] == self.name:
                    more_metadata = {header_line[i]:row[i] for i in keep_columns}
                    break
        self.metadata = self.metadata + more_metadata

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

### HELPER FUNCTIONS and objects ###

 # a function to de-bug input files and extract header lines (if present)
def validate_input_file(input_path,enforced_extension=None,more_context='',extract_header=False,delimiter='\t'):
    if not os.path.exists(input_path):
        raise FileNotFoundError(f'{more_context}{input_path} is not an accurate path to an existing file.')
    if os.path.getsize(input_path) == 0:
        raise EOFError(f'{input_path} is an empty file.')
    if enforced_extension:
        input_file_extension = os.path.splitext(input_path)[1] # extract the file extension
        if input_file_extension != enforced_extension:
            raise ValueError(f'Must have {enforced_extension} extension, not {input_file_extension}.')
    if extract_header:
        with open(input_path, newline='', encoding='utf-8-sig') as infile:
            return next(csv.reader(infile, delimiter=delimiter), None)

 # a function to translate nucleotide sequences to protein
def translated(nucleotide_string):
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
    if not len(nucleotide_string) % 3 == 0:
        raise ValueError("Nucleotide sequence length must be a multiple of three.")
    for nt_pos in range(0, len(nucleotide_string) - 2, 3):
        codon = nucleotide_string[nt_pos:nt_pos+3]
        amino_acid = codon_table.get(codon, '?')  # uses '?' to represent unknown codons
        protein_sequence.append(amino_acid)
    return ''.join(protein_sequence)

# Main PROGGLE script:
if __name__ == "__main__":

       # Custom error and warning format
    def ProggleWarningFormat(message, category, filename, lineno, file=None, line=None):
        return f"{category.__name__}: {message}\n{filename} line {lineno}\n"

    def ProggleErrorFormat(exc_type, exc_value, exc_traceback):
        from_traceback = traceback.extract_tb(exc_traceback)[-1]
        filename,lineno = (from_traceback.filename,from_traceback.lineno)
        print(f"{exc_type.__name__}: {exc_value}\n{filename} line {lineno}")
    
    warnings.formatwarning = ProggleWarningFormat
    sys.excepthook = ProggleErrorFormat

       # Custom help page format
    class ProggleHelpFormatter(argparse.RawTextHelpFormatter):
        def _split_lines(self, text, width):
            return text.splitlines()
    
      # PROGGLE main script functions
    def parse_args():
        parser = argparse.ArgumentParser(
            description="Quickly edit GFF files, search features, and add or extract data in various formats.",
            formatter_class=ProggleHelpFormatter)
        parser.add_argument("gff_input",
            help="Path to your GFF file to search or edit.")
        parser.add_argument("-sf", "--sequence_file",
            default=None,
            help="Path to corresponding FASTA sequence file (if not in GFF file).\nDefault: 'input filepath stem with fasta/fna/fa extension.'")
        parser.add_argument("-o", "--output_file_name",
            help="Output file path.\nDefault: 'None (print to screen) or input filepath stem with new extension (see -f).' ",
            default=None)
        parser.add_argument("-f", "--output_format",
            default=None,
            help="Output format. Used as extension in output file path. Must be one of:\n\
            none (default): The ID stat for each feature\n\
            coords: tab delimited ID, contig no., start and end of each feature\n\
            number: description\n\
            index: A unique index assigned to each\n\
            bed: bed format\n\
            gff: general feature format\n\
            stats: feature statistics\n\
            subset: feature statistics\n\
            tab | table: tabular format, delimited by tab characters\n\
            subtab | subset_table: tabular format, delimited by tab characters\n\
            fa | fna: contig/scaffold nucleotide sequences in multi-FASTA format\n\
            fasta | ffn: feature nucleotide sequences in multi-FASTA format\n\
            protein | faa: feature protein sequences in multi-FASTA format")
        parser.add_argument("-s", "--search",
            default=None,
            help="Feature search term.")
        parser.add_argument("-c", "--contigs",
            default=None,
            nargs="+",
            help="Name or # of contig/scaffolds to search or extract (see README.txt).\n\
            Takes 1-3 ' '-delimited strings in the format 'Contigs' 'Start' 'Stop'.\n\
            'Contigs' takes a ';'-delimited string of contigs names.\n\
            A string starting ';' is read as contigs to exclude.\n\
            Integers± (e.g 200+) are read as min/max. length limits (bp).\n\
            'Start'/'Stop' (optional) take integer coords of sub-regions to extract.\n\
            Negative 'Stop' values are read as X-bp from end pos.\n\
            Start without Stop is read as 'X-bp from each end'.\n\
            Default 'None'")
        parser.add_argument("-n", "--nearby", "--neighborhood",
            default=None,
            nargs="+",
            help="Include nearby syntenic features in search output\n\
            Takes 1-3 ' '-delimited strings in the format 'Method' '#1' '#2'.\n\
            Method takes a single letter or digit from 0-9 (see README.txt).\n\
            #1/#2 (optional) specify the # ±flanking features to include.")
        parser.add_argument("-d", "--output_delimiter",
            default="\t",
            help="Output delimiter.")
        parser.add_argument("-fh", "--fasta_header",
            default="ID",
            help="A ';' delimited subset of stats to display or include in feature FASTA headers.\nDefault: 'ID'.")
        parser.add_argument("-us", "--upstream",
            default=0,
            help="The length (bp) of upstream sequence to include (no. AAs for protein sequences).\nDefault: '0'.")
        parser.add_argument("-ds", "--downstream",
            default=0,
            help="The length (bp) of downstream sequence to include (no. AAs for protein sequences).\nDefault: '0'.")
        parser.add_argument("-fs", "--fasta_split_every",
            default=60,
            help="The no. characters per line in FASTA sequence output (must be an integer).\nDefault: '60'.")
        parser.add_argument("-u", "--update_stats",
            action="store_true",
            help="Calculate additional statistics for feature entries?\nDefault: 'False'")
        parser.add_argument("-m", "--metadata",
            default=None,
            help="Path to corresponding genome metadata file.\nDefault: 'False'")
        parser.add_argument("-a", "--add_data",
            default=None,
            nargs="+",
            help="Add new feature data from tabular input file?\n\
            Takes 2-3 ' ' delimited strings in the format 'File_path' 'Analysis-type' 'Columns'.\n\
            Analysis-type takes one of 'general', 'pangenome' or 'blast_subject'.\n\
            Columns (optional) takes a ','-delimited string of column # ranges, e.g.'1;3;5-8;11'.")
        parser.add_argument("-x", "-ow", "--overwrite",
            action="store_true",
            help="Allow overwriting if named output file exists?\nDefault: 'False'")
        parser.add_argument("-r", "--rename_contigs",
            nargs="?",
            const=True,
            default=False,
            help="Customise how contigs are renamed in output.\n\
            Default: Retain input contig names in output files.\n\
            Input: 'n/number/#', append contig # to existing names (if absent).\n\
            Input: 'NEW_NAME', appends contig # to NEW_NAME provided.\n\
            No Input: appends contig # to input file name (excluding extension).")
        parser.add_argument("-ID", "--ID_stat",
            default='ID',
            help="Unique identifer statistic for features in the GFF file.\nDefault: 'ID'")
        parser.add_argument("-ff", "--filter_features",
            default=[],
            nargs="+",
            type=str,
            help="Arithmetic operation to filter features.\n\
            Takes 3 ' ' delimited strings in the format 'stat' 'operation' 'value'.\n\
            Operation must be one of '==', '!=', '>=', or '<='.")
        parser.add_argument("-wf", "--without_fasta",
            action="store_false",
            help="Do not inclue FASTA sequence in GFF output file.")
        parser.add_argument("-q", "--quiet_mode_off",
            action="store_true",
            help="Print error and warning messages verbosely.")
        return parser.parse_args()

    def write_region(out_region,out_format,rstart=0,rstop=0,output_file=None,rename_in_fasta=False,fasta_split=60):
        def return_string(input_data,output_file=output_file):
            string = str(input_data)
            if output_file:
                output_file.write(string)
                output_file.write('\n')
            else:
                print(string)
        out_contig = gff_input.contigs.get(out_region)
        use_region = (rstart,rstop) if (rstop + rstart) != 0 else False
        out_sequence = out_region.print_sequence(reverse_strand=(rstart>rstop),split_every=fasta_split,region=use_region,rename=rename_in_fasta)
        return_string(out_sequence)

    def write_feature(my_out,out_format,output_file=None,fasta_split=60):
        def return_string(input_data,output_file=output_file):
            string = str(input_data)
            if output_file:
                output_file.write(string)
                output_file.write('\n')
            else:
                print(string)
        if not out_format:
            return_string(my_out.ID)
        elif out_format == 'coords':
            return_string(str(my_out.ID) + '\t' + my_out.coords.replace('~','\t'))
        elif out_format == 'features':
            return_string(my_out.feature_info)
        elif out_format == 'number' or out_format == 'n' or out_format == 'Number' or out_format == 'N' or out_format == 'index' or out_format == 'idx':
            return_string(my_out.idx)
        elif out_format == 'bed':
            return_string(my_out.write(bed=True))
        elif out_format == 'gff':
            return_string(my_out.raw_entry)
        elif out_format == 'stats':
            for key,value in my_out.feature_info.items():
                return_string(f'{out_delimiter}'.join([key,str(value)]))
        elif out_format == 'subset':
            for stat in fasta_header.split(';'):
                if stat in my_out.feature_info:
                    return_string(f'{out_delimiter}'.join([stat,str(my_out.feature_info.get(stat))]))
        elif out_format == 'tab' or out_format == 'table':
            return_string(f'{out_delimiter}'.join(['Genome','Contig','Start','Stop'] + list(my_out.feature_info.keys())))
            return_string(f'{out_delimiter}'.join([gff_input.name] + my_out.coords.split('~') + [str(val) for val in my_out.feature_info.values()]))
        elif out_format == 'subtab' or out_format == 'subset_table':
            stat_subset = {key:my_out.feature_info.get(key) for key in fasta_header.split(';') if key in my_out.feature_info}
            return_string(f'{out_delimiter}'.join(key for key in stat_subset.keys()))
            return_string(f'{out_delimiter}'.join([str(val) for val in stat_subset.values()]))
        elif out_format == 'fasta' or out_format == 'ffn' or out_format == 'fa' or out_format == 'fna':
            return_string(my_out.print_sequence(gff_input.contigs.get(my_out.contig_name),split_every=fasta_split,fasta_name_stats=fasta_header.split(';'),us=us,ds=ds))
        elif out_format == 'protein' or out_format == 'faa':
            return_string(my_out.print_sequence(gff_input.contigs.get(my_out.contig_name),fasta_header.split(';'),split_every=fasta_split,us=us,ds=ds,protein=True))
        elif out_format in my_out.feature_info:
            return_string(my_out.feature_info.get(out_format))
        else:
            raise ValueError(f'No info {out_format} for {my_out}')

    # collect user arguments
    args = parse_args()
    gff_input = args.gff_input
    out_format = args.output_format
    out_file = args.output_file_name
    overwrite = args.overwrite
    search_feature_info = args.search
    include_nearby_features = args.nearby
    out_delimiter = args.output_delimiter
    fasta_header = args.fasta_header
    feature_analysis_data = args.add_data
    genome_metadata = args.metadata
    update_gff_stats = args.update_stats
    corresponding_fasta = args.sequence_file
    feature_ID_stat = args.ID_stat
    rename_contigs = args.rename_contigs
    contig_parameters = args.contigs
    filtering = args.filter_features
    verbose = args.quiet_mode_off
    without_fasta = args.without_fasta
    rstart,rstop = (0,0) # set rstart and rstop to 0 without contig parameters
    # parse contig parameters
    if contig_parameters:
        filter_contigs = contig_parameters[0]
        rstart,rstop = (contig_parameters[1],contig_parameters[-1]) if len(contig_parameters) > 1 else (0,0)
    # parse integers
    try:
        us = int(args.upstream)
        ds = int(args.downstream)
        if contig_parameters:
            rstart,rstop = (int(rstart),int(rstart))
        split_every = int(args.fasta_split_every)
    except:
        raise TypeError("split_every, up/downstream and conting start/stop values can only be numeric integers.")
    # parse contig sub-region extraction details from include_nearby_features input parameter
    if include_nearby_features:
        synteny_method_defaults = {
            'B':(5,0), 'b':(5,0), '-':(5,0), 'A':(0,5), 'a':(0,5), '+':(0,5),  
            'U':(us,0), 'u':(us,0), '_':(us,0), 'D':(0,ds), 'd':(0,ds), '=':(0,ds),
            'n':(0,1), 'p':(1,0), 'N':(0,1), 'P':(1,0), 
            'F':(7,7), 'f':(7,7), 'S':(7,7), 's':(7,7), '~':(7,7)
        }
        synteny_method = include_nearby_features[0][0] # Extract first character only
        # I could add a clause for 'if output_format in ['fa', 'fna'] here!?
        if synteny_method not in list('AaBbUuDdNnPpFfSs-+_=~'):
            synteny_method = 'f'
            include_nearby_features = [synteny_method] + include_nearby_features
        if len(include_nearby_features) == 1 and synteny_method in list('-+_=~'):
            include_nearby_features = [synteny_method] + include_nearby_features
        if len(include_nearby_features) > 1:
            try:
                n_low = include_nearby_features[1]
                n_low = ''.join([i for i in list(n_low) if i in list('0123456789')])
                n_low = int(n_low)
            except:
                raise TypeError(f'Cannot determine contig region: {n_low} is not a usable method or a numeric integer.')
            if len(include_nearby_features) > 2:
                try:
                    n_high = int(include_nearby_features[2])
                except:
                    raise TypeError(f'Cannot determine contig region: {n_high} is not a usable method or a numeric integer.')
            else:
                if synteny_method in list('BbUuPp-_'):
                    n_high = 0 # if there's only a single value, make sure it's n_low
                elif synteny_method in list('AaDdNn+='):
                    n_high = n_low # if there's only a single value, make sure it's n_high
                    n_low = 0
                else:
                    n_high = n_low # if there's only a single value, assume it's for both
        else:
            n_low, n_high = synteny_method_defaults.get(synteny_method) # search for defaults
        n_in_bp = True if max(n_low,n_high) > 500 else False
        n_ftype = 'CDS' if synteny_method.lower() == synteny_method else None
        n_self = False if synteny_method in list('NnPp') else True
        n_stranded = True if synteny_method in list('UuDdSs~_=') else False
    # parse filtering information from filtering input parameter
    if filtering:
        if not filtering[0] in gff_input.all_recorded_stats:
            raise ValueError("{} information is not recorded for gff features.".format(filtering[0]))
        if not filtering[1] in ['==','!=','>=','<=']:
            raise ValueError("Operation must be one of ['==', '!=', '>=', or '<=']".format(filtering[2]))
        if not type(filtering[2]) == int:
            raise TypeError("{} is not a numeric integer".format(filtering[2]))
    # parse data from filtering 'add feature data' parameter
    if feature_analysis_data:
        analysis_types = ['general', 'pangenome', 'blast_subject']
        if len(feature_analysis_data) == 1:
            feature_analysis_data.append('general')
        if not feature_analysis_data[1] in analysis_types:
            raise ValueError("Operation must be one of {}, not {}".format(str(analysis_types),feature_analysis_data[1]))
        if len(feature_analysis_data) < 3:
            feature_analysis_data.append(None)

    # read GFF input file
    validate_input_file(gff_input, enforced_extension='.gff')
    gff_input = GFF(gff_input.rstrip('/'),ID_stat=feature_ID_stat,update_feature_stats=update_gff_stats,alt_fasta_file=corresponding_fasta)
    if rename_contigs:
        gff_input.rename_contigs(rename_contigs=rename_contigs)
    if genome_metadata:
        try:
            gff_input.add_metadata(genome_metadata)
        except:
            warnings.warn("Could not add genome metadata.")
    if feature_analysis_data:
        gff_input.add_feature_data(feature_analysis_data[0],analysis_type=feature_analysis_data[1],only_extract_columns=feature_analysis_data[2],pad_missing=True)
    
    # generate output file name
    if not out_file and not search_feature_info:
        # use input name if none specified (default)
        out_format = out_format if out_format else '.gff'
        out_format_ext = out_format + '.txt' if out_format in gff_input.all_recorded_stats else out_format
        out_file = os.path.splitext(args.gff_input)[0] + '.' + out_format_ext # use input name if none specified (default)
    if out_file:
        out_path = os.path.dirname(out_file)
        if out_path and not os.path.exists(out_path):
            os.makedirs(out_path)
        if os.path.isdir(out_file):
        # append input name to path if a path is specified
            if not os.path.exists(out_file):
                os.makedir(out_file)
            out_format_ext = out_format if out_format else '.gff'
            out_format_ext = out_format_ext + '.txt' if out_format in gff_input.all_recorded_stats else out_format_ext
            out_file = os.path.join(out_file,(os.path.splitext(args.gff_input)[0] + '.' + out_format_ext))
        if os.path.exists(out_file):
            # if the output file exists, the script must be in overwrite mode to overwrite
            if not overwrite:
                raise RuntimeError("Output file name '{}' matches existing file. Use '-x' to force overwrite. Exiting...".format(out_file))

    # perform contig filtering
    if contig_parameters:
        if re.search('^[0-9]*[-+]$',filter_contigs):
            filter_contigs_as_int = int(filter_contigs.rstrip('-').rstrip('+'))
            if filter_contigs.endswith('-'):
                filtered_contig_list = [ID for ID,contig in gff_input.contigs.items() if contig.length <= filter_contigs_as_int]
            else:
                filtered_contig_list = [ID for ID,contig in gff_input.contigs.items() if contig.length >= filter_contigs_as_int]
        else:
            filtered_contig_list = filter_contigs.lstrip(';').split(';')
            if not filter_contigs.startswith(';'):
                filtered_contig_list = [contig for contig in gff_input.contigs if str(contig) not in filtered_contig_list]
    else:
        filtered_contig_list = []
    # gff_input.contigs = {contig_ID:contig for contig in gff_input.contigs if contig_ID not in filtered_contig_list}
    # gff_input.features = {feature_ID:contig for contig in gff_input.features if feature.contig_name not in filtered_contig_list}
    # gff_input.features = {feature_ID:contig for contig in gff_input.features if feature.contig_number not in filtered_contig_list}

    # conduct a single search if the search parameter was set - can then end script early
    if search_feature_info:
        # conduct the search
        found_features = gff_input.feature(search_feature_info)
        if not found_features:
            found_features = gff_input.search_features(search_feature_info,regex=True)
        if not found_features:
            raise Exception("Nothing matched by search.")
        # format the search result as a list - important for downstream output processing
        if type(found_features) != list:
            found_features = [found_features]
        remember_strand = found_features[0].strand
        # parse and apply synteny_method specs - uses found features list extracted above
        if include_nearby_features:
            if len(found_features) != 1:
                raise Exception(f'Search for {search_feature_info} matches multiple features. Can only extract region surrounding a single feature:')
            found_features = gff_input.feature_region(found_features[0].ID,feature_type=n_ftype,bf=n_low,af=n_high,in_bp=n_in_bp,strand_aware=n_stranded,include_self=n_self)
            if synteny_method in list('PpNn'):
                found_features = [found_features[0]]
            # placeholder - add mode to extract region between two searched features 
            if not out_format:
                out_format = 'bed' # return bed format by default
        if out_file:
            out_file = open(out_file,'a')
        if out_format in ['fa','fna']:
            contig_with_nearby_features = found_features[0].contig_name
            n_contig = gff_input.contigs.get(contig_with_nearby_features)
            out_region_low = min(map(int, found_features[0].coords.split('~')[1:]))
            out_region_high = max(map(int, found_features[-1].coords.split('~')[1:]))
            out_region_low,out_region_high = (out_region_high,out_region_low) if remember_strand == '-' else (out_region_low,out_region_high)
            # print([found_features,str(out_region_low),str(out_region_high)])
            # print([found_features[0].ID,n_contig,str(out_region_low),str(out_region_high)])
            out_region_tofasta = gff_input.contigs.get(n_contig)
            new_header = fasta_header if fasta_header not in gff_input.all_recorded_stats else None
            write_region(n_contig,out_format,rstart=out_region_low,rstop=out_region_high,output_file=out_file,rename_in_fasta=new_header,fasta_split=split_every)
        else:
            for feature in found_features:
                write_feature(feature,out_format,output_file=out_file,fasta_split=split_every)
        if out_file:
            out_file.close()
        sys.exit(0)

    if out_format == 'fa' or out_format == 'fna':
        with open(out_file,'a') as outfile:
            for i in range(1,gff_input.contig_count+1):
                current_contig = gff_input.contigs.get(i)
                if current_contig.name in filtered_contig_list or current_contig.number in filtered_contig_list:
                    continue
                # if rstop == 0:
                #     rstop = current_contig.length 
                # if str(rstop).startswith('-') or rstart == rstop:
                #     rstop = current_contig.length + rstop
                # print(rstart,rstop)
                outfile.write(current_contig.print_sequence(reverse_strand=False,split_every=split_every,region=False,rename=False)+'\n')
    elif not out_format or out_format == 'gff':
        gff_input.to_newfile(out_file=out_file,
            rename_contigs=False,update_stats=update_gff_stats,
            add_FASTA_sequence=without_fasta,FASTA_Split_Every=split_every,
            skip_entries=filtering,skip_contigs=filtered_contig_list)
    elif out_format in ['subset', 'subset_table', 'subtab', 'features']:
        raise ValueError("out_format {} cannot be applied.".format(out_format))
    elif out_format in ['index','number','coords','ffn','fasta','stats','protein','faa','bed','tab','table']:
        with open(out_file,'a') as outfile:
            for feature in gff_input.features.values():
                if feature.contig_name in filtered_contig_list or feature.contig_number in filtered_contig_list:
                    continue
                if feature.seq_type != 'CDS' or feature.ID in filtering:
                    continue
                write_feature(feature,out_format,output_file=outfile,fasta_split=split_every)
    else:
        print("No parameters detected - printing file statistics:")
        print(gff_input.feature_type_dict)
        print(gff_input.all_recorded_stats)
