# def load_StoreObj(parameterObj):
#     #storeObj = StoreObj()

#     return StoreObj()

# class StoreObj(object):
#     def __init__(self):
#         self.zarr_dir = None
#         self.dataset = None
#         self.bed_f = None
#         self.sample_ids_by_pop_id = {} # 'pop1', 'pop2', ...
#         self.pop_ids_order = []
#         #self.outprefix = parameterObj.outprefix
#         #if parameterObj.zarr_dir:
#         #    self._load_store(parameterObj)
#         #if parameterObj.vcf_file:
#         #    self._initiate_store(parameterObj)
#         #if parameterObj.bed_file:
#         #    self._parse_intervals(parameterObj)

#     def parse_sample_file(self, parameterObj):
#         # parse sample CSV
#         samples_df = pd.read_csv(\
#             parameterObj.sample_file, \
#             sep=",", \
#             usecols=[0, 1], \
#             names=['sample_id', 'population_id'], \
#             header=None, \
#             dtype={ \
#                 'sample_id': 'category', \
#                 'population_id': 'category' \
#                 } \
#             )
#         # Error if not 2 populations
#         if not len(samples_df.groupby('population_id').count().index) == 2:
#              exit('[X] Invalid number of populations: %s (must be 2)' % len(samples_df.groupby('population_id').count().index))
#         sorted_samples_df = samples_df.sort_values(['population_id', 'sample_id'])
#         print(sorted_samples_df)
#         population_idx = -1
#         sample_ids_by_population_id = collections.defaultdict(list)
#         for sample_idx, (sample_id, population_id) in enumerate(sorted_samples_df.values.tolist()):
#             if not population_id in self.population_idx_by_population_id:
#                 population_idx += 1
#                 self.add_population(population_idx, population_id)
#             self.add_sample_to_populationObj(population_id, sample_id)
#             sample_ids_by_population_id[population_idx].append(sample_id)
#         for pair_idx, pair_id in enumerate([(x) for x in itertools.product(*sorted(sample_ids_by_population_id.values()))]):
#             self.add_pair(pair_idx, pair_id)

#     def parse_genome_file(self, parameterObj):
#         df = pd.read_csv(parameterObj.genome_file, sep="\t", usecols=[0, 1], names=['sequence_id', 'length'], header=None)
#         for sequence_idx, (sequence_id, length_str) in enumerate(df.values.tolist()):
#             try:
#                 length = int(length_str)
#             except TypeError:
#                 exit("[X] Line %s: Second column of --genome_file %s must be integers, not '%s'" % (sequence_idx, parameterObj.genome_file, length_str))
#             self.add_sequence(sequence_idx, sequence_id, length)


#     def tree(self):
#         return self.dataset.tree()

#     def _initiate_store(self, parameterObj):
#         self.zarr_dir = str(pathlib.Path(parameterObj.outprefix).with_suffix('.zarr'))
#         self.dataset = self._get_dataset(overwrite=True)
#         self._parse_data(parameterObj)
#         self.sample_ids_by_pop_id = self._get_sample_ids_by_pop_id()
#         self.pop_ids_order = sorted(set([pop_id for pop_id in self.yield_pop_ids()]))

#     def _parse_intervals(self, parameterObj):
#         logging.info("[#] Processing BED file ...")
#         print(parameterObj.bed_file)
#         bed_df = pd.read_csv( 
#             parameterObj.bed_file, 
#             sep="\t", 
#             usecols=[0, 1, 2, 4], 
#             names=['sequence_id', 'start', 'end', 'samples'], 
#             skiprows=1, 
#             header=None, 
#             dtype={ 
#                 'sequence_id': 'category', 
#                 'start': np.int, 
#                 'end': np.int, 
#                 'samples': 'category'})
#         # MIGHT NOT BE NECESSARY
#         # filter rows based on sequence_ids, sort by sequence_id, start
#         # bed_df = bed_df[bed_df['sequence_id'].isin(entityCollection.sequence_idx_by_sequence_id)].sort_values(['sequence_id', 'start'], ascending=[True, True])
#         # compute length 
#         bed_df['length'] =  bed_df['end'] - bed_df['start'] 
#         print(bed_df)
#         bed_length_total = int(bed_df['length'].sum())
#         print(bed_length_total)
#         # compute cov-matrix (sample_ids are column names)
#         cov_df = bed_df.samples.str.get_dummies(sep=',')
#         # for each combination of sample_ids ...
#         print(self.sample_ids_by_pop_id)
#         for a, b in itertools.combinations(self.sample_ids_by_pop_id['all'], 2):
#             pair_id = frozenset(a, b)
#             print(pair_id)
#             boolean_mask = cov_df[[a, b]].all(axis='columns')
#             interval_space = interval_df[boolean_mask]
#             block_arrays = np.split(interval_space, block_length)

#         # https://stupidpythonideas.blogspot.com/2014/01/grouping-into-runs-of-adjacent-values.html
#         # https://stackoverflow.com/questions/2154249/identify-groups-of-continuous-numbers-in-a-list
#         # https://stackoverflow.com/questions/4494404/find-large-number-of-consecutive-values-fulfilling-condition-in-a-numpy-array
#         a = np.array((bed_df.start, bed_df.end)).T
#         r = create_ranges(a)
#         def consecutive(data, stepsize=1):
#             return np.split(data, np.where(np.diff(data) != stepsize)[0]+1)
#         a = np.array([0, 47, 48, 49, 50, 97, 98, 99])
#         consecutive(a)

#     def _load_store(self, parameterObj):
#         self.zarr_dir = parameterObj.zarr_dir
#         self._add_sample_information(parameterObj)
#         self.dataset = self._get_dataset()
#         self.sample_ids_by_pop_id = self._get_sample_ids_by_pop_id()
#         self.pop_ids_order = sorted(set([pop_id for pop_id in self.yield_pop_ids()]))

#     def yield_pop_ids(self):
#         self._yield('/pop_ids')

#     def yield_sample_ids(self):
#         self._yield('/pop_ids')

#     def yield_seq_ids(self):
#         self._yield('/seq_ids')

#     def yield_seq_lengths(self):
#         self._yield('/seq_lengths')

#     def _yield(self, key):
#         for value in list(self.dataset[key]):
#             yield value

#     def _get_sample_ids_by_pop_id(self):
#         sample_ids_by_pop_id = collections.defaultdict(list)
#         sample_ids_by_pop_id['all'] = list(range(len(self.sample_ids)))
#         for sample_id, pop_id in zip(range(len(self.sample_ids)), self.pop_ids):
#             sample_ids_by_pop_id[pop_id].append(sample_id)
#         return sample_ids_by_pop_id

#     def _get_dataset(self, overwrite=False):
#         if overwrite:
#             if os.path.isdir(self.zarr_dir):
#                 logging.info("[!] ZARR store %r exists. Deleting ..." % self.zarr_dir)
#                 shutil.rmtree(self.zarr_dir)
#             logging.info("[+] Generating ZARR store %r" % self.zarr_dir)
#             return zarr.open(self.zarr_dir, mode='w')
#         return zarr.open(self.zarr_dir, mode='a')

#     def _parse_data(self, parameterObj):
#         self.sample_ids = self._add_sample_ids(allel.read_vcf_headers(parameterObj.vcf_file).samples)
#         seq_ids = self._add_sequences_from_header(allel.read_vcf_headers(parameterObj.vcf_file).headers)
#         self._add_variants(parameterObj.vcf_file, seq_ids)
#         self._add_sample_information(parameterObj, self.sample_ids)

#     def _add_sample_ids(self, sample_ids):
#         self.dataset.create_dataset('sample_ids', data=np.array(sample_ids))
#         return sample_ids

#     def _add_sequences_from_header(self, header_lines):
#         pattern = re.compile(r'##contig=<ID=(\S+),length=(\d+)>')
#         seq_ids, lengths = [], []
#         for header_line in header_lines:
#             if header_line.startswith("##contig"):
#                 seq_id, length = re.match(pattern, header_line).groups()
#                 seq_ids.append(seq_id)
#                 lengths.append(int(length))
#         self.dataset.create_dataset('seq_ids', data=np.array(seq_ids))
#         self.dataset.create_dataset('seq_lengths', data=np.array(lengths))
#         return seq_ids

#     def _add_variants(self, vcf_file, seq_ids):
#         for seq_id in tqdm(seq_ids, total=len(seq_ids), desc="[%] ", ncols=100):
#             allel.vcf_to_zarr( 
#                 vcf_file, 
#                 self.zarr_dir, 
#                 group="%s" % seq_id,
#                 region=seq_id,
#                 fields=[
#                     'variants/POS',
#                     'variants/REF',
#                     'variants/ALT',
#                     'variants/QUAL',
#                     'variants/is_snp',
#                     'variants/numalt',
#                     'variants/DP',
#                     'calldata/GT',
#                     'calldata/DP'
#                 ], 
#                 overwrite=True)

#     def _add_sample_information(self, parameterObj, sample_ids):
#         pop_id_by_sample_id = get_pop_ids_by_sample_ids_from_csv(parameterObj)
#         self.dataset.create_dataset('pop_ids', data=np.array([pop_id_by_sample_id[sample_id] for sample_id in sample_ids]))

# def mutype_counter_to_mutuple(mutype_counter, full=True):
#     if full:
#         return tuple([mutype_counter[mutype] for mutype in FULL_MUTYPE_ORDER])
#     return tuple([mutype_counter[mutype] for mutype in MUTYPE_ORDER])

# def mutype_counter_to_dict(mutype_counter, full=True):
#     order = MUTYPE_ORDER
#     if full:
#         order = FULL_MUTYPE_ORDER
#     return dict({mutype: mutype_counter[mutype] for mutype in order})

# def gt_counter_to_dict(gt_counter): 
#     return dict({gt: gt_counter[gt] for gt in GT_ORDER})

# class EntityCollection(object):
#     def __init__(self):
#         # populations
#         self.population_idx_by_population_id = {}
#         self.populationObjs = []
#         # pairs
#         self.pair_idx_by_pair_id = {}
#         self.pairObjs = []
#         # sequences
#         self.sequence_idx_by_sequence_id = {}
#         self.sequenceObjs = []
#         # blocks
#         self.block_idx_by_block_id = {}
#         self.blockObjs = []
#         # samples
#         self.population_id_by_sample_id = {}
#         # windows
#         self.windowObjs = []
        
#     def __str__(self):
#         output = []
#         output.append("# %s" % str(self.population_idx_by_population_id))
#         for populationObj in self.populationObjs:
#             output.append(str(populationObj))
#         output.append("# %s" % str(self.pair_idx_by_pair_id))
#         for pairObj in self.pairObjs:
#             output.append(str(pairObj))
#         output.append("# %s" % str(self.sequence_idx_by_sequence_id))
#         for sequenceObj in self.sequenceObjs:
#             output.append(str(sequenceObj))
#         output.append("# %s" % str(self.block_idx_by_block_id))
#         for blockObj in self.blockObjs:
#             output.append(str(blockObj))
#         return "%s\n" % "\n".join(output)

#     def count(self, entity_type):
#         if entity_type == 'populations':
#             return len(self.populationObjs)
#         elif entity_type == 'samples':
#             return sum([len(populationObj.sample_ids) for populationObj in self.populationObjs])
#         elif entity_type == 'pairs':
#             return len(self.pairObjs)
#         elif entity_type == 'sequences':
#             return len(self.sequenceObjs)
#         elif entity_type == 'bases':
#             return sum([sequenceObj.length for sequenceObj in self.sequenceObjs])
#         elif entity_type == 'blocks':
#             return len(self.blockObjs)
#         elif entity_type == 'windows':
#             return len(self.windowObjs)
#         else:
#             return 0
    
#     def sample_string_to_sample_ids(self, sample_string):
#         return (sample_id for sample_id in sample_string.split(","))

#     def sample_string_to_pair_idxs(self, sample_string):
#         # works only for two pops ...
#         pair_idxs = frozenset(filter(lambda x: x >= 0, [self.pair_idx_by_pair_id.get(frozenset(x), -1) for x in itertools.combinations(sample_string.split(","), 2)])) 
#         if pair_idxs:
#             return pair_idxs
#         return np.nan

#     def count_pair_idxs(self, pair_idxs):
#         if not pair_idxs is np.nan:
#             return len(list(pair_idxs))
#         return 0

#     def pair_idxs_to_sample_ids(self, pair_idxs):
#         return set(itertools.chain.from_iterable([self.pairObjs[pair_idx].id for pair_idx in pair_idxs]))

#     def add_population(self, population_idx, population_id):
#         self.populationObjs.append(PopulationObj(population_idx, population_id))
#         self.population_idx_by_population_id[population_id] = population_idx

#     def add_pair(self, pair_idx, pair_id):
#         self.pairObjs.append(PairObj(pair_idx, pair_id))
#         self.pair_idx_by_pair_id[frozenset(pair_id)] = pair_idx        

#     def add_sequence(self, sequence_idx, sequence_id, length):
#         self.sequenceObjs.append(SequenceObj(sequence_idx, sequence_id, length))
#         self.sequence_idx_by_sequence_id[sequence_id] = sequence_idx

#     def add_blockObjs(self, blockObjs):
#         blockObjs.sort(key=lambda x: (x.sequence_id, x.start)) 
#         for idx, blockObj in enumerate(blockObjs):
#             if blockObj.span < 64 or blockObj.span > 80:
#                 print(blockObj.id, blockObj.span)
#             self.blockObjs.append(blockObj)
#             self.block_idx_by_block_id[blockObj.id] = len(self.blockObjs) - 1
    
#     def add_sample_to_populationObj(self, population_id, sample_id):
#         population_idx = self.population_idx_by_population_id[population_id]
#         self.populationObjs[population_idx].sample_ids.append(sample_id)
#         self.population_id_by_sample_id[sample_id] = population_id

#     def add_windowObjs(self, windowObjs):
#         self.windowObjs = windowObjs

#     def load_blocks(self, parameterObj, purpose='variants'):
#         block_cols = ['block_id', 'length', 'span', 'sample_ids', 'pair_idxs', 'distance'] # maybe other colums need not be saved ...
#         blocks_hdf5_store = pd.HDFStore(parameterObj.blocks_file)
#         block_bed_df = pd.read_hdf(blocks_hdf5_store, key='bed')   
#         block_df = pd.read_hdf(blocks_hdf5_store, 'block').reindex(columns=block_cols)
#         blocks_hdf5_store.close()
#         block_lengths = block_df.length.unique().tolist()
#         if not len(block_lengths) == 1:
#             exit('[X] Non uniform block length found : %s' % block_lengths)
#         parameterObj.block_length = block_lengths[0]
#         bed_tuples_by_block_id = collections.defaultdict(list)
#         for block_id, sequence_id, start, end in block_bed_df.values.tolist():
#             bed_tuples_by_block_id[block_id].append((sequence_id, start, end))
#         blockObjs = []
#         # parallelisation possible
#         for block_id, length, span, sample_ids, _pair_idxs, distance in tqdm(block_df[block_cols].values.tolist(), total=len(block_df.index), desc="[%] ", ncols=100):
#             if block_id in bed_tuples_by_block_id: # only Blocks in BED file are instanciated!
#                 pair_idxs = [int(pair_idx) for pair_idx in _pair_idxs.split(",")]
#                 blockObj = BlockObj(block_id, parameterObj.block_length)
#                 for bed_tuple in bed_tuples_by_block_id[block_id]:
#                     bedObj = BedObj(bed_tuple[0], bed_tuple[1], bed_tuple[2], pair_idxs, bed_tuple[2] - bed_tuple[1])
#                     blockObj.add_bedObj(bedObj, parameterObj, self)
#                 blockObj.mutype_counter_by_pair_idx = {pair_idx: collections.Counter() for pair_idx in pair_idxs}
#                 blockObjs.append(blockObj)
#             else:
#                 print("[-] Block %s is missing from %s. Block will be ignored." % (block_id, parameterObj.block_bed))
#         self.add_blockObjs(blockObjs)

#         if purpose == 'windows':
#             block_region_ids = (block_df["distance"].fillna(parameterObj.max_block_distance + 1).shift() > float(parameterObj.max_block_distance)).cumsum() # generate indices for splitting
#             block_id_batches = []
#             for idx, block_region_df in block_df.groupby(block_region_ids):
#                 #print("####", block_region_df)
#                 if len(block_region_df) > parameterObj.window_size: # remove regions below window_size
#                     #block_region_df = block_region_df.drop(columns=['length', 'span', 'sample_ids', 'pair_idxs', 'distance'], axis=1) # remove distance/sample_idsx columns
#                     block_id_batches.append(block_region_df['block_id'].tolist())
#             if len(block_id_batches) == 0:
#                 exit("[X] Insufficient consecutive blocks with parameters '--window_size' (%s) and '--max_block_distance' (%s)" % (parameterObj.window_size, parameterObj.max_block_distance) )
#             return block_id_batches

#     def get_mutype_counters_by_block_id(self, parameterObj):
#         variants_hdf5_store = pd.HDFStore(parameterObj.variants_file)
#         variant_cols = ['block_id', 'pair_idx', 'hetA', 'fixed', 'hetB', 'hetAB']
#         variants_block_df = pd.read_hdf(variants_hdf5_store, key='blocks').reindex(columns=variant_cols)
#         variants_hdf5_store.close()
#         mutype_counters_by_block_id = collections.defaultdict(list)
#         for block_id, pair_idx, hetA, fixed, hetB, hetAB in tqdm(variants_block_df.values.tolist(), total=len(variants_block_df.index), desc="[%] ", ncols=100):
#             mutype_counters_by_block_id[block_id].append(collections.Counter({'hetA': hetA, 'fixed': fixed, 'hetB': hetB, 'hetAB': hetAB}))
#         return mutype_counters_by_block_id

#     def transform_coordinates(self, parameterObj, coordinateTransformObj):
#         params = [(blockObj, coordinateTransformObj) for blockObj in self.blockObjs]
#         if parameterObj.threads < 2:
#             with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
#                 for param in params:
#                     #print(blockObj, blockObj.void) 
#                     self.transform_coordinates_blockObj(param)
#                     #print(blockObj, blockObj.void)
#                     pbar.update()
#         else:
#             with poolcontext(processes=parameterObj.threads) as pool:
#                 with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
#                     for blockObj in pool.imap_unordered(self.transform_coordinates_blockObj, params):
#                         pbar.update()

#     def transform_coordinates_blockObj(self, params):
#         blockObj, coordinateTransformObj = params
#         # print(">", blockObj.sequence_id, blockObj.start, blockObj.end)
#         new_sequence_id, new_start, new_end = coordinateTransformObj.transform_coordinates(blockObj.sequence_id, blockObj.start, blockObj.end)
#         # print("<", new_sequence_id, new_start, new_end)
#         new_bed_tuples = []
#         if new_sequence_id is None:
#             blockObj.void = True
#         else:
#             blockObj.sequence_id = new_sequence_id
#             blockObj.start = new_start
#             blockObj.end = new_end
#             for bed_tuple in blockObj.bed_tuples:
#                 new_bed_sequence_id, new_bed_start, new_bed_end = coordinateTransformObj.transform_coordinates(bed_tuple[0], bed_tuple[1], bed_tuple[2])
#                 new_bed_tuples.append((new_bed_sequence_id, new_bed_start, new_bed_end))
#             new_bed_tuples.sort(key=lambda i: (i[0], i[1]))
#             blockObj.bed_tuples = new_bed_tuples

#     def parse_sample_file(self, parameterObj):
#         if not parameterObj.sample_file:
#             exit('[X] File cannot be read: %s' % parameterObj.args['sample_file'])
#         # parse sample CSV
#         samples_df = pd.read_csv(\
#             parameterObj.sample_file, \
#             sep=",", \
#             usecols=[0, 1], \
#             names=['sample_id', 'population_id'], \
#             header=None, \
#             dtype={ \
#                 'sample_id': 'category', \
#                 'population_id': 'category' \
#                 } \
#             )
#         if not len(samples_df.groupby('population_id').count().index) == 2:
#              exit('[X] Invalid number of populations: %s (must be 2)' % len(samples_df.groupby('population_id').count().index))
#         sorted_samples_df = samples_df.sort_values(['population_id', 'sample_id'])
#         population_idx = -1
#         sample_ids_by_population_id = collections.defaultdict(list)
#         for sample_idx, (sample_id, population_id) in enumerate(sorted_samples_df.values.tolist()):
#             if not population_id in self.population_idx_by_population_id:
#                 population_idx += 1
#                 self.add_population(population_idx, population_id)
#             self.add_sample_to_populationObj(population_id, sample_id)
#             sample_ids_by_population_id[population_idx].append(sample_id)
#         for pair_idx, pair_id in enumerate([(x) for x in itertools.product(*sorted(sample_ids_by_population_id.values()))]):
#             self.add_pair(pair_idx, pair_id)

#     def parse_genome_file(self, parameterObj):
#         df = pd.read_csv(parameterObj.genome_file, sep="\t", usecols=[0, 1], names=['sequence_id', 'length'], header=None)
#         for sequence_idx, (sequence_id, length_str) in enumerate(df.values.tolist()):
#             try:
#                 length = int(length_str)
#             except TypeError:
#                 exit("[X] Line %s: Second column of --genome_file %s must be integers, not '%s'" % (sequence_idx, parameterObj.genome_file, length_str))
#             self.add_sequence(sequence_idx, sequence_id, length)

#     #def calculate_variation(self, mutype_counter, sites):
#     #    # if missing, assume invariant
#     #    pi_1 = float("%.8f" % ((mutype_counter['hetA'] + mutype_counter['hetAB']) / sites))
#     #    pi_2 = float("%.8f" % ((mutype_counter['hetB'] + mutype_counter['hetAB']) / sites))
#     #    d_xy = float("%.8f" % ((((mutype_counter['hetA'] + mutype_counter['hetB'] + mutype_counter['hetAB']) / 2.0) + mutype_counter['fixed']) / sites))
#     #    mean_pi = (pi_1 + pi_2) / 2.0
#     #    total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
#     #    f_st = np.nan
#     #    if (total_pi):
#     #        f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
#     #    return pi_1, pi_2, d_xy, f_st

#     def calculate_variation_from_df(self, df, sites):
#         # print(df)
#         pi_1 = float("%.8f" % ((df.hetA.sum() + df.hetAB.sum()) / sites))
#         pi_2 = float("%.8f" % ((df.hetB.sum() + df.hetAB.sum()) / sites))
#         d_xy = float("%.8f" % ((((df.hetA.sum() + df.hetB.sum() + df.hetAB.sum()) / 2.0) + df.fixed.sum()) / sites))
#         mean_pi = (pi_1 + pi_2) / 2.0
#         total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
#         f_st = np.nan
#         if (total_pi):
#             f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
#         return (pi_1, pi_2, d_xy, f_st)

#     def generate_modify_output(self, parameterObj):
#         print("[#] Generating output...")
#         self.generate_block_output(parameterObj, mode='modify')

#     def generate_window_output(self, parameterObj, entityCollection, mutype_counters_by_block_id):
#         mutuples = []
#         #variants_df.set_index(keys=['block_id'], drop=False, inplace=True)
#         mutuple_counters_by_window_id = collections.defaultdict(collections.Counter)
#         window_vals = []
#         for windowObj in tqdm(self.windowObjs, total=len(self.windowObjs), desc="[%] ", ncols=100, unit_scale=True):
#             window_mutuples = []
#             for block_id in windowObj.block_ids:
#                 mutype_counters = mutype_counters_by_block_id[block_id]
#                 for mutype_counter in mutype_counters:
#                     mutype_tuple = mutype_counter_to_mutuple(mutype_counter, full=False)
#                     # window_tally
#                     mutuple_counters_by_window_id[windowObj.id][mutype_tuple] += 1
#                     # global tally
#                     mutuples.append(mutype_tuple)
#                     # window list
#                     window_mutuples.append(mutype_tuple)

#             window_df = pd.DataFrame(window_mutuples, columns=MUTYPE_ORDER)
#             #print(window_df)
#             pi_1, pi_2, dxy, fst = self.calculate_variation_from_df(window_df, sites=(len(window_df.index) * parameterObj.block_length))
#             mean_sample_fraction = (sum(windowObj.sample_counts) / len(self.population_id_by_sample_id)) / parameterObj.window_size
#             window_vals.append([windowObj.id, windowObj.sequence_id, windowObj.start, windowObj.end, windowObj.span, windowObj.centre, pi_1, pi_2, dxy, fst, mean_sample_fraction])
        
#         print("[#] Creating dataframe of window mutuple tallies...")
#         window_mutuple_tally_cols = ['window_id', 'count'] + MUTYPE_ORDER 
#         window_mutuple_tally_vals = []
#         #print(mutuple_counters_by_window_id)
#         for window_id, mutype_counters in mutuple_counters_by_window_id.items():
#             #print(window_id, mutype_counters)
#             for mutype, count in mutype_counters.most_common():
#                 window_mutuple_tally_vals.append([window_id, count] + list(mutype))
#         window_mutuple_tally_df = pd.DataFrame(window_mutuple_tally_vals, columns=window_mutuple_tally_cols)

#         print("[#] Creating dataframe of window metrics...")
#         window_cols = [
#             'window_id',
#             'sequence_id',
#             'start',
#             'end', 
#             'span',
#             'centre',
#             'pi_%s' % (self.populationObjs[0].id),
#             'pi_%s' % (self.populationObjs[1].id),
#             'dxy',
#             'fst',
#             'sample_cov'
#             ]
#         window_df = pd.DataFrame(window_vals, columns=window_cols)
#         plot_fst_genome_scan(window_df, '%s.fst_genome_scan.png' % parameterObj.dataset, self.sequenceObjs)
#         plot_pi_genome_scan(window_df, '%s.pi_genome_scan.png' % parameterObj.dataset, self.sequenceObjs)
#         plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)
#         # storage
#         window_hdf5_store = create_hdf5_store(
#             out_f='%s.windows.h5' % (parameterObj.prefix), 
#             path=parameterObj.path
#             )
#         window_df.to_hdf(window_hdf5_store, 'window_metrics', append=True)
#         window_mutuple_tally_df.to_hdf(window_hdf5_store, 'mutypes', append=True)
#         # tabulate_df(window_df, columns=window_cols, title="Windows")
#         window_hdf5_store.close()

#     def generate_variant_output(self, parameterObj):
#         # collate data
#         block_vals = []
#         counts_by_site_by_sample_id_by_population_id = collections.defaultdict(dict)
#         mutuples = []
#         print("[#] Analysing variants in blocks...")
#         for blockObj in tqdm(self.blockObjs, total=len(self.blockObjs), desc="[%] ", ncols=100, unit_scale=True):
#             for sample_id in blockObj.sample_ids:
#                 population_id = self.population_id_by_sample_id[sample_id]
#                 gt_counter = blockObj.gt_counter_by_sample_id.get(sample_id, collections.Counter())
#                 if not sample_id in counts_by_site_by_sample_id_by_population_id[population_id]:
#                     counts_by_site_by_sample_id_by_population_id[population_id][sample_id] = collections.Counter()
#                 counts_by_site_by_sample_id_by_population_id[population_id][sample_id] += gt_counter
#             for pair_idx in blockObj.pair_idxs:
#                 mutype_counter = blockObj.mutype_counter_by_pair_idx.get(pair_idx, collections.Counter())
#                 mutuples.append(mutype_counter_to_mutuple(mutype_counter, full=False))
#                 block_vals.append([blockObj.id, pair_idx] + list(mutype_counter_to_mutuple(mutype_counter, full=True)))

#         global_mutype_counter = collections.Counter(mutuples)
#         print("[+] Monomorphic blocks: %s (%s of blocks)" % (global_mutype_counter[(0, 0, 0, 0)], format_percentage(global_mutype_counter[(0, 0, 0, 0)] / len(mutuples) )))
#         print("[#] Creating dataframe of global mutuple tallies...")
#         global_mutuple_tally_cols = ['count'] + MUTYPE_ORDER 
#         global_mutuple_tally_vals = []
#         for mutype, count in global_mutype_counter.most_common():
#             global_mutuple_tally_vals.append([count] + list(mutype))
#         global_mutuple_tally_df = pd.DataFrame(global_mutuple_tally_vals, columns=global_mutuple_tally_cols)
#         # Mutype plot
#         plot_mutuple_barchart(
#             '%s.mutuple_barchart.png' % parameterObj.dataset,
#             global_mutype_counter
#             )
#         print("[#] Creating dataframe of blocks...")
#         block_idx = ['block_id', 'pair_idx']
#         block_df = pd.DataFrame(block_vals, columns=block_idx + FULL_MUTYPE_ORDER)
#         # add four gamete violations (FGV)
#         block_df['FGV'] = np.where(( (block_df['fixed'] > 1) & (block_df['hetAB'] > 1) ), True, False) 

#         print("[#] Creating dataframe for samples...")
#         sample_df = pd.concat(
#             {population_id: pd.DataFrame(dict_of_dicts).T for population_id, dict_of_dicts in counts_by_site_by_sample_id_by_population_id.items()}, 
#             axis=0, 
#             names=['population_id', 'sample_id'], 
#             sort=False).fillna(0).astype(int)
#         sample_df['sites'] = sample_df['blocks'] * parameterObj.block_length
#         sample_df['heterozygosity'] = sample_df['HET'] / sample_df['sites']
#         sample_df['missingness'] = sample_df['MISS'] / sample_df['sites']
#         sample_df.drop(columns=['HET', 'HOM', 'MISS'], axis=1, inplace=True)
#         tabulate_df(sample_df, columns = (sample_df.index.names + list(sample_df.columns)), title="Sample metrics: %s" % parameterObj.args.get('--prefix', 'gimble'))        

#         print("[#] Creating dataframe for populations...")
#         population_df = sample_df.groupby('population_id').mean()
#         tabulate_df(population_df, columns=['population_id'] + list(population_df.columns), title="Population metrics: %s" % parameterObj.args.get('--prefix', 'gimble'))        

#         print("[#] Creating dataframe for dataset...")
#         #print(block_df)
#         mutype_df = block_df.drop(columns=MUTYPE_OTHER + ['FGV'], axis=1)[block_idx + MUTYPE_ORDER].set_index(['block_id', 'pair_idx'])
#         #print(mutype_df)
#         pi_1, pi_2, dxy, fst = self.calculate_variation_from_df(mutype_df, sites=(parameterObj.block_length * len(mutype_df.index)))
#         #print(parameterObj.block_length*len(mutype_df))
#         #print((len(self.blockObjs) * parameterObj.block_length))
#         sites = (len(self.blockObjs) * parameterObj.block_length)
#         dataset_cols = [
#             'blocks',
#             'sites',
#             'FGVs',
#             'pi_%s' % (self.populationObjs[0].id),
#             'pi_%s' % (self.populationObjs[1].id),
#             'dxy',
#             'fst'
#             ]
#         dataset_vals = [
#             len(self.blockObjs),
#             sites,
#             (len(block_df[block_df['FGV'] == True]) / sites),
#             pi_1,
#             pi_2,
#             dxy,
#             fst
#         ]
#         dataset_df = pd.DataFrame.from_dict(dict(zip(dataset_cols, dataset_vals)), orient='index').T
#         tabulate_df(dataset_df, columns=dataset_cols, title="Dataset metrics: %s" % parameterObj.args.get('--prefix', 'gimble'))        
        
#         # storage
#         variant_hdf5_store = create_hdf5_store(
#             out_f='%s.variants.h5' % (parameterObj.prefix), 
#             path=parameterObj.path, 
#             )
#         population_df.to_hdf(variant_hdf5_store, 'populations', append=True)
#         sample_df.to_hdf(variant_hdf5_store, 'samples', append=True)
#         block_df.to_hdf(variant_hdf5_store, 'blocks', append=True)
#         global_mutuple_tally_df.to_hdf(variant_hdf5_store, 'mutypes', append=True)
#         dataset_df.to_hdf(variant_hdf5_store, 'dataset', append=True)
#         variant_hdf5_store.close()

        
#         #for mutype in FULL_MUTYPE_ORDER:
#             #print(mutype, np.unique(block_df[mutype], return_index=False, return_inverse=False, return_counts=True, axis=0))
#         #non_missing_mutypes = block_df.loc[block_df['missing'] <= 4,:][MUTYPE_ORDER]
#         #values, counts = np.unique(non_missing_mutypes[MUTYPE_ORDER].values, return_index=False, return_inverse=False, return_counts=True, axis=0)
        
#     def generate_block_output(self, parameterObj, mode='blocks'):
#         block_bed_cols = ['block_id', 'sequence_id', 'bed_start', 'bed_end']
#         block_bed_vals = []
#         block_cols = ['block_id', 'sequence_id', 'block_start', 'block_end', 'length', 'span', 'sample_ids', 'pair_idxs', 'count_samples', 'count_pairs']
#         block_vals = []
#         bases_blocked_by_sequence_id = collections.Counter()
#         bases_blocked_by_sample_id = collections.Counter()
#         bases_blocked_by_pair_id = collections.Counter()
#         bases_blocked_by_sample_count = collections.Counter()
#         void_count = 0
#         for blockObj in tqdm(self.blockObjs, total=len(self.blockObjs), desc="[%] ", ncols=100, unit_scale=True):
#             if blockObj.void:
#                 void_count += 1
#             else:
#                 block_vals.append([
#                     blockObj.id, 
#                     blockObj.sequence_id,
#                     blockObj.start,
#                     blockObj.end,
#                     blockObj.length, 
#                     blockObj.span, 
#                     ",".join([str(x) for x in sorted(blockObj.sample_ids)]), 
#                     ",".join([str(x) for x in sorted(blockObj.pair_idxs)]),
#                     len(blockObj.sample_ids), 
#                     len(blockObj.pair_idxs)
#                     ])
#                 for bed_tuple in blockObj.bed_tuples:
#                     block_bed_vals.append([blockObj.id, bed_tuple[0],bed_tuple[1],bed_tuple[2]])
#                 bases_blocked_by_sequence_id[blockObj.sequence_id] += parameterObj.block_length
#                 # collect base counts
#                 for sample_count in range(len(blockObj.sample_ids), 0, -1):
#                     bases_blocked_by_sample_count[sample_count] += parameterObj.block_length
#                 for sample_id in blockObj.sample_ids:
#                     bases_blocked_by_sample_id[sample_id] += parameterObj.block_length
#                 for pair_idx in blockObj.pair_idxs:
#                     bases_blocked_by_pair_id[self.pairObjs[pair_idx].id] += parameterObj.block_length
#         block_bed_df = pd.DataFrame(block_bed_vals, columns=block_bed_cols)
#         #tabulate_df(block_bed_df, columns=block_bed_cols, title="Block BED intervals")
#         block_df = pd.DataFrame(block_vals, columns=block_cols)
#         #tabulate_df(block_df, columns=block_cols, title="Blocks")
#         out_f = '%s.blocks.h5' % parameterObj.prefix
#         if mode == 'modify':
#             print("[#] %s of %s blocks (%s) are being retained ..." % (
#                 (len(self.blockObjs) - void_count), 
#                 (len(self.blockObjs)), 
#                 (format_percentage((len(self.blockObjs) - void_count) / len(self.blockObjs)))))
#             block_bed_df.sort_values(['sequence_id', 'bed_start'], ascending=[False, True], inplace=True)
#             block_df.sort_values(['sequence_id', 'block_start'], ascending=[False, True], inplace=True)
#             out_f = '%s.blocks.h5' % parameterObj.dataset
#         block_hdf5_store = create_hdf5_store(
#             out_f=out_f, 
#             path=parameterObj.path
#             )
#         block_bed_df.to_hdf(block_hdf5_store, 'bed', append=True)
#         block_df['distance'] = pd.to_numeric(np.where((block_df['sequence_id'] == block_df['sequence_id'].shift(-1)), block_df['block_start'].shift(-1) - block_df['block_end'], np.nan)) # compute distance to next interval)
#         plot_distance_scatter(
#             '%s.distance.png' % parameterObj.dataset,
#             "Distance between blocks (in b)", 
#             block_df['distance'].dropna().tolist()
#             )
#         barchart_y_vals, barchart_x_vals, barchart_labels, barchart_colours, barchart_populations = [], [], [], [], []
#         for idx, (sample_id, population_id) in enumerate(self.population_id_by_sample_id.items()):
#             barchart_populations.append(self.populationObjs[self.population_idx_by_population_id[population_id]].id)
#             barchart_colours.append(self.populationObjs[self.population_idx_by_population_id[population_id]].colour)
#             barchart_y_vals.append(bases_blocked_by_sample_id[sample_id])
#             barchart_x_vals.append(idx)
#             barchart_labels.append(sample_id)
#         plot_sample_barchart(
#             '%s.blocks_per_sample.png' % parameterObj.dataset,
#             "Bases in blocks by sample", 
#             barchart_y_vals, barchart_x_vals, barchart_labels, barchart_colours, barchart_populations
#             )
# #        #plot_shared_blocks(
#         #   '%s.shared_blocks_between_samples.png' % parameterObj.dataset,
#         #   "Distance between blocks (in b)", 
#         #   block_df['distance'].dropna().tolist()
#          #  )
#         block_df.to_hdf(block_hdf5_store, 'block', append=True)
#         block_hdf5_store.close()
#         # call plot from here
#         #print(bases_blocked_by_sequence_id)
#         #print(bases_blocked_by_sample_id)
#         #print(bases_blocked_by_pair_id)
#         #print(bases_blocked_by_sample_count)
#         #print(bases_blocked_by_pair_count)

# class BedObj(object):
#     __slots__ = ['sequence_id', 'start', 'end', 'pair_idxs', 'length']
    
#     def __init__(self, sequence_id, start, end, pair_idxs, length):
#         self.sequence_id = sequence_id
#         self.start = int(start)
#         self.end = int(end)
#         self.length = int(length)
#         self.pair_idxs = set(pair_idxs)

#     def __str__(self):
#         return "\t".join([self.sequence_id, str(self.start), str(self.end), str(self.length), str(self.pair_idxs)]) 

# class WindowObj(object):
#     #__slots__ = ["sequence_id", 
#     #             "start", 
#     #             "end", 
#     #             "id", 
#     #             "span", 
#     #             "centre", 
#     #             "block_ids", 
#     #             "sample_counts"
#     #             ]

#     def __init__(self, blockObjs):
#         self.sequence_id = blockObjs[0].sequence_id
#         self.start = blockObjs[0].start
#         self.end = blockObjs[-1].end
#         self.id = "%s_%s_%s" % (self.sequence_id, self.start, self.end)
#         self.span = self.end - self.start
#         self.centre = self.start + (self.span / 2)
#         self.block_ids = [blockObj.id for blockObj in blockObjs]
#         #self.variant_weights = [len(blockObj.pair_idxs) for blockObj in blockObjs]
#         self.sample_counts = [len(blockObj.sample_ids) for blockObj in blockObjs]

# class CoordinateTransformObj(object):
#     def __init__(self):
#         self.tuples_by_sequence_id = collections.defaultdict(list)

#     def add_tuple(self, sequence_id, sequence_start, sequence_end, sequence_orientation, chrom_id, chrom_start, chrom_end):
#         self.tuples_by_sequence_id[sequence_id].append((sequence_id, sequence_start, sequence_end, sequence_orientation, chrom_id, chrom_start, chrom_end))

#     def transform_coordinates(self, old_id, old_start, old_end):
#         new_id, new_start, new_end = None, None, None
#         if old_id in self.tuples_by_sequence_id:
#             for interval in self.tuples_by_sequence_id[old_id]:
#                 sequence_id, sequence_start, sequence_end, sequence_orientation, chrom_id, chrom_start, chrom_end = interval
#                 if old_start >= sequence_start and old_end <= sequence_end:
#                     new_id = chrom_id 
#                     if sequence_orientation == "+":
#                         new_start = chrom_start + (old_start - sequence_start)
#                         new_end = chrom_start + (old_end - sequence_start)
#                         return (new_id, new_start, new_end)
#                     else:
#                         new_start = chrom_start + (sequence_end - old_end)
#                         new_end = chrom_start + (sequence_end - old_start)
#                         return (new_id, new_start, new_end)
#         return (new_id, new_start, new_end)

# class BlockObj(object):

#     __slots__ = [
#         "sequence_id", 
#         "id", 
#         "pair_idxs", 
#         "sample_ids", 
#         "start", 
#         "end", 
#         "length", 
#         "span", 
#         "score", 
#         "needed", 
#         "bed_tuples", 
#         "mutype_counter_by_pair_idx", 
#         "gt_counter_by_sample_id", 
#         "void" 
#         ]

#     def __init__(self, block_id, block_length):
#         self.sequence_id = block_id.split(".")[0]
#         self.id = block_id
#         self.pair_idxs = None
#         self.sample_ids = None
#         self.start = None
#         self.end = None
#         self.length = 0
#         self.span = 0 
#         self.score = 0.0
#         self.needed = block_length
#         self.bed_tuples = [] # list of tuples (sequence_id, start, end) of consecutive regions

#         self.mutype_counter_by_pair_idx = {} # dict of counters 
#         self.gt_counter_by_sample_id = {} # dict of counters
#         self.void = False

#     def __str__(self):
#         return "[B] ID=%s %s %s %s LEN=%s SPAN=%s SCORE=%s [P]=%s\n%s\n%s" % (
#             self.id,
#             self.sequence_id, 
#             self.start, 
#             self.end, 
#             format_bases(self.length), 
#             format_bases(self.span), 
#             format_fraction(self.score, 1), 
#             self.pair_idxs, 
#             str(self.mutype_counter_by_pair_idx),
#             str(self.gt_counter_by_sample_id)
#             )

#     def __nonzero__(self):
#         if self.length:
#             return True
#         return False

#     def add_bedObj(self, bedObj, parameterObj, entityCollection):
#         '''
#         Function for adding a bedObj to the blockObj
#         [parameters]
#             - bedObj to be added
#             - parameterObj
#         [returns]
#             a) None (if bedObj has been consumed)
#             b) bedObj
#                 b1) original bedObj (if span-violation)
#                 b2) remainder bedObj (if bedObj.length > blockObj.needed)
#         [comments]
#             - span-violation:
#                 if blockObj.span > parameterObj.max_block_length:
#                     - original bedObj is returned, blockObj.score is set to 0.0
#             - blockObj.needed: allows distinction between 
#                 a) finished block: blockObj.needed = 0
#                 b) virgin block: blockObj.needed = parameterObj.block_length
#                 c) started block: 0 < blockObj.needed < parameterObj.block_length
#             - blockObj.score: 
#                 a) if blockObj.needed == parameterObj.block_length (virgin block):
#                     blockObj.score = (bedObj.pair_idxs / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
#                 b) if blockObj.needed < parameterObj.block_length (non-virgin block):
#                     blockObj.score = (len(blockObj.pair_idxs.intersection(bedObj.pair_idxs)) / parameterObj.pairs_count) * (min(bedObj.length, required_length) / parameterObj.block_length)
#                 c) if span-violation:
#                     blockObj.score = 0.0
#         '''
#         interval_length = min(self.needed, bedObj.length)
#         block_end = bedObj.start + interval_length
#         try:
#             _span = block_end - self.start # TypeError: int() argument must be a string, a bytes-like object or a number, not 'NoneType' 
#         except TypeError:
#             self.start = bedObj.start
#             _span = block_end - self.start
#         if parameterObj.max_block_length and _span > parameterObj.max_block_length:
#             self.score = 0.0
#             return bedObj
#         try:
#             self.pair_idxs = self.pair_idxs.intersection(bedObj.pair_idxs) # AttributeError: 'NoneType' object has no attribute 'intersection'
#         except AttributeError:
#             self.pair_idxs = bedObj.pair_idxs
#         self.end = block_end 
#         self.span = _span
#         self.length += interval_length
#         self.needed -= interval_length
#         self.score = (len(self.pair_idxs) / len(entityCollection.pairObjs)) * (self.length / parameterObj.block_length)
#         self.sample_ids = entityCollection.pair_idxs_to_sample_ids(self.pair_idxs)   
#         try:
#             last_tuple = self.bed_tuples.pop()
#             if (last_tuple[2] - bedObj.start) == 0: # no gap
#                 self.bed_tuples.append((bedObj.sequence_id, last_tuple[1], block_end)) 
#             else: # gap
#                 self.bed_tuples.append(last_tuple)
#                 self.bed_tuples.append((bedObj.sequence_id, bedObj.start, block_end)) 
#         except IndexError:
#             self.bed_tuples.append((bedObj.sequence_id, bedObj.start, block_end))
#         self.sequence_id = bedObj.sequence_id
#         if interval_length < bedObj.length:
#             return BedObj(bedObj.sequence_id, (bedObj.start + interval_length), bedObj.end, bedObj.pair_idxs, (bedObj.length - interval_length))    
#         else:
#             return None

# class PopulationObj(object):
#     def __init__(self, idx, population_id):
#         self.idx = idx
#         self.id = population_id
#         self.sample_ids = []
#         self.colour = COLOURS[idx]
        
#     def __str__(self):
#         return "[Population] %s %s %s" % (self.idx, self.id, ", ".join(self.sample_ids))

# class PairObj(object):
#     def __init__(self, idx, pair_id):
#         self.idx = idx
#         self.id = pair_id # tuple, use this for sample_ids since order is preserved!

#     def __str__(self):
#         return "[Pair] %s %s" % (self.idx, self.id)

# class SequenceObj(object):
#     def __init__(self, sequence_idx, sequence_id, length):
#         self.idx = sequence_idx
#         self.id = sequence_id
#         self.length = length

#     def __str__(self):
#         return "[Sequence] %s %s %s" % (self.idx, self.id, format_bases(self.length))

#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################
#########################################################################################################################

# def _write_block_bed(self, parameterObj, sample_sets='X'):
#         '''[needs fixing]
#         - use getters
#         - append in batches or https://xarray-extras.readthedocs.io/en/latest/api/csv.html
#         - or write only for specified sequences

#         '''
#         meta_seqs = self._get_meta('seqs')
#         meta_blocks = self._get_meta('blocks')
#         sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
#         blocks_count_total = sum([meta_blocks['count_by_sample_set_idx'][idx] for idx in sample_set_idxs])
#         starts = np.zeros(blocks_count_total, dtype=np.int64)
#         ends = np.zeros(blocks_count_total, dtype=np.int64)
#         # dynamically set string dtype for sequence names
#         MAX_SEQNAME_LENGTH = max([len(seq_name) for seq_name in meta_seqs['seq_names']])
#         sequences = np.zeros(blocks_count_total, dtype='<U%s' % MAX_SEQNAME_LENGTH) 
#         sample_sets = np.zeros(blocks_count_total, dtype=np.int64) 
#         variation = np.zeros((blocks_count_total, meta_seqs['mutypes_count']), dtype=np.int64)
#         missing = np.zeros(blocks_count_total, dtype=np.int64) 
#         multiallelic = np.zeros(blocks_count_total, dtype=np.int64) 
#         with tqdm(total=(len(meta_seqs['seq_names']) * len(sample_set_idxs)), desc="[%] Preparing data...", ncols=100, unit_scale=True) as pbar: 
#             offset = 0
#             for seq_name in meta_seqs['seq_names']: 
#                 for sample_set_idx in sample_set_idxs:
#                     start_key = 'blocks/%s/%s/starts' % (seq_name, sample_set_idx)
#                     end_key = 'blocks/%s/%s/ends' % (seq_name, sample_set_idx)
#                     if start_key in self.data:
#                         start_array = np.array(self.data[start_key])
#                         block_count = start_array.shape[0]
#                         starts[offset:offset+block_count] = start_array
#                         ends[offset:offset+block_count] = np.array(self.data[end_key])
#                         sequences[offset:offset+block_count] = np.full_like(block_count, seq_name, dtype='<U%s' % MAX_SEQNAME_LENGTH)
#                         sample_sets[offset:offset+block_count] = np.full_like(block_count, sample_set_idx)
#                         variation_key = 'blocks/%s/%s/variation' % (seq_name, sample_set_idx)
#                         missing_key = 'blocks/%s/%s/missing' % (seq_name, sample_set_idx)
#                         multiallelic_key = 'blocks/%s/%s/missing' % (seq_name, sample_set_idx)
#                         variation[offset:offset+block_count] = np.array(self.data[variation_key])
#                         missing[offset:offset+block_count] = np.array(self.data[missing_key]).flatten()
#                         multiallelic[offset:offset+block_count] = np.array(self.data[multiallelic_key]).flatten()
#                         offset += block_count
#                     pbar.update()
#         columns = ['sequence', 'start', 'end', 'sample_set']
#         int_bed = np.vstack([starts, ends, sample_sets, missing, multiallelic, variation.T]).T
#         mutypes_count = ["m_%s" % str(x+1) for x in range(meta_seqs['mutypes_count'])]
#         columns += ['missing', 'multiallelic'] + mutypes_count    
#         # header
#         header = ["# %s" % parameterObj._VERSION]
#         header += ["# %s = %s" % (sample_set_idx, ", ".join(meta_seqs['sample_sets'][int(sample_set_idx)])) for sample_set_idx in sample_set_idxs] 
#         header += ["# %s" % "\t".join(columns)]  
#         out_f = '%s.blocks.bed' % self.prefix
#         with open(out_f, 'w') as out_fh:
#             out_fh.write("\n".join(header) + "\n")
#         # bed
#         bed_df = pd.DataFrame(data=int_bed, columns=columns[1:])
#         bed_df['sequence'] = sequences
#         bed_df.sort_values(['sequence', 'start'], ascending=[True, True]).to_csv(out_f, na_rep='NA', mode='a', sep='\t', index=False, header=False, columns=columns, float_format='%.5f')

    # old
    #def _make_blocks(self, parameterObj, debug=False):
    #    meta_seqs = self._get_meta('seqs')
    #    meta_blocks = self._get_meta('blocks')
    #    meta_blocks['length'] = parameterObj.block_length
    #    meta_blocks['span'] = parameterObj.block_span
    #    meta_blocks['gap_run'] = parameterObj.block_gap_run
    #    meta_blocks['max_missing'] = parameterObj.block_max_missing
    #    meta_blocks['max_multiallelic'] = parameterObj.block_max_multiallelic
    #    blocks_raw_by_sample_set_idx = collections.Counter()   # all possible blocks
    #    blocks_by_sample_set_idx = collections.Counter()       # all valid blocks => only these get saved to store
    #    with tqdm(total=(len(meta_seqs['seq_names']) * len(meta_seqs['sample_sets'])), desc="[%] Building blocks ", ncols=100, unit_scale=True) as pbar:        
    #        for seq_name in meta_seqs['seq_names']:
    #            pos_key = "seqs/%s/variants/pos" % (seq_name)
    #            gt_key = "seqs/%s/variants/matrix" % (seq_name)
    #            pos = np.array(self.data[pos_key], dtype=np.int64) if pos_key in self.data else None
    #            sa_genotype_array = allel.GenotypeArray(self.data[gt_key].view(read_only=True)) if gt_key in self.data else None
    #            for sample_set_idx, sample_set in enumerate(meta_seqs['sample_sets']):
    #                start_end = self._get_interval_coordinates(seq_name=seq_name, sample_set=sample_set)
    #                if not start_end is None:
    #                    # Cut sample-set specific blocks based on intervals and block-algoritm parameters
    #                    starts, ends = start_end
    #                    #print("\n")
    #                    #print(sample_set)
    #                    #print(starts)
    #                    #print(ends)
    #                    block_sites = cut_blocks(starts, ends, meta_blocks['length'], meta_blocks['span'], meta_blocks['gap_run'])
    #                    if not block_sites is None and np.any(block_sites):
    #                        # Allocate starts/ends before overwriting position ints
    #                        block_starts = np.array(block_sites[:, 0], dtype=np.int64)
    #                        block_ends = np.array(block_sites[:, -1] + 1, dtype=np.int64)
    #                        # variants take longer than blocking
    #                        if np.any(pos) or pos is not None:
    #                            ##print('pos', pos.shape, pos)
    #                            idx_pos_in_block_sites = np.isin(pos, block_sites, assume_unique=True)
    #                            #print('idx_pos_in_block_sites', idx_pos_in_block_sites)
    #                            if np.any(idx_pos_in_block_sites):
    #                                sample_set_vcf_idxs = [meta_seqs['variants_idx_by_sample'][sample] for sample in sample_set]
    #                                idx_block_sites_in_pos = np.isin(block_sites, pos, assume_unique=True) 
    #                                sa_sample_set_genotype_array = sa_genotype_array.subset(idx_pos_in_block_sites, sample_set_vcf_idxs)
    #                                block_sites = genotype_to_mutype_array(sa_sample_set_genotype_array, idx_block_sites_in_pos, block_sites, debug)
    #                            else:
    #                                block_sites[:] = 2 # if no variants, set all to invariant    
    #                        else:
    #                            block_sites[:] = 2 # if no variants, set all to invariant
    #                        multiallelic, missing, monomorphic, variation = block_sites_to_variation_arrays(block_sites)
    #                        valid = (np.less_equal(missing, meta_blocks['max_missing']) & np.less_equal(multiallelic, meta_blocks['max_multiallelic'])).flatten()
    #                        blocks_raw_by_sample_set_idx[sample_set_idx] += valid.shape[0]
    #                        blocks_by_sample_set_idx[sample_set_idx] += valid[valid==True].shape[0]
    #                        blocks_starts_key = 'blocks/%s/%s/starts' % (seq_name, sample_set_idx)
    #                        self.data.create_dataset(blocks_starts_key, data=block_starts[valid], overwrite=True)
    #                        blocks_ends_key = 'blocks/%s/%s/ends' % (seq_name, sample_set_idx)
    #                        self.data.create_dataset(blocks_ends_key, data=block_ends[valid], overwrite=True)
    #                        blocks_variation_key = 'blocks/%s/%s/variation' % (seq_name, sample_set_idx)
    #                        self.data.create_dataset(blocks_variation_key, data=variation[valid], overwrite=True)
    #                        blocks_missing_key = 'blocks/%s/%s/missing' % (seq_name, sample_set_idx)
    #                        self.data.create_dataset(blocks_missing_key, data=missing[valid], overwrite=True)
    #                        blocks_multiallelic_key = 'blocks/%s/%s/multiallelic' % (seq_name, sample_set_idx)
    #                        self.data.create_dataset(blocks_multiallelic_key, data=multiallelic[valid], overwrite=True)
    #                pbar.update(1)
    #    meta_blocks['count_by_sample_set_idx'] = dict(blocks_by_sample_set_idx) # keys are strings
    #    meta_blocks['count_raw_by_sample_set_idx'] = dict(blocks_raw_by_sample_set_idx) # keys are strings
    #    meta_blocks['count'] = sum([count for count in blocks_by_sample_set_idx.values()])
    
#### old code below

    # def _make_windows(self, parameterObj, sample_sets='X'):
    #     meta_seqs = self._get_meta('seqs')
    #     meta_windows = self._get_meta('windows')
    #     meta_windows['size'] = parameterObj.window_size
    #     meta_windows['step'] = parameterObj.window_step
    #     meta_windows['count'] = 0
    #     sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
    #     with tqdm(meta_seqs['seq_names'], total=(len(meta_seqs['seq_names']) * len(sample_set_idxs)), desc="[%] Constructing windows ", ncols=100, unit_scale=True) as pbar: 
    #         for seq_name in meta_seqs['seq_names']:
    #             variation, starts, ends = [], [], []
    #             for sample_set_idx in sample_set_idxs:
    #                 variation_key = 'blocks/%s/%s/variation' % (seq_name, sample_set_idx)
    #                 if variation_key in self.data:
    #                     variation.append(np.array(self.data[variation_key]))
    #                     start_key = 'blocks/%s/%s/starts' % (seq_name, sample_set_idx)
    #                     starts.append(np.array(self.data[start_key]))
    #                     end_key = 'blocks/%s/%s/ends' % (seq_name, sample_set_idx)
    #                     ends.append(np.array(self.data[end_key]))
    #                 pbar.update()
    #             variation_array = np.concatenate(variation, axis=0)
    #             start_array = np.concatenate(starts, axis=0)
    #             end_array = np.concatenate(ends, axis=0)
    #             # window_variation : shape = (windows, blocklength, 4)
    #             window_variation, window_starts, window_ends, window_pos_mean, window_pos_median = cut_windows(variation_array, sample_set_idxs, start_array, end_array, num_blocks=parameterObj.window_size, num_steps=parameterObj.window_step)
    #             #b, counts = np.unique(variation, return_counts=True, axis=0)
    #             self.data.create_dataset("windows/%s/variation" % seq_name, data=window_variation, overwrite=True)
    #             self.data.create_dataset("windows/%s/starts" % seq_name, data=window_starts, overwrite=True)
    #             self.data.create_dataset("windows/%s/ends" % seq_name, data=window_ends, overwrite=True)
    #             self.data.create_dataset("windows/%s/pos_mean" % seq_name, data=window_pos_mean, overwrite=True)
    #             self.data.create_dataset("windows/%s/pos_median" % seq_name, data=window_pos_median, overwrite=True)
    #             meta_windows['count'] += window_variation.shape[0]

        #window_info_rows = []
        #window_mutuple_tally = []
        ## window bsfs
        #for seq_name in tqdm(meta['seq_names'], total=len(meta['seq_names']), desc="[%] Generating output ", ncols=100):
        #    variations = self.data["seqs/%s/windows/variation" % seq_name]
        #    print('variations', variations.shape, variations.info, variations.hexdigest())
        #    for window_id, variation, midpoint_mean, midpoint_median in zip(window_ids, variations, midpoint_means, midpoint_medians):
        #        pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (meta['blocks_length'] * variation.shape[0]))
        #        blocks = variation.shape[0]
        #        window_info_rows.append([window_id, seq_name, midpoint_mean, midpoint_median, pi_1, pi_2, d_xy, f_st, fgv, blocks])
        #        # bsfs for each window
        #        mutuple, counts = np.unique(variation, return_counts=True, axis=0)
        #        tally = np.concatenate([counts[:, np.newaxis], mutuple], axis =-1)
        #        windows = np.array([window_id] * tally.shape[0])
        #        window_mutuple_tally.append(np.concatenate([windows[:, np.newaxis], tally], axis =-1))
        #window_bsfs_cols = ['window_id', 'count'] + [x+1 for x in range(meta['mutypes_count'])]
        #print(window_bsfs_cols)
        #print(window_mutuple_tally)
        #window_bsfs_df = pd.DataFrame(np.vstack(window_mutuple_tally), columns=window_bsfs_cols)
        #print("[+] Made %s windows" % window_bsfs_df['window_id'].nunique()) 
        #window_bsfs_f = "%s.window_bsfs.tsv" % self.prefix
        #window_bsfs_df.to_csv(window_bsfs_f, sep='\t', index=False)
        #print("[>] Created: %r" % str(window_bsfs_f))
        #window_info_cols = ['window_id', 'seq_name', 'midpoint_mean', 'midpoint_median', 'pi_%s' % meta['population_ids'][0], 'pi_%s' % meta['population_ids'][1], 'd_xy', 'f_st', 'f']
        #window_info_df = pd.DataFrame(window_info_rows, columns=window_info_cols)
        #window_info_f = "%s.window_info.tsv" % self.prefix
        #window_info_df.to_csv(window_info_f, sep='\t', index=False)
        #print("[>] Created: %r" % str(window_info_f))
        #self.plot_fst_genome_scan(window_info_df)
        #self.plot_pi_genome_scan(window_info_df)
        ##     plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)


    #def _make_blocks_threaded(self, parameterObj, debug=False, threads=2):
    #    '''there might be some speed improvement here ... has to be finished...'''
    #    meta = self.data['seqs'].attrs
    #    meta['blocks_length'] = parameterObj.block_length
    #    meta['blocks_span'] = parameterObj.block_span
    #    meta['blocks_gap_run'] = parameterObj.block_gap_run
    #    meta['blocks_max_missing'] = parameterObj.block_max_missing
    #    meta['blocks_max_multiallelic'] = parameterObj.block_max_multiallelic
    #    blocks_raw_by_sample_set_idx = collections.Counter()   # all possible blocks
    #    blocks_by_sample_set_idx = collections.Counter()       # all valid blocks => only these get saved to store
    #    params = [(seqs, str(sample_seq_idx)) for seqs, sample_seq_idx in itertools.product(meta['seq_names'], range(0, len(meta['sample_sets'])))]
    #    print(params)
    #    results = []
    #    with poolcontext(processes=threads) as pool:
    #        with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
    #            for blockObjs in pool.imap_unordered(block_algorithm, params):
    #                results.append(blockObjs)
    #                pbar.update()
#
    #    with tqdm(total=(len(meta['seq_names']) * len(meta['sample_sets'])), desc="[%] Building blocks ", ncols=100, unit_scale=True) as pbar:        
    #        for seq_name in meta['seq_names']:        
    #            pos_key = "seqs/%s/variants/pos" % (seq_name)
    #            gt_key = "seqs/%s/variants/matrix" % (seq_name)
    #            for sample_set_idx, (sample_set, sample_set_cartesian) in enumerate(zip(meta['sample_sets'], meta['sample_sets_inter'])):
    #                params.append(seq_name, sample_set_idx, pos_key, gt_key)
    #                pbar.update(1)
    #    meta['blocks_by_sample_set_idx'] = dict(blocks_by_sample_set_idx) # keys are strings
    #    meta['blocks_raw_by_sample_set_idx'] = dict(blocks_raw_by_sample_set_idx) # keys are strings

    #def plot_bsfs_pcp(self, sample_set='X'):
    #    '''plots a bsfs pcp for a given sample_set. 
    #    '''
    #    # https://stackoverflow.com/questions/8230638/parallel-coordinates-plot-in-matplotlib
    #    # http://www.shengdongzhao.com/publication/tracing-tuples-across-dimensions-a-comparison-of-scatterplots-and-parallel-coordinate-plots/
    #    #meta = self.data['seqs'].attrs
    #    bsfs = bsfs_to_2d(self._get_block_bsfs(sample_sets='X'))
    #    print(bsfs)
    #    freq = bsfs[:,0] / np.sum(bsfs[:,0])
    #    meta = self.data['seqs'].attrs
    #    x = ['m_1', 'm_2', 'm_3', 'm_4']
    #    bins = np.linspace(0, 1, 9)
    #    freq = bins[np.digitize(freq, bins)]
    #    data = bsfs[:,1:] 
    #    cmap = plt.get_cmap("Greys")
    #    print('data', data.shape, data)
    #    print('freq', freq.shape, freq)
    #    fig, axes = plt.subplots(1, 3, sharey=False, figsize=(15,5))
    #    axes[0].plot(x,data.T, c=cmap(freq))
    #    axes[1].plot(x,data.T, c=cmap(freq))
    #    axes[2].plot(x,data.T, c=cmap(freq))
    #    plt.subplots_adjust(wspace=0)
    #    #     min_max_range[mutype] = [np.min(, df[col].max(), np.ptp(df[col])]
    #    #         df[col] = np.true_divide(df[col] - df[col].min(), np.ptp(df[col]))
    #    # ynames = ['m_%s' for idx in range(1, meta['mutypes_count'] + 1)]
    #    # import pandas as pd
    #    # from pandas.tools.plotting import parallel_coordinates
    #    # ax = pd.tools.plotting.parallel_coordinates(mutypes)
#
    #    # #ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
    #    # #ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
    #    # #y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
    #    # #ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
    #    # #ax.set_ylim(y_lim)
    #    # #ax.spines['right'].set_visible(False)
    #    # #ax.spines['top'].set_visible(False)
    #    # #ax.legend(numpoints=1)
    #    # #plt.ylabel('Pi')
    #    # #plt.xlabel("Genome coordinate")
    #    out_f = '%s.pcp.png' % self.prefix
    #    plt.savefig(out_f, format="png")
        #print("[>] Created: %r" % str(out_f))
        #plt.close(fig)

        #fig, host = plt.subplots()
        #
        ## create some dummy data
        #ynames = ['P1', 'P2', 'P3', 'P4', 'P5']
        #N1, N2, N3 = 10, 5, 8
        #N = N1 + N2 + N3
        #category = np.concatenate([np.full(N1, 1), np.full(N2, 2), np.full(N3, 3)])
        #y1 = np.random.uniform(0, 10, N) + 7 * category
        #y2 = np.sin(np.random.uniform(0, np.pi, N)) ** category
        #y3 = np.random.binomial(300, 1 - category / 10, N)
        #y4 = np.random.binomial(200, (category / 6) ** 1/3, N)
        #y5 = np.random.uniform(0, 800, N)
        #
        ## organize the data
        #ys = np.dstack([y1, y2, y3, y4, y5])[0]
        #ymins = ys.min(axis=0)
        #ymaxs = ys.max(axis=0)
        #dys = ymaxs - ymins
        #ymins -= dys * 0.05  # add 5% padding below and above
        #ymaxs += dys * 0.05
        #dys = ymaxs - ymins
        #
        ## transform all data to be compatible with the main axis
        #zs = np.zeros_like(ys)
        #zs[:, 0] = ys[:, 0]
        #zs[:, 1:] = (ys[:, 1:] - ymins[1:]) / dys[1:] * dys[0] + ymins[0]
        #
        #
        #axes = [host] + [host.twinx() for i in range(ys.shape[1] - 1)]
        #for i, ax in enumerate(axes):
        #    ax.set_ylim(ymins[i], ymaxs[i])
        #    ax.spines['top'].set_visible(False)
        #    ax.spines['bottom'].set_visible(False)
        #    if ax != host:
        #        ax.spines['left'].set_visible(False)
        #        ax.yaxis.set_ticks_position('right')
        #        ax.spines["right"].set_position(("axes", i / (ys.shape[1] - 1)))
        #
        #host.set_xlim(0, ys.shape[1] - 1)
        #host.set_xticks(range(ys.shape[1]))
        #host.set_xticklabels(ynames, fontsize=14)
        #host.tick_params(axis='x', which='major', pad=7)
        #host.spines['right'].set_visible(False)
        #host.xaxis.tick_top()
        #host.set_title('Parallel Coordinates Plot', fontsize=18)
        #
        #colors = plt.cm.tab10.colors
        #for j in range(N):
        #    # to just draw straight lines between the axes:
        #    # host.plot(range(ys.shape[1]), zs[j,:], c=colors[(category[j] - 1) % len(colors) ])
        #
        #    # create bezier curves
        #    # for each axis, there will a control vertex at the point itself, one at 1/3rd towards the previous and one
        #    #   at one third towards the next axis; the first and last axis have one less control vertex
        #    # x-coordinate of the control vertices: at each integer (for the axes) and two inbetween
        #    # y-coordinate: repeat every point three times, except the first and last only twice
        #    verts = list(zip([x for x in np.linspace(0, len(ys) - 1, len(ys) * 3 - 2, endpoint=True)],
        #                     np.repeat(zs[j, :], 3)[1:-1]))
        #    # for x,y in verts: host.plot(x, y, 'go') # to show the control points of the beziers
        #    codes = [Path.MOVETO] + [Path.CURVE4 for _ in range(len(verts) - 1)]
        #    path = Path(verts, codes)
        #    patch = patches.PathPatch(path, facecolor='none', lw=1, edgecolor=colors[category[j] - 1])
        #    host.add_patch(patch)
        #plt.tight_layout()
        #plt.show()

    
#     def plot_pi_genome_scan(self, window_df):
#         offset_by_sequence_id = {}
#         offset = 0
#         x_boundaries = []
#         for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
#             offset_by_sequence_id[sequence_id] = offset
#             x_boundaries.append(offset)
#             offset += sequence_length
#         x_boundaries.append(offset)
#         #print([(sequenceObj.id, sequenceObj.length) for sequenceObj in sequenceObjs])
#         #print(x_boundaries)
#         fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
#         #connecting dots
#         ax = fig.add_subplot(111)  
#         window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
#         window_df.sort_values(['rel_pos'], inplace=True)
#         #print(window_df)
#         pi_A_key = list(window_df.columns)[4]
#         pi_B_key = list(window_df.columns)[5]
#         ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
#         ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
#         y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
#         ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
#         ax.set_ylim(y_lim)
#         ax.spines['right'].set_visible(False)
#         ax.spines['top'].set_visible(False)
#         ax.legend(numpoints=1)
#         plt.ylabel('Pi')
#         plt.xlabel("Genome coordinate")
#         out_f = '%s.pi_genome_scan.png' % self.prefix
#         plt.tight_layout()
#         fig.savefig(out_f, format="png")
#         print("[>] Created: %r" % str(out_f))
#         plt.close(fig)

#     def plot_fst_genome_scan(self, window_df):
#         offset_by_sequence_id = {}
#         offset = 0
#         x_boundaries = []
#         for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
#             offset_by_sequence_id[sequence_id] = offset
#             x_boundaries.append(offset)
#             offset += sequence_length
#         x_boundaries.append(offset)
#         fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
#         #connecting dots
#         ax = fig.add_subplot(111)  
#         y_lim = (0.0, 1.0)
#         window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
#         window_df.sort_values(['rel_pos'], inplace=True)
#         ax.plot(window_df['rel_pos'], window_df['f_st'], color='lightgrey', alpha=0.8, linestyle='-', linewidth=1)
#         scatter = ax.scatter(window_df['rel_pos'], window_df['f_st'], c=window_df['d_xy'], alpha=1.0, cmap='PiYG_r', edgecolors='white', marker='o', s=40, linewidth=0.2)
#         cbar = fig.colorbar(scatter, ax=ax)
#         cbar.ax.set_title('D_xy')
#         ax.vlines(x_boundaries, 0.0, 1.0, colors=['lightgrey'], linestyles='dashed', linewidth=1)
#         ax.set_ylim(y_lim)
#         ax.spines['right'].set_visible(False)
#         ax.spines['top'].set_visible(False)
#         plt.ylabel('F_st')
#         plt.xlabel("Genome coordinate")
#         ax.autoscale_view(tight=None, scalex=True, scaley=True)
#         out_f = '%s.fst_genome_scan.png' % self.prefix
#         fig.savefig(out_f, format="png")
#         plt.close(fig)
#         print("[>] Created: %r" % str(out_f))

    #def dump_blocks(self, parameterObj, cartesian_only=True):
    #    meta = self.data['seqs'].attrs
    #    sample_set_idxs = [idx for (idx, is_cartesian) in enumerate(meta['sample_sets_inter']) if is_cartesian] if cartesian_only else range(len(meta['sample_sets']))
    #    variation_global = []
    #    with tqdm(total=(len(meta['seq_names']) * len(sample_set_idxs)), desc="[%] Writing bSFSs ", ncols=100, unit_scale=True) as pbar: 
    #        for seq_name in meta['seq_names']: 
    #            for sample_set_idx in sample_set_idxs:
    #                variation_key = 'seqs/%s/blocks/%s/variation' % (seq_name, sample_set_idx)
    #                variation_global.append(np.array(self.data[variation_key]))#[valid]
    #                pbar.update()
    #    variation_global_array = np.concatenate(variation_global, axis=0)
    #    # popgen
    #    variation_global = []
        #metrics_rows = []
        # is order (pi_1, pi_2, d_xy, f_st, fgv) correct?
        # for sample_set_idx in data_by_key_by_sample_set_idx:
        #     sample_set_ids = self.data.attrs['sample_sets'][sample_set_idx]
        #     #print(data_by_key_by_sample_set_idx)
        #     block_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites'], axis=0))
        #     interval_sites = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['interval_sites'], axis=0))
        #     block_sites_valid = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['block_sites_valid'], axis=0))
        #     variation_array = np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['variation'], axis=0)
        #     missing_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['missing'], axis=0))
        #     multiallelic_count = np.sum(np.concatenate(data_by_key_by_sample_set_idx[sample_set_idx]['multiallelic'], axis=0))
        #     hetB_count, hetA_count, hetAB_count, fixed_count = np.sum(variation_array, axis=0)
        #     #pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, (self.data.attrs['block_length'] * variation_array.shape[0]))    
        #     pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_array, block_sites_valid)    
        #     metrics_rows.append([
        #         sample_set_ids[0], 
        #         sample_set_ids[1],     
        #         block_sites,
        #         interval_sites,
        #         block_sites_valid,
        #         np.divide(block_sites_valid, self.data.attrs['block_length']),
        #         fgv,
        #         missing_count,
        #         multiallelic_count,
        #         hetA_count, 
        #         hetB_count, 
        #         hetAB_count, 
        #         fixed_count,
        #         pi_1, 
        #         pi_2, 
        #         d_xy, 
        #         f_st
        #         ])
        # # output metrics 
        # header = [
        #     self.data.attrs['pop_ids'][0], 
        #     self.data.attrs['pop_ids'][1], 
        #     'block_sites', 
        #     'interval_sites', 
        #     'block_sites_valid', 
        #     'blocks', 
        #     'fgv', 
        #     'missing', 
        #     'multiallelic', 
        #     'hetA', 
        #     'hetB', 
        #     'hetAB', 
        #     'fixed', 
        #     'piA', 
        #     'piB', 
        #     'dxy', 
        #     'fst'
        #     ]
        # pd.DataFrame(data=metrics_rows, columns=header, dtype='int64').to_hdf("%s.block_stats.h5" % self.prefix, 'bsfs', format='table')

        #pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation_global_array, (self.data.attrs['block_length'] * variation_global_array.shape[0]))
        #print("[+] Pi_%s = %s; Pi_%s = %s; D_xy = %s; F_st = %s; FGVs = %s / %s blocks (%s)" % (self.data.attrs['pop_ids'][0], pi_1, self.data.attrs['pop_ids'][1], pi_2, d_xy, f_st, fgv, variation_global_array.shape[0], format_percentage(fgv / variation_global_array.shape[0]))) 


#    def test_dask(self, grids=None, data=None):
#        if grids is None or data is None:
#            raise ValueError('gridsearch: needs grid and data') 
#        from dask.distributed import LocalCluster, Client
#        cluster= LocalCluster(
#            n_workers=8, 
#            threads_per_worker=1,
#            dashboard_address='localhost:8787', 
#            interface='lo0', 
#            **{'local_directory': 'dasktemp', 'memory_limit': '2G'})
#        client = Client(cluster)
#        array = dask.array.ones((100000,100000), chunks=(10000,10000))
#        a = client.submit(dask.array.mean, array).result()
#        return a.compute()
        
#   def gridsearch(self, grids=None, data=None):
#       '''returns 2d array of likelihoods of shape (windows, grids)'''
#       if grids is None or data is None:
#           raise ValueError('gridsearch: needs grid and data') 
#       from dask.distributed import Client, LocalCluster
#       cluster= LocalCluster(
#           n_workers=8, 
#           threads_per_worker=1,
#           dashboard_address='localhost:8787', 
#           interface='lo0', 
#           **{'local_directory': 'dasktemp', 'memory_limit': '3G'})
#       client = Client(cluster)
#       data_array = dask.array.from_array(data, chunks=(10000, *data.shape[1:]))
#       print('# data_array', type(data_array))
#       grid_array = dask.array.from_array(grids, chunks=(1, *grids.shape[1:]))
#       print('# grid_array', type(grid_array))
#       grids_masked = dask.array.ma.masked_array(grids, grids==0, fill_value=np.nan)
#       grids_log_temp = client.submit(dask.array.log, grids_masked).result()

#       print('# grids_log_temp', type(grids_log_temp))
#       grids_log = client.submit(dask.array.nan_to_num, grids_log_temp).result().compute()
#       print('grids_log', type(grids_log))
#       print('done')
#       # print('data', type(data))
#       # print('grids', type(grids))
#       # data = dask.array.from_array(data, chunks=(500, *data.shape[1:]))
#       # # print('data', type(data), data.chunks)
#       # grids = dask.array.from_array(grids, chunks=(10, *grids.shape[1:]))
#       # # print('grids', type(grids), grids.chunks)
#       # #np.log(grids, where=grids>0, out=grids_log)
#       # grids[grids==0] = np.nan
#       # 
#       # grids_log = client.submit(dask.array.log, grids).result().compute()
#       # print('grids_log', type(grids_log))
#       # grids_log[grids_log==np.nan] = 0
#       # data_scattered = client.scatter(data[:, None])
#       # grids_log_scattered = client.scatter(grids_log)
#       # m_res = client.submit(dask.array.multiply, data_scattered, grids_log_scattered)
#       # #m_res.rechunk({0:100, 1: 4, 2: 4, 3: 4, 4: 4})
#       # #res_scattered = client.scatter(m_res)
#       # #r_res = client.submit(dask.array.apply_over_axes, dask.array.sum, res_scattered, axes=[2,3,4,5])
#       # x_res = client.submit(dask.array.apply_over_axes, dask.array.sum, m_res, axes=[2,3,4,5])
#       # #m1_res.rechunk(100000)
#       # #r_res_scattered = client.scatter(r_res)
#       # y_res = client.submit(dask.array.squeeze, x_res)
#       # y_res.result().compute()
#       # #res.visualize()
#       # print('y_res', type(y_res))
#       # return y_res
#       #return a
#       return True

    #def _process_config(self):
#
    #    if self._MODULE in ['makegrid', 'inference', 'simulate', 'gridsearch']:
    #        self.config['mu']['blockslength'] = self._get_blocks_length(self.zstore)
    #        self.config['parameters']['mu'] = self.config['mu']['mu']
    #        self._expand_params()
    #        self.reference, self.toBeSynced = self._get_pops_to_sync()
    #        self._remove_pop_from_dict(self.toBeSynced)
    #        self.parameter_combinations = self._dict_product()
    #        self._sync_pop_sizes(self.reference, self.toBeSynced)
    #    elif self._MODULE in ['optimize', 'optimize']:
    #        #TO BE CHECKED: which bits are we still using
    #        #determine parameters that are fixed:
    #        self.fixed_params = self._get_fixed_params()
    #        #self.config['mu']['blockslength'] = self._get_blocks_length(self.zstore)
    #        #self.config['parameters']['mu'] = self.config['mu']['mu']
    #        reference_pop=self.config['populations']['reference_pop']
    #        #syncing pop sizes
    #        self.reference, self.toBeSynced = self._get_pops_to_sync()
    #        if self.toBeSynced:
    #            if reference_pop in self.toBeSynced:
    #                sys.exit(f"[X] Set reference pop to {self.reference}.")
    #        toBeSynced_pops = [f'Ne_{s}' for s in self.toBeSynced] if self.toBeSynced!=None else []
    #        self.fixed_params = [pop for pop in self.fixed_params if pop not in toBeSynced_pops]
    #        #verify if any Ne fixed, whether one of those Ne is self.reference
    #        fixed_Nes = [p for p in self.fixed_params if p.startswith('Ne')]
    #        if len(fixed_Nes)>0:
    #            if not f"Ne_{reference_pop}" in fixed_Nes:
    #                sys.exit("[X] No. No. No. It would make much more sense to set a population with a fixed size as reference.")
    #        #self._sync_pop_sizes_optimize(self.reference, self.toBeSynced)
    #        self.parameter_combinations = self._return_boundaries()
    #    else:
    #        sys.exit("[X] gimble.py_processing_config: Not implemented yet.")

# def calculate_popgen_from_array(mutype_array, sites):
#     # print('# Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
#     pi_1 = float("%.8f" % np.divide(np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), sites)) # average heterozygosity
#     pi_2 = float("%.8f" % np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 2]), sites)) # average heterozygosity
#     d_xy = float("%.8f" % np.divide(np.divide(np.sum(mutype_array[:, 0]) + np.sum(mutype_array[:, 1]) + np.sum(mutype_array[:, 2]), 2.0) + np.sum(mutype_array[:, 3]), sites))
#     mean_pi = (pi_1 + pi_2) / 2.0
#     total_pi = (d_xy + mean_pi) / 2.0 # special case of pairwise Fst
#     f_st = np.nan
#     if (total_pi):
#         f_st = float("%.8f" % ((total_pi - mean_pi) / total_pi)) # special case of pairwise Fst
#     fgv = len(mutype_array[(mutype_array[:, 2] > 0) & (mutype_array[:, 3] > 0)])
#     return (pi_1, pi_2, d_xy, f_st, fgv)

    #def dump_windows(self, parameterObj):
    #    window_info_rows = []
    #    window_mutuple_tally = []
    #    for sequence_id in tqdm(self.data.attrs['sequence_ids'], total=len(self.data.attrs['sequence_ids']), desc="[%] Generating output ", ncols=100):
    #        variations = self.data["%s/windows/variation" % sequence_id]
    #        #print(self.data["%s/windows/starts" % sequence_id][:])
    #        #print(self.data["%s/windows/pos_mean" % sequence_id][:])
    #        #print(self.data["%s/windows/pos_median" % sequence_id][:])
    #        window_ids = np.array(["_".join([sequence_id, _start, _end]) for (_start, _end) in zip(
    #            np.array(self.data["%s/windows/starts" % sequence_id]).astype(str), 
    #            np.array(self.data["%s/windows/ends" % sequence_id]).astype(str))])
    #        #window_ids = self.data["%s/windows/window_id" % sequence_id]
    #        midpoint_means = self.data["%s/windows/pos_mean" % sequence_id]
    #        midpoint_medians = self.data["%s/windows/pos_median" % sequence_id]
    #        for window_id, variation, midpoint_mean, midpoint_median in zip(window_ids, variations, midpoint_means, midpoint_medians):
    #            pi_1, pi_2, d_xy, f_st, fgv = calculate_popgen_from_array(variation, (self.data.attrs['block_length'] * variation.shape[0]))
    #            window_info_rows.append([window_id, sequence_id, midpoint_mean, midpoint_median, pi_1, pi_2, d_xy, f_st, fgv/variation.shape[0]])
    #            # mutuple barchart
    #            mutypes, counts = np.unique(variation, return_counts=True, axis=0)
    #            tally = np.concatenate([counts[:, np.newaxis], mutypes], axis =-1)
    #            windows = np.array([window_id] * tally.shape[0])
    #            window_mutuple_tally.append(np.concatenate([windows[:, np.newaxis], tally], axis =-1))
    #    window_bsfs_cols = ['window_id', 'count'] + [x+1 for x in range(self.data.attrs['mutypes_count'])]
    #    window_bsfs_df = pd.DataFrame(np.vstack(window_mutuple_tally), columns=window_bsfs_cols)
    #    print("[+] Made %s windows" % window_bsfs_df['window_id'].nunique()) 
    #    window_bsfs_f = "%s.window_bsfs.tsv" % self.prefix
    #    window_bsfs_df.to_csv(window_bsfs_f, sep='\t', index=False)
    #    print("[>] Created: %r" % str(window_bsfs_f))
    #
    #    window_info_cols = ['window_id', 'sequence_id', 'midpoint_mean', 'midpoint_median', 'pi_%s' % self.data.attrs['pop_ids'][0], 'pi_%s' % self.data.attrs['pop_ids'][1], 'd_xy', 'f_st', 'fgv']
    #    window_info_df = pd.DataFrame(window_info_rows, columns=window_info_cols)
    #    window_info_f = "%s.window_info.tsv" % self.prefix
    #    window_info_df.to_csv(window_info_f, sep='\t', index=False)
    #    print("[>] Created: %r" % str(window_info_f))        
    #    self.plot_fst_genome_scan(window_info_df)
    #    self.plot_pi_genome_scan(window_info_df)
    #    #     plot_pi_scatter(window_df, '%s.pi_scatter.png' % parameterObj.dataset)
    
    #def plot_pi_genome_scan(self, window_df):
    #    offset_by_sequence_id = {}
    #    offset = 0
    #    x_boundaries = []
    #    for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
    #        offset_by_sequence_id[sequence_id] = offset
    #        x_boundaries.append(offset)
    #        offset += sequence_length
    #    x_boundaries.append(offset)
    #    #print([(sequenceObj.id, sequenceObj.length) for sequenceObj in sequenceObjs])
    #    #print(x_boundaries)
    #    fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
    #    #connecting dots
    #    ax = fig.add_subplot(111)  
    #    window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
    #    window_df.sort_values(['rel_pos'], inplace=True)
    #    #print(window_df)
    #    pi_A_key = list(window_df.columns)[4]
    #    pi_B_key = list(window_df.columns)[5]
    #    ax.plot(window_df['rel_pos'], window_df[pi_A_key], color='orange', alpha=0.8, linestyle='-', linewidth=1, label=pi_A_key.replace('pi_', ''))
    #    ax.plot(window_df['rel_pos'], window_df[pi_B_key], color='dodgerblue', alpha=0.8, linestyle='-', linewidth=1, label=pi_B_key.replace('pi_', ''))
    #    y_lim = (min(window_df[pi_A_key].min(), window_df[pi_B_key].min()), max(window_df[pi_A_key].max(), window_df[pi_B_key].max()))
    #    ax.vlines(x_boundaries, y_lim[0], y_lim[1], colors=['lightgrey'], linestyles='dashed', linewidth=1)
    #    ax.set_ylim(y_lim)
    #    ax.spines['right'].set_visible(False)
    #    ax.spines['top'].set_visible(False)
    #    ax.legend(numpoints=1)
    #    plt.ylabel('Pi')
    #    plt.xlabel("Genome coordinate")
    #    out_f = '%s.pi_genome_scan.png' % self.prefix
    #    plt.tight_layout()
    #    fig.savefig(out_f, format="png")
    #    print("[>] Created: %r" % str(out_f))
    #    plt.close(fig)
    #
    #def plot_fst_genome_scan(self, window_df):
    #    offset_by_sequence_id = {}
    #    offset = 0
    #    x_boundaries = []
    #    for sequence_id, sequence_length in zip(self.data.attrs['sequence_ids'], self.data.attrs['sequence_length']):
    #        offset_by_sequence_id[sequence_id] = offset
    #        x_boundaries.append(offset)
    #        offset += sequence_length
    #    x_boundaries.append(offset)
    #    fig = plt.figure(figsize=(18,4), dpi=200, frameon=True)
    #    #connecting dots
    #    ax = fig.add_subplot(111)  
    #    y_lim = (0.0, 1.0)
    #    window_df['rel_pos'] = window_df['midpoint_median'] + window_df['sequence_id'].map(offset_by_sequence_id)
    #    window_df.sort_values(['rel_pos'], inplace=True)
    #    ax.plot(window_df['rel_pos'], window_df['f_st'], color='lightgrey', alpha=0.8, linestyle='-', linewidth=1)
    #    scatter = ax.scatter(window_df['rel_pos'], window_df['f_st'], c=window_df['d_xy'], alpha=1.0, cmap='PiYG_r', edgecolors='white', marker='o', s=40, linewidth=0.2)
    #    cbar = fig.colorbar(scatter, ax=ax)
    #    cbar.ax.set_title('D_xy')
    #    ax.vlines(x_boundaries, 0.0, 1.0, colors=['lightgrey'], linestyles='dashed', linewidth=1)
    #    ax.set_ylim(y_lim)
    #    ax.spines['right'].set_visible(False)
    #    ax.spines['top'].set_visible(False)
    #    plt.ylabel('F_st')
    #    plt.xlabel("Genome coordinate")
    #    ax.autoscale_view(tight=None, scalex=True, scaley=True)
    #    out_f = '%s.fst_genome_scan.png' % self.prefix
    #    fig.savefig(out_f, format="png")
    #    plt.close(fig)
    #    print("[>] Created: %r" % str(out_f))

### OLD CODE 

# old
#def genotype_to_mutype_array(sa_genotype_array, idx_block_sites_in_pos, block_sites, debug=True):
#    # DEPRECATED
#    warnings.warn("lib.gimble.genotype_to_mutype_array() is deprecated. Use gt2fmac() ...", DeprecationWarning)
#    '''
#    - possible errors:
#        - if whole sequence has only missing GTs in genotypes for a sample set, then np_allele_count_array will be empty ... (should never happen)
#    '''
#    np_genotype_array = np.array(sa_genotype_array)
#    #print('np_genotype_array', np_genotype_array.shape, np_genotype_array)
#    np_allele_count_array = np.ma.masked_equal(sa_genotype_array.count_alleles(), 0, copy=True) 
#    #print('np_allele_count_array', np_allele_count_array.shape, np_allele_count_array)
#    #print('np.any(np_allele_count_array)', np.any(np_allele_count_array))
#    allele_map = np.ones((np_allele_count_array.shape), dtype=np.int64) * np.arange(np_allele_count_array.shape[-1], dtype=np.int64)
#    #print('allele_map', allele_map.shape, allele_map)
#    #print('np.any(allele_map)', np.any(allele_map))
#    if np.any(np_allele_count_array) and np.any(allele_map):
#        idx_max_global_allele_count = np.nanargmax(np_allele_count_array, axis=1)
#        idx_min_global_allele_count = np.nanargmin(np_allele_count_array, axis=1)
#        has_major_allele = (idx_max_global_allele_count != idx_min_global_allele_count)
#        idx_min_prime_allele = np.amin(np_genotype_array[:,0], axis=1)
#        idx_min_global_allele = np.amin(np.amin(np_genotype_array, axis=1), axis=1)
#        idx_max_global_allele = np.amax(np.amax(np_genotype_array, axis=1), axis=1)
#        idx_major_allele = np.where(
#            has_major_allele, 
#            idx_max_global_allele_count, 
#            idx_min_prime_allele)
#        idx_minor_allele = np.where(
#            has_major_allele, 
#            idx_min_global_allele_count, 
#            np.where((
#                idx_min_global_allele == idx_min_prime_allele),
#                np.max((idx_min_global_allele, idx_max_global_allele), axis=0), 
#                np.min((idx_min_global_allele, idx_max_global_allele), axis=0)))
#        # for each genotype (np.arange(allele_map.shape[0])), set minor allele to 1 (1st do minor, so that overwritten if monomorphic)
#        allele_map[np.arange(allele_map.shape[0]), idx_minor_allele] = 1 
#        # for each genotype (np.arange(allele_map.shape[0])), set major allele to 0
#        allele_map[np.arange(allele_map.shape[0]), idx_major_allele] = 0
#    folded_minor_allele_counts = sa_genotype_array.map_alleles(allele_map).to_n_alt(fill=-1)
#    #print(folded_minor_allele_counts)
#    folded_minor_allele_counts[np.any(sa_genotype_array.is_missing(), axis=1)] = np.ones(2) * -1        # -1, -1 for missing => -1
#    folded_minor_allele_counts[(np_allele_count_array.count(axis=1) > 2)] = np.ones(2) * (-1, -2)       # -1, -2 for multiallelic => -2
#    block_sites[idx_block_sites_in_pos] = szudzik_pairing(folded_minor_allele_counts) + 2               # add 2 so that not negative for bincount
#    block_sites[~idx_block_sites_in_pos] = 2                                                            # monomorphic = 2 (0 = multiallelic, 1 = missing)
#    if debug == True:
#        print("# block_sites as mutype array", block_sites)
#    if debug == True:
#        block_sites_pos = block_sites.flatten()
#        pos_df = pd.DataFrame(block_sites_pos[idx_block_sites_in_pos.flatten()], dtype='int64', columns=['pos'])
#        genotypes_df = pd.DataFrame(np_genotype_array.reshape(np_genotype_array.shape[0], 4), dtype='i4', columns=['a1', 'a2', 'b1', 'b2'])        
#        block_sites_df = pos_df.join(genotypes_df)
#        folded_minor_allele_count_df = pd.DataFrame(folded_minor_allele_counts, dtype='int8', columns=['fmAC_a', 'fmAC_b'])
#        block_sites_df = block_sites_df.join(folded_minor_allele_count_df)
#        variants = pd.DataFrame(block_sites[idx_block_sites_in_pos], dtype='int', columns=['SVar'])
#        block_sites_df = block_sites_df.join(variants)
#        print('# Mutypes: 0=MULTI, 1=MISS, 2=MONO, 3=HetB, 4=HetA, 5=HetAB, 6=Fixed')
#        print(block_sites_df)
#    return block_sites

# old
#def cut_blocks(interval_starts, interval_ends, block_length, block_span, block_gap_run):
#    #print("\n")
#    #print('interval_starts', type(interval_starts), interval_starts.shape, interval_starts)
#    #print('interval_ends', type(interval_ends), interval_ends.shape, interval_ends)
#    sites = create_ranges(np.array((interval_starts, interval_ends), dtype=np.int64).T)
#    if sites is None: 
#        return None
#    block_sites = np.concatenate([
#        x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length) 
#            for x in np.split(sites, np.where(np.diff(sites) > block_gap_run)[0] + 1)
#        ]) 
#    block_span_valid_mask = (((block_sites[:, -1] - block_sites[:, 0]) + 1) <= block_span)
#    return block_sites[block_span_valid_mask]

# old
#def create_ranges(aranges):
#    # does the transformation form 0-based (BED) to 1-based (VCF) coordinate system 
#    # https://stackoverflow.com/a/47126435
#    l = (aranges[:, 1] - aranges[:, 0])
#    clens = l.cumsum()
#    if np.any(clens):
#        ids = np.ones(clens[-1], dtype=np.int64)
#        ids[0] = aranges[0, 0]
#        ids[clens[:-1]] = aranges[1:, 0] - aranges[:-1, 1] + 1
#        return ids.cumsum()
#    return None

##################################################################

    # def o_set_intervals(self, bed_f):
    #     # [needs fixing]
    #     meta = self._get_meta('seqs')
    #     query_sequences = set(meta['seq_names'])
    #     df = parse_csv(
    #         csv_f=bed_f, 
    #         sep="\t", 
    #         usecols=[0, 1, 2, 4], 
    #         dtype={'sequence': 'category', 'start': 'int64', 'end': 'int64', 'samples': 'category'},
    #         header=None)
    #     print(df)
    #     intervals_df = df[df['sequence'].isin(set(meta['seq_names']))].sort_values(['sequence', 'start'], ascending=[True, True]).reset_index(drop=True)
    #     intervals_df = pd.concat([intervals_df, intervals_df.samples.str.get_dummies(sep=',').filter(meta['samples'])], axis=1).drop(columns=['samples'])
    #     intervals_df_samples = [sample for sample in intervals_df.columns[3:]]
    #     print(intervals_df)
    #     query_samples = ordered_intersect(a=intervals_df_samples, b=meta['samples'], order='a')
    #     print(query_samples)
    #     intervals_df['length'] = (intervals_df['end'] - intervals_df['start'])
    #     # Check if all samples were found
    #     if set(query_samples) != set(meta['samples']):
    #             sys.exit("[X] The following samples in SAMPLE_FILE were not found in BED_FILE: %s" % (
    #                 ", ".join(list(set(meta['samples_sorted']).difference(set(query_samples))))))
    #     # Set up counts arrays
    #     count_shape = (len(meta['seq_names']), len(query_samples))
    #     count_bases_samples = np.zeros(count_shape, dtype=np.int64)
    #     print(intervals_df)
    #     for idx, (sequence, _df) in tqdm(enumerate(intervals_df.groupby(['sequence'], observed=True)), total=len(query_sequences), desc="[%] Reading intervals", ncols=100):
    #         interval_matrix = _df[query_samples].to_numpy()
    #         #length_matrix = np.repeat(_df['length'].to_numpy(), interval_matrix.shape[1]).reshape(interval_matrix.shape)        
    #         #length_matrix[interval_matrix == 0] = 0 # sets length to 0 if interval not present in interval_matrix 
    #         print('interval_matrix', interval_matrix)
    #         length_matrix = interval_matrix * _df['length'].to_numpy().reshape(-1, 1)
    #         print('length_matrix', length_matrix)
    #         count_bases_samples[idx,:] = np.sum(length_matrix, axis=0)
    #         self.data.create_dataset("seqs/%s/intervals/matrix" % sequence, data=interval_matrix)
    #         self.data.create_dataset("seqs/%s/intervals/starts" % sequence, data=_df['start'].to_numpy())
    #         self.data.create_dataset("seqs/%s/intervals/ends" % sequence, data=_df['end'].to_numpy())
    #     meta['intervals_span_sample'] = [int(x) for x in np.sum(count_bases_samples, axis=0)] # JSON encoder does not like numpy dtypes   
    #     meta['intervals_count'] = len(intervals_df.index)
    #     meta['intervals_span'] = int(intervals_df['length'].sum())
    #     meta['intervals_idx_by_sample'] = {sample: idx for idx, sample in enumerate(query_samples)}
    #     meta['bed_f'] = bed_f
    #     # QC plots
    #     #intervals_df['distance'] = np.where((intervals_df['sequence'] == intervals_df['sequence'].shift(-1)), (intervals_df['start'].shift(-1) - intervals_df['end']) + 1, np.nan)
    #     #distance_counter = intervals_df['distance'].dropna(how="any", inplace=False).value_counts()
    #     #length_counter = intervals_df['length'].value_counts()
    #     #distance_f = "%s.intervals.distance.png" % parameterObj.outprefix
    #     #plot_loglog(distance_counter, 'Distance to downstream BED interval', distance_f)
    #     #length_f = "%s.intervals.length.png" % parameterObj.outprefix
    #     #plot_loglog(length_counter, 'Length of BED interval', length_f)
    #     #count_sequences = intervals_df['sequence'].nunique()
    #     #count_intervals = len(intervals_df.index)
    #     #count_samples = len(query_samples)


    #     def variation_to_2d(variation, kmax_by_mutype=None):
    #     '''max_k is not capped by default (since it does not matter for variation arrays)'''
    #     max_k = np.array(list(kmax_by_mutype.values())) + 1 if kmax_by_mutype else None
    #     if variation.ndim == 2: 
    #         mutuples = np.clip(variation, 0, max_k)
    #     elif variation.ndim == 3:
    #         index = np.repeat(np.arange(variation.shape[0]), variation.shape[1]).reshape(-1, 1)
    #         mutuples = np.concatenate((index, np.clip(variation.reshape(-1, variation.shape[-1]), 0, max_k)), axis=-1)
    #     else:
    #         raise ValueError('variation.ndim is %r, should either be 2 (blocks) or 3 (windows)' % variation.ndim)
    #     mutuples_unique, counts = np.unique(mutuples, return_counts=True, axis=0)

    # def variation_to_bsfs(variation, kmax_by_mutype=None):
    #     '''max_k capped at (8,8,8,8) if not specified
    #     # we should use scipy/sparse in the long run ... 
    #     https://github.com/benbovy/xarray-simlab/issues/165
    #     https://github.com/zarr-developers/zarr-python/issues/152'''
    #     max_k = np.array(list(kmax_by_mutype.values())) + 1 if kmax_by_mutype else np.array([8,8,8,8]) 
    #     if variation.ndim == 2: 
    #         mutuples = np.clip(variation, 0, max_k)
    #     elif variation.ndim == 3:
    #         index = np.repeat(np.arange(variation.shape[0]), variation.shape[1]).reshape(-1, 1)
    #         mutuples = np.concatenate((index, np.clip(variation.reshape(-1, variation.shape[-1]), 0, max_k)), axis=-1)
    #     else:
    #         raise ValueError('variation.ndim is %r, should either be 2 (blocks) or 3 (windows)' % variation.ndim)
    #     try:
    #         mutuples_unique, counts = np.unique(mutuples, return_counts=True, axis=0)
    #         out = np.zeros(tuple(np.max(mutuples_unique, axis=0) + 1), np.uint64) 
    #         out[tuple(mutuples_unique.T)] = counts
    #     except MemoryError as e:
    #         sys.exit('[X] variation_to_bsfs() ran out of memory. %s. Try specifying lower k-max values.' % str(e))
    #     return out

# def szudzik_pairing(folded_minor_allele_counts):
#     # adapted from: https://drhagen.com/blog/superior-pairing-function/
#     return np.where(
#             (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
#             np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
#             folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
#             )
#     #if isinstance(folded_minor_allele_counts, np.ndarray):
#     #    # assumes folded_minor_allele_counts is array with shape (n,2)
#     #    return np.where(
#     #        (folded_minor_allele_counts[:,0] >= folded_minor_allele_counts[:,1]),
#     #        np.square(folded_minor_allele_counts[:,0]) + folded_minor_allele_counts[:,0] + folded_minor_allele_counts[:,1],
#     #        folded_minor_allele_counts[:,0] + np.square(folded_minor_allele_counts[:,1])
#     #        )
#     #elif isinstance(folded_minor_allele_counts, tuple):
#     #    # assumes folded_minor_allele_counts is tuple of size 2
#     #    a, b = folded_minor_allele_counts
#     #    if a >= b:
#     #        return (a**2) + a + b 
#     #    return a + (b**2)
#     #else:
#     #    pass

############## Only needed once we have demand for multidimensional pairing function
# def multidimensional_box_pairing(lengths: List[int], indexes: List[int]) -> int:
#     n = len(lengths)
#     index = 0
#     dimension_product = 1
#     for dimension in reversed(range(n)):
#         index += indexes[dimension] * dimension_product
#         dimension_product *= lengths[dimension]
#     return index
# def multidimensional_szudzik_pairing(*indexes: int) -> int:
#     n = len(indexes)
#     if n == 0:
#         return 0
#     shell = max(indexes)
#     def recursive_index(dim: int):
#         slice_dims = n - dim - 1
#         subshell_count = (shell + 1) ** slice_dims - shell ** slice_dims
#         index_i = indexes[dim]
#         if index_i == shell:
#             return subshell_count * shell + multidimensional_box_pairing([shell + 1] * slice_dims, indexes[dim + 1:])
#         else:
#             return subshell_count * index_i + recursive_index(dim + 1)
#     return shell ** n + recursive_index(0)