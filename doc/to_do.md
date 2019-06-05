# To do

# [General]
- Implement central hdf5 storage
	- Define entry points

- Changes to HDF5 storage:
	- deleting table with 'key' : store.remove('key')

- Compression
	- http://alimanfoo.github.io/2016/09/21/genotype-compression-benchmark.html
	- TL;DR 
	from zarr import Blosc
	compressors = (Blosc(cname='lz4', clevel=9, shuffle=0))

- Plots 


# [META]
- Command: ./gimble meta -s sample_file -g genome_file -p prefix
- Creates HDF5 storage which is used in all subsequent modules

## HDF5 interactions:
--------------------
- [C] sampleTable 			= ['sample_id', 'population_id'] 		
- [C] populationTable 		= ['population_id'] 		
- [C] datasetTable			= ['dataset_id', 'population_count', 'sample_count']
- [C] pairTable 			= ['pair_idx', 'sample_1', sample_2]		
- [C] sequenceTable 		= ['sequence_id', 'length']	 	
- [C] prefixVar 			= prefix		
- [C] pathVar 				= path			
- [C] countSequenceVar 		= sequence_count
- [C] countSampleVar 		= sample_count
- [C] countPairVar 			= pair_count	
--------------------

create_entityCollection(meta):
	- setup of entityCollection based on dfs

# [BLOCKS]

- fail if bed samples are not in SAMPLES
- fail if no blocks are made

## HDF5 interactions:
--------------------
- [R] sampleTable			| \ 				
- [R] pairTable				| => parse to entityCollection			
- [R] sequenceTable			| /  				
--------------------
- [C] blockableSitesVar		= sum of bases in bed intervals (i.e. blockable span) 	
- [C] blockTable 			= ['block_id', 'sequence_id', 'block_start', 'block_end', 'length', 'span', 'sample_ids', 'pair_idxs', 'count_samples', 'count_pairs']						
- [C] bedTable				= ['block_id', 'sequence_id', 'bed_start', 'bed_end']
--------------------
- [A] sampleTable			+= ['blocks', 'sites', 'percBlockableSites']
- [A] populationTable		+= ['blocksMean', 'sitesMean', 'percBlockableSites']
- [A] datasetTable			+= ['blocksMean', 'sitesMean', 'percBlockableSites']
--------------------

# [VARIANTS]

## HDF5 interactions:
--------------------
- [R] sampleTable			| \ 				
- [R] pairTable				| => parse to entityCollection			
- [R] sequenceTable			| /  				
- [R] blockTable			| \	
- [R] bedTable				| => parse blockObjs to entityCollection; parallelise creation of BlockObjs (should be way faster)  	
--------------------	
- [C] mutypeCounter			= ['count'] + MUTYPE_ORDER 
- [C] mutypeBlocks			= ['block_id', 'pair_idx'] + FULL_MUTYPE_ORDER 
--------------------
- [A] sampleTable			+= ['heterozygosity', 'missingness']
- [A] populationTable		+= ['heterozygosityMean', 'missingnessMean']
- [A] datasetTable			+= ['heterozygosityMean', 'missingnessMean']

# [MODIFY BLOCKS]

## HDF5 interactions:
--------------------
- [R] sampleTable			| \ 				
- [R] pairTable				| => parse to entityCollection			
- [R] sequenceTable			| /  				
- [R] blockTable			| \	
- [R] bedTable				| => parse blockObjs to entityCollection; parallelise creation of BlockObjs (should be way faster)  	
--------------------	
- [C] mutypeCounter			= ['count'] + MUTYPE_ORDER 
- [C] mutypeBlocks			= ['block_id', 'pair_idx'] + FULL_MUTYPE_ORDER 
--------------------
- [A] sampleTable			+= ['heterozygosity', 'missingness']
- [A] populationTable		+= ['heterozygosityMean', 'missingnessMean']
- [A] datasetTable			+= ['heterozygosityMean', 'missingnessMean']
--------------------
- [A] sampleTable			| => values for 'blocks', 'sites', 'percBlockableSites', heterozygosity, missingness have to be recalculated
- [A] populationTable		| 	
- [A] datasetTable			| 
--------------------
# [MODIFY VARIANTS]

# [WINDOWS]

## comments:
- Creating dataframe of window metrics is slow 
	- part of it will be access 

- barplots 
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set_style('white')
sns.set_style('ticks')
sns.set_context('paper')

mn = data.min()
mx = data.max()
x = np.arange(mn, mx + 1)
# split bincount in two to avoid memory error
y = np.bincount(data[:1000000].reshape(-1) - mn)
y += np.bincount(data[1000000:2000000].reshape(-1) - mn)

fig, ax = plt.subplots()
sns.despine(ax=ax, offset=10)
ax.bar(x + .1, y * 100 / y.sum(), width=.8)
ax.set_xticks(x + .5)
ax.set_xticklabels(x)
ax.set_xlabel('Value')
ax.set_ylabel('Frequency (%)')
ax.set_title('Distribution of data values');