# To do

#conda create -n gimble -y && source activate gimble && conda install -c conda-forge more-itertools tqdm scipy numpy matplotlib sympy giac networkx psutil pandas docopt pytables tabulate htop git -y && conda install -c bioconda pysam -y

./gIMble blocks -s test/test.samples.csv -g test/test.genomefile -b test/test.bed -o test_output/ && \
./gIMble variants -s test/test.samples.csv -g test/test.genomefile -v test/test.vcf.gz -b test_output/gimble.blocks.h5 -o test_output/ && \
./gIMble modify blocks -s test/test.samples.csv -g test/test.genomefile -b test_output/gimble.blocks.h5 -o test_output/ -c test/test.coord.txt && \
./gIMble modify variants -o test_output/ -v test_output/gimble.variants.h5 -m 4 && \
./gIMble windows -s test/test.samples.csv -g test/test.coord.genomefile -b test_output/gimble.blocks.modified.h5 -v test_output/gimble.variants.modified.h5 -w 3 -l 1 -o test_output/ 

./gIMble blocks -o hmel.chr18.Hmel218019 -s input/hmel.samples.csv -b input/hmel.chr18.Hmel218019.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.Hmel218019.genomefile -t 4 -a 5
./gIMble variants -o hmel.chr18.Hmel218019 -s input/hmel.samples.csv -b hmel.chr18.Hmel218019.blocks.h5 -g input/hmel.chr18.Hmel218019.genomefile -v input/hmel.chr18.vcf.gz -t 4
./gIMble modify blocks -o hmel.chr18.Hmel218019 -s input/hmel.samples.csv -c input/hmel.chr18.Hmel218019.chrom_coordinates.txt -b hmel.chr18.Hmel218019.blocks.h5  -g input/hmel.chr18.Hmel218019.genomefile
./gIMble modify variants -v hmel.chr18.Hmel218019.variants.h5 -o hmel.chr18.Hmel218019 -m 4 -M 100
./gIMble windows -o hmel.chr18.Hmel218019 -s input/hmel.samples.csv -b hmel.chr18.Hmel218019.blocks.modified.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.Hmel218019.variants.modified.h5 -w 50 -l 10 -m 1000000 -t 4

conda create -n gimble -y && source activate gimble && conda install -c conda-forge more-itertools tqdm scipy numpy matplotlib sympy giac networkx psutil pandas docopt pytables tabulate htop git -y && conda install -c bioconda pysam -y && cd gIMble/ && python lib/probs_sympy.py -p output/master.2_pop.2_ploidy.E.paths.txt -t 1 -A 1.0 -D 1.0 -M 2.34 -m 1.2 -T 1.4 -P 10

./gIMble windows -o hmel.min_5_samples -s ../hmel_data/hmel.samples.csv -g ../hmel_data/hmel.autosomes.new_coordinates.genomefile -b hmel.min_5_samples.blocks.modified.h5 -v ../hmel.min_5_samples.variants.modified.h5 -w 500 -l 100 -t 24 -m 10000

# [VARIANTS]
- [+] Read parameters in 0.000s (66.26MB) REMOVE
- [+] Read 20 samples from 2 populations and generated 100 pairs in 0.009s. INDENT
- [+] Read 1 sequences with total length of 16803890 b in 0.011s INDENT
- [#] Loading blocks ... from h5
- Fetching variants ... from vcf
- [+] VCF parsed in 315.415s (4481.38MB) REMOVE
- [#] Generating output ... REMOVE
- [+] Monomorphic blocks: 2097277 (19.68% of blocks) format thousands
- [>] Creating dataframe of global mutuple tallies...
- [>] Creating dataframe of blocks...
- [>] Creating dataframe for samples...
- Missingness -> Missing_Genotypes
- Plot 20 most common mutypes and their counts as barcharts

# [MODIFY blocks]
- [+] Read parameters in 0.000s (66.26MB) REMOVE
- [#] Parse coordinate transformation file ... from FILE
- [+] Read file in 0.019s (67.98MB) REMOVE
- [+] Read 20 samples from 2 populations and generated 100 pairs in 0.009s. INDENT
- [+] Read 1 sequences with total length of 16803890 b in 0.011s INDENT
- [#] Loading blocks ... from FILE
- [+] Transformed coordinates in 0.469s (4467.70MB) add how many sequence not transformed
- [>] Created: '/Users/dlaetsch/git/gIMble/hmel.chr18.min_5_sample.distance.png' change to "modified"
- [>] Created: '/Users/dlaetsch/git/gIMble/hmel.chr18.min_5_sample.blocks_per_sample.png' change to "modified"

# [MODIFY variants]
- [+] Read parameters in 0.000s (66.26MB) REMOVE
- [+] Read 20 samples from 2 populations and generated 100 pairs in 0.009s. INDENT
- [+] Read 1 sequences with total length of 16803890 b in 0.011s INDENT
- default are None for filtering, exit if None
- Loading variants ... from file, load into blocks 
- Plot: for the 20 most common mutuples plot barchart of counts

# [WINDOWs]
- [#] Parsing parameters ... REMOVE
- [+] Read parameters in 0.001s (66.43MB).. REMOVE
- [+] Read 20 samples from 2 populations and generated 100 pairs in 0.009s. INDENT
- [+] Read 1 sequences with total length of 16803890 b in 0.011s INDENT, format thousand
- [#] Loading blocks ... from file
- [#] Loading variants ... from file, how many monomorphic?
- [#] Making windows .. how many are made?
- flip genomescan colourmap, make dots bigger, lines between dot
- pi_scatter, make dots small, increase alpha, remove lines from colorbar

#### prep for course

# full Hmel
- run complete to generate window dxy vs pi

# PLOTs

# Histogram of variants (sum of mutype counts) per block
# distribution of window spans
# fix offset of genome scan

#### Commands course

# ensure tabix is younger than VCF

# min_N = 5 Samples (go ahead for 4th part)
./gIMble blocks -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile -t 4 -a 5 && \
./gIMble variants -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -b hmel.chr18.min_5_sample.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -c input/hmel.chr18.chrom_coordinates.txt -b hmel.chr18.min_5_sample.blocks.h5 -g input/hmel.chr18.genomefile && \
./gIMble modify variants -v hmel.chr18.min_5_sample.variants.h5 -o hmel.chr18.min_5_sample -m 4 -M 100 && \
./gIMble windows -o hmel.chr18.min_5_sample -s input/hmel.samples.csv -b hmel.chr18.min_5_sample.modified.blocks.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.min_5_sample.modified.variants.h5 -w 500 -l 100 
./gIMble gridsearch -s ../gIMble/input/hmel.samples.csv -g input/hmel.chr18.genomefile -l models/model.IM.M_D2A.MM_D2A.txt -A "chi" -v hmel.chr18.min_5_sample.variants.modified.h5 -k 2 -o hmel.chr18.min_5_sample. -t 4 --mu 1.9e-9 --block_size 64 --migration_MLE 3.8866e-7 --time_MLE 4e6 --derived_MLE 0.5 --theta_low 0.4 --theta_high 1.2
# print probabilities for each mutation configuration

# which seeds work? what happens when A and D are switched

# ./gIMble likelihood -s input/hmel.chr18.samples.csv -g input/hmel.chr18.genomefile -l models/model.divergence.txt -A "chi" -v hmel.chr18/gimble.variants.modified.h5 -k 2 --derived_Ne 1.0 --migration 2.34 --theta 1.2 --time 1.4 -o hmel.chr18/ -t 1
# ./gIMble likelihood -s input/hmel.chr18.samples.csv -g input/hmel.chr18.genomefile -l models/model.divergence.txt -A "chi" -v hmel.chr18/gimble.variants.modified.h5 -k 2 --derived_Ne 1.0 --theta_low 0.2 --theta_high 2.0 --time_low 0.1 --time_high 2.0 -o hmel.chr18/ -t 1


# min_N = 1 Samples
./gIMble blocks -o hmel.chr18.min_1_sample -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile -t 4 -a 1 && \
./gIMble variants -o hmel.chr18.min_1_sample -s input/hmel.samples.csv -b hmel.chr18.min_1_sample.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks -o hmel.chr18.min_1_sample -s input/hmel.samples.csv -c input/hmel.chr18.chrom_coordinates.txt -b hmel.chr18.min_1_sample.blocks.h5 -g input/hmel.chr18.genomefile && \
./gIMble modify variants -v hmel.chr18.min_1_sample.variants.h5 -o hmel.chr18.min_1_sample -m 4 -M 100 && \
./gIMble windows -o hmel.chr18.min_1_sample -s input/hmel.samples.csv -b hmel.chr18.min_1_sample.modified.blocks.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.min_1_sample.modified.variants.h5 -w 500 -l 100

# min_N = 10 Samples
./gIMble blocks -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile -t 4 -a 10 && \
./gIMble variants -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -b hmel.chr18.min_10_sample.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -c input/hmel.chr18.chrom_coordinates.txt -b hmel.chr18.min_10_sample.blocks.h5 -g input/hmel.chr18.genomefile && \
./gIMble modify variants -v hmel.chr18.min_10_sample.variants.h5 -o hmel.chr18.min_10_sample -m 4 -M 100 && \
./gIMble windows -o hmel.chr18.min_10_sample -s input/hmel.samples.csv -b hmel.chr18.min_10_sample.blocks.modified.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.min_10_sample.variants.modified.h5 -w 500 -l 100






# Hmel all (SORTED BED)
# 10 samples
/gIMble blocks -o hmel.min_10_samples -s ../hmel_data/hmel.samples.csv -b ../hmel_data/hmel.multiinter.samples_as_string.only_intergenic.sorted.bed -g ../hmel_data/hmel.autosomes.genomefile -t 60 -a 10
[+] Made 916,620 blocks covering 58,663,680 b (39.46% of BED intervals, 22.44% of genome) (12304.05MB)
./gIMble variants -o hmel.min_10_samples -s ../hmel_data/hmel.samples.csv -b hmel.min_10_samples.blocks.h5 -g ../hmel_data/hmel.autosomes.genomefile -t 24 -v ../hmel_data/ros10_chi10.DP8MIN2MAC1.vcf.gz

# 5 samples
./gIMble blocks -o hmel.min_5_samples -s ../hmel_data/hmel.samples.csv -b ../hmel_data/hmel.multiinter.samples_as_string.only_intergenic.sorted.bed -g ../hmel_data/hmel.autosomes.genomefile -t 32 -a 5
[+] Made 1,827,522 blocks covering 116,961,408 b (78.67% of BED intervals, 44.75% of genome)
./gIMble variants -o hmel.min_5_samples -s ../hmel_data/hmel.samples.csv -b hmel.min_5_samples.blocks.h5 -g ../hmel_data/hmel.autosomes.genomefile -t 24 -v ../hmel_data/ros10_chi10.DP8MIN2MAC1.vcf.gz

# 1 sample
./gIMble blocks -o hmel.min_1_sample -s ../hmel_data/hmel.samples.csv -b ../hmel_data/hmel.multiinter.samples_as_string.only_intergenic.bed -g ../hmel_data/hmel.autosomes.genomefile -t 60 -a 10
./gIMble variants -o hmel.chr18.min_5_sample/ -s input/hmel.samples.csv -b hmel.chr18.min_5_sample/gimble.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz

# [MAIN]
- generate sample.h5
- add number of samples to header
- implement min interval sample/pair filter 
	- all samples
	- one sample in each pop
	- half
- overlap in windows always 1/5
- fix plots
	- make dots bigger, invert colormap for genomescan

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

MODIFY:
- uncouple from blockObj parsing, just dataframes (or avoid plotting distance.png again) 
VARIANTS:
- print possible values for k_max and how much would be filtered

- Plots 

Course:
- What happens if one does not deal with missing data
- Plot sample coverage of windows
[To do]                      
plot:
    - improve setting of figure dimensions based on graph complexity
    - legend with edge types
    - account for overlapping edges (using width/styles/colours)
    - adjust alpha/colour of edge labels

[To Do]
- output distribution of distance between blocks
- output barchart of bases blocked per sample
- output bases blocked vs per sequence
[ To Do ]
 - Variants: filter blocks on 4 missing both for global counts and window counts
    # PLOTs

# Histogram of variants (sum of mutype counts) per block
# distribution of window spans
# fix offset of genome scan

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
- improve parallelisation, pandas is inefficient

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