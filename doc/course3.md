# [Course]

# Introduction

- We will analyse variation within and between two populations of Heliconius melpomene ('chi' and 'ros') along chromosome 18 ('Chr18') of the reference genome Hmel2 (Davey et al, 2017). Each population is composed of 10 samples. 
- Variation data is typically encoded in [VCF]{http://www.internationalgenome.org/wiki/Analysis/vcf4.0/} (Variant Call Format) files. 
- VCF files are the result of DNA sequencing data being subjected to an complex analytical (and often, non-linear!) process composed of read mapping (against a reference genome), variant calling and data QC.   
- Typically, only sites that show variation in at least one sample compared to the reference are recorded in a VCF file. 
- This is based on the assumption that all sites in the reference were uniformly covered in all sequenced samples.
- In reality 
  this is that  is However, this means that all positions not recorded in a VCF files  monomorphic 
Variation data based on DNA sequencing data is generated via stochastic process. Sequencing success of a particular genomic reads are generated depending on template availability

- Dealing with multiallelic sites properly
# Running gIMble

talk about Missingness -> Missing_Genotypes

[Concepts]
- interspecies pair:

- block = genomic region covered by at least on sample in each population
	- fixed length (sum of bases) 
	- maximum span (end - start)
	- composed of intervals covered by interpopulational pairs (IPs)


- genotype: sample genotype in VCF GT field
	- 
- mutuple = mutation configuration
	- a block has a mutuple for each IP. 
	- with format: [hetA, fixed, hetB, hetAB]

- kmax

## Analyses
- A) Only use regions present in ALL samples from both populations
- B) Only use regions present in at least 5 samples in each population
- B) Only use regions present in at least 1 sample in each population

# Enter gIMble directory
cd gIMble

## Get data
cd gIMble && mkdir -p input && wget https://www.dropbox.com/s/hhz959alniylfpm/hmel.chr18.tar.gz && tar zxvf hmel.chr18.tar.gz -C input/ && rm hmel.chr18.tar.gz

# a = 5
# blocks : 1 min
# variants : <12 min
# modify blocks 1 min
# modify variants 1 min
# windows 4 min
# gridsearch 1 min
# likelihood 12 min : 674.1598691130057s using 30 iterations (Composite Likelihood = -502631354.145): theta=0.5623, Time=0.3828

./gIMble blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-b /data/hmel.chr18/hmel.chr18.bed \
	-o hmel.chr18.n_5 \
	-n 5 \
	-t 4 && \
./gIMble variants \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v /data/hmel.chr18/hmel.chr18.vcf.gz \
	-o hmel.chr18.n_5 \
	-b hmel.chr18.n_5.blocks.h5 \
	-t 4 && \
./gIMble modify blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-c /data/hmel.chr18/hmel.chr18.chrom_coordinates.txt \
	-o hmel.chr18.n_5 \
	-b hmel.chr18.n_5.blocks.h5 && \
./gIMble modify variants \
	-v hmel.chr18.n_5.variants.h5 \
	-o hmel.chr18.n_5 \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-b hmel.chr18.n_5.modified.blocks.h5 \
	-v hmel.chr18.n_5.modified.variants.h5 \
	-o hmel.chr18.n_5 \
	-w 500 \
	-l 100 \
	-t 4 && \
./gIMble gridsearch \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-o hmel.chr18.n_5 \
	-w hmel.chr18.n_5.windows.h5 \
	-l models/model.IM.M_D2A.MM_D2A.txt \
	--grid grid/grid.csv \
	-A "chi" \
	-k 2 \
	-t 4 && \
./gIMble likelihood \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v hmel.all_mutype_counts.h5 \
	-l models/model.divergence.txt \
	-A "chi" \
	-k 2 \
	--theta_low 0.5 \
	--theta_high 1.4 \
	--time_low 0.2 \
	--time_high 0.6 \
	--derived_Ne 1.0 \
	 -t 4 && \ 
#./gIMble likelihood \
#	-s input/hmel.samples.csv \
#	-g input/hmel.chr18.genomefile \
#	-l models/model.divergence.txt \
#	-A "chi" \
#	-v hmel.all_mutype_counts.h5 \
#	-k 2 \
#	--theta_low 0.5 \
#	--theta_high 1.4 \
#	--time_low 0.2 \
#	--time_high 0.6 \
#	--derived_low 0.5 \
#	--derived_high 2.0 \
#	-t 4


# n= 10
# blocks : 1 min
# variants : <6 min
# modify blocks 1 min
# modify variants 1 min
# windows 2 min

# likelihood 12 min : 674.1598691130057s using 30 iterations (Composite Likelihood = -502631354.145): theta=0.5623, Time=0.3828
# gridsearch 1 min

677.0936897310021s using 30 iterations (Composite Likelihood = -502631354.145): theta=0.5623, Time=0.3828

./gIMble blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-b /data/hmel.chr18/hmel.chr18.bed \
	-o hmel.chr18.n_10 \
	-n 10 \
	-t 4 && \
./gIMble variants \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v /data/hmel.chr18/hmel.chr18.vcf.gz \
	-o hmel.chr18.n_10 \
	-b hmel.chr18.n_10.blocks.h5 \
	-t 4 && \
./gIMble modify blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-c /data/hmel.chr18/hmel.chr18.chrom_coordinates.txt \
	-o hmel.chr18.n_10 \
	-b hmel.chr18.n_10.blocks.h5 && \
./gIMble modify variants \
	-v hmel.chr18.n_10.variants.h5 \
	-o hmel.chr18.n_10 \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-b hmel.chr18.n_10.modified.blocks.h5 \
	-v hmel.chr18.n_10.modified.variants.h5 \
	-o hmel.chr18.n_10 \
	-w 500 \
	-l 100 \
	-t 4 && \
./gIMble gridsearch \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-o hmel.chr18.n_10 \
	-w hmel.chr18.n_10.windows.h5 \
	-l models/model.IM.M_D2A.MM_D2A.txt \
	--grid grid/grid.csv \
	-A "chi" \
	-k 2 \
	-t 4 && \
./gIMble likelihood \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v hmel.all_mutype_counts.h5 \
	-l models/model.divergence.txt \
	-A "chi" \
	-k 2 \
	--theta_low 0.5 \
	--theta_high 1.4 \
	--time_low 0.2 \
	--time_high 0.6 \
	--derived_Ne 1.0 \
	 -t 4


Full commands:

# n=5
jovyan@jupyter-drl:~/gIMble$ ./gIMble blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-b /data/hmel.chr18/hmel.chr18.bed \
	-o hmel.chr18.min_10_sample \
	-a 10 \
	-t 4 && \
./gIMble variants \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v /data/hmel.chr18/hmel.chr18.vcf.gz \
	-b hmel.chr18.min_10_sample.blocks.h5 \
	-o hmel.chr18.min_10_sample \
	-t 4 && \
./gIMble modify blocks \
	-o hmel.chr18.min_10_sample \
	-s /data/hmel.chr18/hmel.samples.csv \
	-c /data/hmel.chr18/hmel.chr18.chrom_coordinates.txt \
	-b hmel.chr18.min_10_sample.blocks.h5 \
	-g /data/hmel.chr18/hmel.chr18.genomefile && \
./gIMble modify variants \
	-v hmel.chr18.min_10_sample.variants.h5 \
	-o hmel.chr18.min_10_sample \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-o hmel.chr18.min_10_sample \
	-s /data/hmel.chr18/hmel.samples.csv \
	-b hmel.chr18.min_10_sample.modified.blocks.h5 \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-v hmel.chr18.min_10_sample.modified.variants.h5 \
	-w 500 \
	-l 100 \
	-t 4 && \
./gIMble gridsearch \
	-s input/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-l models/model.IM.M_D2A.MM_D2A.txt \
	-A "chi" \
	-w hmel.chr18.min_10_sample.windows.h5 \
	-k 2 \
	-o hmel.chr18.min_10_sample \
	-t 4 \
	--mu 1.9e-9 \
	--block_size 64 \
	--time_MLE 4e6 \
	--grid grid/grid.csv && \
./gIMble likelihood \
	-s input/hmel.samples.csv \
	-g input/hmel.chr18.genomefile \
	-l models/model.divergence.txt \
	-A "chi" \
	-v hmel.all_mutype_counts.h5 \
	-k 2 \
	--theta_low 0.5 \
	--theta_high 1.4 \
	--time_low 0.2 \
	--time_high 0.6 \
	--derived_Ne 1.0
	-t 4
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.009s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.011s
[#] Processing BED file ...
[#] Splitting bed file into chunks for downstream processing (this might take a while) ...
[+] Found 1,023,304 BED intervals adding up to 10,627,982 b (63.25% of genome) ...
[+] BED intervals processed (232.92MB)
[#] Generating blocks ...
[+] Analysing 19698 BED regions using 4 threads ...
[%] : 100%|█████████████████████████████████████████████████████| 19.7k/19.7k [00:25<00:00, 763it/s]
[+] Made 71,814 blocks covering 4,596,096 b (43.25% of BED intervals, 27.35% of genome) (937.11MB)
[#] Generating output ...
[%] : 100%|███████████████████████████████████████████████████| 71.8k/71.8k [00:06<00:00, 10.4kit/s]
[>] Created hdf5 store: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.blocks.h5'
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.distance.png'
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.blocks_per_sample.png'
[+] Total runtime: 72.048s
[#] Parsing parameters ...
[+] Read parameters in 0.001s (95.10MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.010s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.014s
[#] Loading blocks ...
[%] : 100%|█████████████████████████████████████████████████| 71814/71814 [00:25<00:00, 2803.26it/s]
[+] Read 71,814 blocks covering 4,596,096 b (27.35% of genome) (3008.52MB)
[#] Fetching variants for blocks...
[%] : 100%|█████████████████████████████████████████████████████| 71.8k/71.8k [04:01<00:00, 297it/s]
[+] VCF parsed in 243.495s (4364.04MB)
[#] Generating output ...
[#] Analysing variants in blocks...
[%] : 100%|███████████████████████████████████████████████████| 71.8k/71.8k [00:52<00:00, 1.36kit/s]
[+] Monomorphic blocks: 1531242 (21.32% of blocks)
[#] Creating dataframe of global mutuple tallies...
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.mutuple_barchart.png'
[#] Creating dataframe of blocks...
[#] Creating dataframe for samples...
### Sample metrics: hmel.chr18.min_10_sample
    population_id    sample_id       blocks    sites    heterozygosity    missingness
--  ---------------  ------------  --------  -------  ----------------  -------------
 0  chi              chi.CJ564        69968  4477952         0.0174278     0.00447124
 1  chi              chi.CAM25137     69968  4477952         0.0172467     0.00518116
 2  chi              chi.CAM585       69968  4477952         0.0172802     0.00531113
 3  chi              chi.CJ553        69968  4477952         0.016468      0.00451903
 4  chi              chi.CAM25091     69968  4477952         0.0169508     0.00543195
 5  chi              chi.CAM580       69968  4477952         0.0168794     0.00580377
 6  chi              chi.CJ565        69968  4477952         0.0178278     0.00440715
 7  chi              chi.CJ560        69968  4477952         0.0174238     0.00448352
 8  chi              chi.CAM586       69968  4477952         0.0169897     0.0056575
 9  chi              chi.CAM582       69968  4477952         0.0171469     0.00538974
10  ros              ros.CAM1880      69968  4477952         0.016966      0.00323719
11  ros              ros.CJ531        69968  4477952         0.0169736     0.00343371
12  ros              ros.CAM2552      69968  4477952         0.0162594     0.00383903
13  ros              ros.CJ546        69968  4477952         0.0165725     0.0035811
14  ros              ros.CAM2519      69968  4477952         0.0170149     0.00349267
15  ros              ros.CAM2059      69968  4477952         0.0162577     0.00357619
16  ros              ros.CJ2071       69968  4477952         0.0170891     0.00305296
17  ros              ros.CAM1841      69968  4477952         0.016723      0.00338213
18  ros              ros.CJ533        69968  4477952         0.0168081     0.00355698
19  ros              ros.CAM2045      69968  4477952         0.0167742     0.00383256
[#] Creating dataframe for populations...
### Population metrics: hmel.chr18.min_10_sample
    population_id      blocks    sites    heterozygosity    missingness
--  ---------------  --------  -------  ----------------  -------------
 0  chi                 69968  4477952         0.0171641     0.00506562
 1  ros                 69968  4477952         0.0167439     0.00349845
[#] Creating dataframe for dataset...
### Dataset metrics: hmel.chr18.min_10_sample
      blocks       sites         FGVs     pi_chi     pi_ros        dxy       fst
--  --------  ----------  -----------  ---------  ---------  ---------  --------
 0     71814  4.5961e+06  0.000560041  0.0155308  0.0148279  0.0218143  0.179353
[>] Created hdf5 store: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.variants.h5'
[+] Total runtime: 348.116s
[+] Read parameters in 0.000s (91.49MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.011s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.014s
[#] Parse coordinate transformation file ...
[%] : 100%|██████████████████████████████████████████████████████| 19/19 [00:00<00:00, 55379.97it/s]
[+] Read file in 0.013s (94.47MB)
[#] Loading blocks ...
[%] : 100%|█████████████████████████████████████████████████| 71814/71814 [00:24<00:00, 2891.81it/s]
[+] Read 71,814 blocks covering 4,596,096 b (27.35% of genome) (3004.80MB)
[#] Transforming coordinates ...
[%] : 100%|████████████████████████████████████████████████████| 71.8k/71.8k [00:00<00:00, 295kit/s]
[+] Transformed coordinates in 0.261s (3004.80MB)
[#] Generating output...
[%] : 100%|███████████████████████████████████████████████████| 71.8k/71.8k [00:06<00:00, 11.8kit/s]
[#] 71814 of 71814 blocks (100.00%) are being retained ...
[>] Created hdf5 store: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.modified.blocks.h5'
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.modified.distance.png'
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.modified.blocks_per_sample.png'
[+] Total runtime: 35.291s
[+] Read parameters in 0.000s (91.26MB)
[#] Loading variants ...
[+] Excluded 67260 out of 7181400 mutuples (0.94%)
[#] Analysing variants in blocks...
[%] : 100%|███████████████████████████████████████████████████| 7.11M/7.11M [00:03<00:00, 2.35Mit/s]
[#] Creating dataframe of global mutuple tallies...
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.modified.mutuple_barchart.png'
[>] Created hdf5 store: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.modified.variants.h5'
[+] Total runtime: 36.926s
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.010s.
[+] Read 1 sequences with total length of 16,803,890 b in 0.013s
[#] Loading blocks ...
[%] : 100%|█████████████████████████████████████████████████| 71814/71814 [00:25<00:00, 2860.07it/s]
[+] Read 71,814 blocks covering 4,596,096 b (27.35% of genome) (3008.16MB)
[#] Loading variants ...
[%] : 100%|███████████████████████████████████████████| 7114140/7114140 [00:30<00:00, 235582.93it/s]
[#] Making windows ...
[%] : 100%|███████████████████████████████████████████████████████| 40.0/40.0 [00:00<00:00, 522it/s]
[#] Generating output ...
[%] : 100%|████████████████████████████████████████████████████████| 341/341 [00:40<00:00, 8.32it/s]
[#] Creating dataframe of window mutuple tallies...
[#] Creating dataframe of window metrics...
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.genome_scan.png'
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.pi_scatter.png'
[>] Created hdf5 store: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.windows.h5'
[+] Total runtime: 116.648s
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.010s.
[+] Read 1 sequences with total length of 16803890 b in 0.013s
[+] Parsing grid parameter file: /home/jovyan/gIMble/grid/grid.csv ...
[%] : 100%|█████████████████████████████████████████████| 306432/306432 [00:00<00:00, 668487.22it/s]
[+] Generating all mutuples for kmax = 2 ...
[=] Generated 45 mutuples (with kmax = 2)
[+] Ancestor is chi ...
[%] : 100%|█████████████████████████████████████████████| 151072/151072 [00:00<00:00, 238901.18it/s]
[+] Calculating likelihoods of model parameters ...
[%] : 100%|███████████████████████████████████████████████████████| 341/341 [00:08<00:00, 37.99it/s]
[+] Generating output ...
[%] : 100%|██████████████████████████████████████████████████████| 341/341 [00:03<00:00, 111.70it/s]
[>] Created: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.parameter_scan.png'
[>] Created hdf5 store: '/home/jovyan/gIMble/hmel.chr18.min_10_sample.composite_likelihoods.h5'
[+] Total runtime: 19.83807822699964 seconds

# References
- Davey, J.W., Barker, S.L., Rastas, P.M., Pinharanda, A., Martin, S.H., Durbin, R., McMillan, W.O., Merrill, R.M. and Jiggins, C.D., 2017. No evidence for maintenance of a sympatric Heliconius species barrier by chromosomal inversions. Evolution Letters, 1(3), pp.138-154.



########################local

# n = 10
./gIMble blocks \
	-o hmel.chr18.n_10 \
	-s input/hmel.samples.csv \
	-b input/hmel.chr18.bed \
	-g input/hmel.chr18.genomefile \
	-t 4 \
	-n 10 && \
./gIMble variants \
	-o hmel.chr18.n_10 \
	-s input/hmel.samples.csv \
	-b hmel.chr18.n_10.blocks.h5 \
	-g input/hmel.chr18.genomefile \
	-t 4 \
	-v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks \
	-o hmel.chr18.n_10 \
	-s input/hmel.samples.csv \
	-c input/hmel.chr18.chrom_coordinates.txt \
	-b hmel.chr18.n_10.blocks.h5 \
	-g input/hmel.chr18.genomefile && \
./gIMble modify variants \
	-v hmel.chr18.n_10.variants.h5 \
	-o hmel.chr18.n_10 \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-o hmel.chr18.n_10 \
	-s input/hmel.samples.csv \
	-b hmel.chr18.n_10.blocks.modified.h5 \
	-g input/hmel.chr18.new_coordinates.genomefile \
	-v hmel.chr18.n_10.variants.modified.h5 \
	-w 500 \
	-l 100

# n = 5

./gIMble blocks \
	-o hmel.chr18.n_5 \
	-s input/hmel.samples.csv \
	-b input/hmel.chr18.bed \
	-g input/hmel.chr18.genomefile \
	-t 4 \
	-n 5 && \
./gIMble variants \
	-o hmel.chr18.n_5 \
	-s input/hmel.samples.csv \
	-b hmel.chr18.n_5.blocks.h5 \
	-g input/hmel.chr18.genomefile \
	-t 4 \
	-v input/hmel.chr18.vcf.gz && \
./gIMble modify blocks \
	-o hmel.chr18.n_5 \
	-s input/hmel.samples.csv \
	-c input/hmel.chr18.chrom_coordinates.txt \
	-b hmel.chr18.n_5.blocks.h5 \
	-g input/hmel.chr18.genomefile && \
./gIMble modify variants \
	-v hmel.chr18.n_5.variants.h5 \
	-o hmel.chr18.n_5 \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-o hmel.chr18.n_5 \
	-s input/hmel.samples.csv \
	-b hmel.chr18.n_5.modified.blocks.h5 \
	-g input/hmel.chr18.new_coordinates.genomefile \
	-v hmel.chr18.n_5.modified.variants.h5 \
	-w 500 \
	-l 100
