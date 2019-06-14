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

## gIMble blocks (<8 min)
./gIMble blocks \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-b /data/hmel.chr18/hmel.chr18.bed \
	-o hmel.chr18.min_5_sample \
	-a 5 \
	-t 4

## gIMble variants (<15 min)
./gIMble variants \
	-s /data/hmel.chr18/hmel.samples.csv \
	-g /data/hmel.chr18/hmel.chr18.genomefile \
	-v /data/hmel.chr18/hmel.chr18.vcf.gz \
	-b hmel.chr18.min_5_sample.blocks.h5 \
	-o hmel.chr18.min_5_sample \
	-t 4 

# 60.5s
./gIMble modify blocks \
	-o hmel.chr18.min_5_sample \
	-s /data/hmel.chr18/hmel.samples.csv \
	-c /data/hmel.chr18/hmel.chr18.chrom_coordinates.txt \
	-b hmel.chr18.min_5_sample.blocks.h5 \
	-g /data/hmel.chr18/hmel.chr18.genomefile && \
./gIMble modify variants \
	-v hmel.chr18.min_5_sample.variants.h5 \
	-o hmel.chr18.min_5_sample \
	-m 4 \
	-M 64 && \
./gIMble windows \
	-o hmel.chr18.min_5_sample \
	-s /data/hmel.chr18/hmel.samples.csv \
	-b hmel.chr18.min_5_sample.modified.blocks.h5 \
	-g /data/hmel.chr18/hmel.chr18.new_coordinates.genomefile \
	-v hmel.chr18.min_5_sample.modified.variants.h5 \
	-w 500 \
	-l 100 \
	-t 4


# References
- Davey, J.W., Barker, S.L., Rastas, P.M., Pinharanda, A., Martin, S.H., Durbin, R., McMillan, W.O., Merrill, R.M. and Jiggins, C.D., 2017. No evidence for maintenance of a sympatric Heliconius species barrier by chromosomal inversions. Evolution Letters, 1(3), pp.138-154.