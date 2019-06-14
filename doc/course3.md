# [Course]

# Introduction

- We will analyse variation within and between two populations of Heliconius melpomene ('chi' and 'ros') along chromosome 18 ('Chr18') of the reference genome Hmel2 (Davey et al, 2017). Each population is composed of 10 samples. 
- Variation data is typically encoded in [VCF]{http://www.internationalgenome.org/wiki/Analysis/vcf4.0/} (Variant Call Format) files. 
- VCF files are the result of DNA sequencing data being subjected to an complex analytical (and often, non-linear!) process composed of read mapping (against a reference genome), variant calling and data QC.   
- Typically, only sites that show variation in at least one sample compared to the reference are recorded in a VCF file. This is based on the assumption that all sites in the reference were uniformly covered in all sequenced samples.
- In reality 
  this is that  is However, this means that all positions not recorded in a VCF files  monomorphic 
Variation data based on DNA sequencing data is generated via stochastic process. Sequencing success of a particular genomic reads are generated depending on template availability

- Dealing with multiallelic sites properly
# Running gIMble

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

# References
- Davey, J.W., Barker, S.L., Rastas, P.M., Pinharanda, A., Martin, S.H., Durbin, R., McMillan, W.O., Merrill, R.M. and Jiggins, C.D., 2017. No evidence for maintenance of a sympatric Heliconius species barrier by chromosomal inversions. Evolution Letters, 1(3), pp.138-154.