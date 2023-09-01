gimble
=========

[![DOI](https://zenodo.org/badge/176883840.svg)](https://zenodo.org/badge/latestdoi/176883840)

# Table of contents

- [Installation](#installation)
- [Workflow](#workflow)
- [Usage](#usage)
- [Gimble modules](#gimble-modules)

# Installation

```
# 1a. clone gimble repository
>>> git clone https://github.com/LohseLab/gimble.git
# 1b. clone gimbleprep repository
>>> git clone https://github.com/LohseLab/gimbleprep.git
# 2. Install miniconda from https://conda.io/miniconda.html
# ...
# 3. Create the following conda environment 
>>> conda create -n gimble python=3.7.12 agemo bedtools bcftools samtools vcflib mosdepth=0.3.2 pysam numpy docopt tqdm pandas tabulate zarr scikit-allel parallel matplotlib msprime demes dask numcodecs python-newick nlopt -c conda-forge -c bioconda
# 4. Load the environment (needs to be activated when using gimble)
>>> conda activate gimble
# 5a. Start gimbleprep'ing ...
>>> (gimble) gimbleprep/gimbleprep --help
# 5b. Start gimble'ing ...
>>> (gimble) gIMble/gimble --help
```

# Workflow

![gIMble workflow](/docs/gIMble.workflow.jpg?raw=true "gIMble workflow")
**gIMble workflow**. `preprocess` (**0**) assures input data conforms to requirements; `parse` (**1**) reads data into a `gIMble`
store, the central data structure that holds all subsequent analysis. The modules `blocks` (**2**) and
`windows` (**3**) partition the data which is summarised as a tally (**4**) of blockwise mutation
configurations (bSFSs) either across all pair-blocks (blocks tally) or for pair-blocks in windows
(windows tally). Tallies may be used either in a bounded search of parameter space via the
module `optimize` (**5**) or to evaluate likelihoods over a parameter grid (which is precomputed using
`makegrid` (**6**)) via the `gridsearch` (**7**) module. The `simulate` (**8**) module allows coalescent
simulation of tallies (simulate tally) based on inferred parameter estimates (either global
estimates or gridsearch results of window-wise data). Simulated data can be analysed to quantify
the uncertainty (and potential bias) of parameter estimates. The results held within a `gIMble` store
can be described, written to column-based output files or removed using the modules `info` (**9**),
`query` (**10**), and `delete` (**11**).

# Usage
```
usage: gimble <module> [<args>...] [-V -h]

  [Input]
    preprocess            Install gimbleprep instead
    parse                 Parse files into GimbleStore
    blocks                Generate blocks from parsed data in GimbleStore (requires 'parse')
    windows               Generate windows from blocks in GimbleStore (requires 'blocks')
    tally                 Tally variation for inference (requires 'blocks' or 'windows')

  [Simulation]
    simulate              Simulate data based on specific parameters or gridsearch results  
    
  [Inference]
    optimize              Perform global parameter optimisation on tally/simulation
    makegrid              Precalculate grid of parameters
    gridsearch            Evaluate tally/simulation against a precomputed grid (requires 'makegrid')

  [Info]
    info                  Print metrics about data in GimbleStore
    list                  List information saved in GimbleStore
    query                 Extract information from GimbleStore
    delete                Delete information in GimbleStore

  [Experimental]
    partitioncds          Partition CDS sites in BED file by degeneracy in sample GTs

  [Options]
    -h, --help            Show this screen
    -V, --version         Show version
```

# Gimble modules
  * [preprocess](#preprocess)
    + [VCF processing details](#vcf-processing-details)
    + [BAM processing details](#bam-processing-details)
    + [Manually modify preprocessed files](#manually-modify-preprocessed-files)
  * [parse](#parse)
  * [blocks](#blocks)
  * [windows](#windows)
  * [info](#info)
  * [tally](#tally)
  * [optimize](#optimize)
  * [makegrid](#makegrid)
  * [gridsearch](#gridsearch)
  * [simulate](#simulate)


## preprocess
Note: The `preprocess` module has been replaced by `gimbleprep`. Everything else is identical.
The `preprocess` module assures that input files are adequately filtered and processed so that the `gimble` workflow can be completed successfully. 
While this processing of input files could be done more efficiently with other means, it has the advantage of generating a VCF file complies with `gimble` data requirements but which can also be used in alternative downstream analyses.

```
./gimbleprep -f FASTA -b BAM_DIR/ -v RAW.vcf.gz -k
```

Based on the supplied input files:
- `-f` : FASTA of reference
- `-b`: directory of BAM files, composed of readgroup-labelled reads mapped to reference 
- `-v`: compressed+indexed Freebayes VCF file 

the module produces the following output files:

- **genome file** (sequence_id, length) based on FASTA file
- **sample file** (sample_id) based on ReadGroupIDs in BAM files
- **coverage threshold report** for each BAM file
- **gimble VCF file** (see [VCF processing details](#vcf-processing-details))
- **gimble BED file** (see [BAM processing details](#bam-processing-details))
- log of executed commands

After running, output files require manual user input (see [Manually modify files](#manually-modify-preprocessed-files))

### VCF processing details
- MNPs are decomposed into SNPs
- variant sets are defined as
  ```
  - {RAW_VARIANTS} := all variants in VCF
  - {NONSNP}       := non-SNP variants 
  - {SNPGAP}       := all variants within +/- X b of {NONSNP} variants
  - {QUAL}         := all variants with QUAL below --min_qual
  - {BALANCE}      := all variants with any un-balanced allele observation (-e 'RPL<1 | RPR<1 | SAF<1 | SAR<1') 
  - {FAIL}         := {{NONSNP} U {SNPGAP} U {QUAL} U {BALANCE}} 
  - {VARIANTS}     := {RAW_VARIANTS} - {FAIL}```
The processed VCF file 
  - only contains variants from set `{VARIANTS}` 
  - only contains sample genotypes with sample read depths within coverage thresholds (others are set to missing, i.e. `./.`)

### BAM processing details 
- definition of BED region sets
  ```
  - {CALLABLE_SITES}   := for each BAM, regions along which sample read depths lie within coverage thresholds.
  - {CALLABLES}        := bedtools multiintersect of all {CALLABLE_SITES} regions across all BAMs/samples
  - {FAIL_SITES}       := sites excluded due to {FAIL} variant set (during VCF processing)
  - {SITES}            := {CALLABLES} - {FAIL_SITES}
  ```
Resulting BED file 
- only contains BED regions from set `{SITES}` 
- lists which samples are adequately covered along each region

#### Manually modify preprocessed files
+ `gimble.genomefile`:
    + **[Optional]** remove sequence IDs to ignore them in the analyses
+ `gimble.samples.csv` 
    + **[Required]** add population IDs the second column. Must be exactly 2 populations
    + **[Optional]** remove sample IDs to ignore them in the analyses 
+ `gimble.bed`
    + **[Recommended]** intersect with BED regions of interest to analyse particular genomic regions, e.g:  
    ```
    bedtools intersect -a gimble.bed -b my_intergenic_regions.bed > gimble.intergenic.bed
    ``` 
    
## parse
+ reads input data into GimbleStore 
```
./gimble parse -v gimble.vcf.gz -b gimble.intergenic.bed -g gimble.genomefile -s gimble.samples.csv -z analysis
```

## blocks
+ The output of this module determines which parts of the sampled sequences are available for analysis.
+ It uses the "callable" regions specified in the BED file and the variants contained in the VCF file to define sequence blocks of a fixed number of callable sites.
+ The blocking of genomic data is controlled by the parameters `--block_length` (number of callable sites in each block) and `--block_span` (maximum distance between the first and last site in a block)
+ Blocks are constructed independently for each sample pair, which ameliorates the asymmetry in coverage profiles among the samples due to stochastic variation in sequencing depth between samples and/or reference bias.
+ Optimal block length will be different for each dataset. The user is encouraged to explore parameter space. 
```
./gimble blocks -z analysis.z -l 64
```

## windows
+ Windows are constructed by traversing each sequence of the reference from start to end, incorporating the heterospecific pair-blocks (X) as they appear based on their start positions.
+ The parameter `--blocks` controls how many blocks are incorporated into each window and the parameter `--steps` by how many blocks the next window is shifted
+ `--blocks` should be chosen so that, given the number of interspecific pairs, enough blocks from each pair can be placed in a window. 
```
./gimble windows -z analysis.z -w 500 -s 100
```

## info
+ Lists basic metrics about the GimbleStore
+ Computes standard population genetic summary statistics ($\pi$, $d_{xy}$ and $H$ mean heterozygosity) based on blocks sampled between (X) and within species/populations (A and B). For for details see [gIMble_info.pdf](/docs/gIMble_info.pdf)
```
./gimble info -z analysis.z
```

## tally
+ Tallies variation for blocks or for blocks in windows into bSFSs
+ The bSFS is a tally of the mutation configurations of blocks which are themselves described by vectors of the form $`\underline{k}_i`$, which count the four possible mutation types found within a pair-block $i$.
+ parameter k-max is the max count per mutation type beyond which counts are treated as marginals. Order of mutation types is (hetB, hetA, hetAB, fixed)
```
./gimble tally -z analysis.z -k 2,2,2,2 -l blocks_kmax2 -t blocks
./gimble tally -z analysis.z -k 2,2,2,2 -l windows_kmax2 -t windows
```

## optimize
+ Searches parameter space for model parameters under a given model for a given data tally (based on parsed or simulated tallies) using an optimization algorithm with bounded constraints.
+ Given a set of bounds, optimization can be initiated either at the midpoint of the bounded parameter space or at a random start point. 
+ Optimizations finalize after user-defined stopping criteria are met
+ The user can assess convergence of optimizations by consulting the log-file. 
```
./gimble optimize -z analysis.z -l IM_BA_optimize -d tally/windows_kmax2 \
    -w -m IM_BA -r A -u 2.9e-09 -A 10_000,2_500_000 -B 10_000,1_000_000 \
    -C 1_000_000,5_000_000 -M 0,1e-5 -T 0,5_000_000 -g CRS2 -p 1 -e 19 -i 10_000
```

## makegrid
+ Computes a grid of parameter combinations for a list of parameters under a given model
+ Pre-computing the probabilities of bSFS configurations in a grid across parameter space is computationally efficient (relative to using an optimization algorithm) and therefore useful when we wish to interrogate data in replicate, i.e. across windows or simulation replicates.
+ Grid searches may be used either for an initial exploration of the likelihood surface (i.e. prior to defining the parameter space across which to run `optimize`), or to fit models to window-wise data.
```
./gimble makegrid -z analysis.z -m IM_AB \
    -b 64 -r A -u 2.9e-9 -k 2,2,2,2 \
    -A=200_000,3_000_000,12,lin \
    -B=100_000,2_000_000,12,lin \
    -C 100_000,2_000_000,12,lin \
    -T 4_256_034 -M 0,2.21E-06,16,lin \
    -p 48 -e 19 -l IM_BA_grid
```

## gridsearch
+ Searches a tally or simulation against a grid made with `makegrid`
```
gimble gridsearch -z analysis.z \
    -g makegrid/IM_BA_grid -d tally/windows_kmax2 -p 50
```

## simulate
+ Simulates tallies based on fixed parameters under a model
+ Simulates tallies based on results of gridsearches
```
# based on fixed parameters
./gimble simulate -z analysis.z \
    --seed 19 --replicates 100 --windows 11217 --blocks 500 \
    --block_length 64 -a 10 -b 10 \
    -k 2,2,2,2 -s IM_BA_sims_100reps -p 55 -m IM_BA \
    -A 1_407_027 -B 545_041 -C 923_309 -M 7.3778010583E-07 -u 2.9e-9 -T 4_256_034 

# based on gridsearch result
./gimble simulate -z analysis.z \
    --seed 19 --replicates 100 --windows 11217 --blocks 500 \
    --block_length 64 -a 10 -b 10 \
    --gridsearch_key gridsearch/windows_kmax2/IM_BA_grid \
    -k 2,2,2,2 -s IM_BA_grid -p 55 -u 2.9e-9 
```
