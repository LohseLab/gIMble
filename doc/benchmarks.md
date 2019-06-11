# benchmarks

## ENVIRONMENT

install conda https://docs.conda.io/en/latest/miniconda.html
conda create -n gimble
source activate gimble
conda install -c conda-forge more-itertools tqdm scipy numpy matplotlib sympy giac networkx psutil pandas docopt pytables tabulate git
conda install -c bioconda pysam

git clone https://github.com/DRL/gIMble.git
git checkout dev

./gIMble blocks -p hmel.chr18.Hmel218018/ -s input/hmel.samples.csv -b input/hmel.chr18.Hmel218018.multiinter.samples_as_string.only_intergenic.bed -g input/hmel.chr18.Hmel218018.genomefile -t 4
./gIMble variants -p hmel.chr18.Hmel218018/ -s input/hmel.samples.csv -b hmel.chr18.Hmel218018/gimble.blocks.h5 -g input/hmel.chr18.Hmel218018.genomefile -v input/hmel.chr18.vcf.gz -t 4
./gIMble modify blocks -p hmel.chr18.Hmel218018/ -s input/hmel.samples.csv -c input/hmel.chr18.Hmel218018.chrom_coordinates.txt -b hmel.chr18.Hmel218018/gimble.blocks.h5  -g input/hmel.chr18.Hmel218018.genomefile
./gIMble modify variants -v hmel.chr18.Hmel218018/gimble.variants.h5 -p hmel.chr18.Hmel218018/ -m 4 -M 100
./gIMble windows -p hmel.chr18.Hmel218018/ -s input/hmel.samples.csv -b hmel.chr18.Hmel218018/gimble.blocks.modified.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18.Hmel218018/gimble.variants.modified.h5 -w 50 -l 10 -m 1000000 -t 4

## 2019-06-02

### minimal test
./gIMble blocks -s input/test.minimal.samples.shuffled.csv -b input/test.minimal.bed -g input/test.minimal.genomefile -p minimal
./gIMble variants -s input/test.minimal.samples.shuffled.csv -b /Users/dlaetsch/git/gIMble/minimal.blocks.h5 -g input/test.minimal.genomefile -p minimal -v input/test.vcf.gz
./gIMble modify -s input/test.minimal.samples.shuffled.csv -b /Users/dlaetsch/git/gIMble/minimal.blocks.h5 -g input/test.minimal.genomefile -p minimal -c input/test.minimal.coordinates.txt

### MacBook+Full

#### 1. [BLOCKS]

##### macbook
./gIMble blocks -p hmel/ -s input/hmel.samples.csv -b input/hmel.multiinter.samples_as_string.only_intergenic.bed -g input/hmel.autosomes.genomefile -t 4
[#] Parsing parameters ...
[+] Read parameters in 0.001s (53.06MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.008s.
[+] Read 747 sequences with total length of 261,379,282 b in 0.012s
[#] Processing BED file ...
[#] Splitting bed file into chunks for downstream processing (this might take a while) ...
[+] Found 15,632,487 BED intervals adding up to 148,668,444 b (56.88% of genome) ...
[+] BED intervals processed (3240.20MB)
[#] Generating blocks ...
[+] Analysing 91309 BED regions using 4 threads ...
[%] : 100%|█████████████████████████████████████████████████████| 91.3k/91.3k [06:45<00:00, 225it/s]
[+] Made 2,002,637 blocks covering 128,168,768 b (86.21% of BED intervals, 49.04% of genome) (5214.23MB)
[#] Generating output ...
[%] : 100%|███████████████████████████████████████████████████| 2.00M/2.00M [04:32<00:00, 7.34kit/s]
[>] Created hdf5 store: '/Users/dlaetsch/git/gIMble/hmel/gimble.blocks.h5'
[+] Total runtime: 940.337s

##### bigwig
./gIMble blocks -p hmel/ -s input/hmel.samples.csv -b input/hmel.multiinter.samples_as_string.only_intergenic.bed -g input/hmel.autosomes.genomefile -t 60
[#] Parsing parameters ...
[+] Read parameters in 0.001s (52.29MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.014s.
[+] Read 747 sequences with total length of 261,379,282 b in 0.021s
[#] Processing BED file ...
[#] Splitting bed file into chunks for downstream processing (this might take a while) ...
[+] Found 15,632,487 BED intervals adding up to 148,668,444 b (56.88% of genome) ...
[+] BED intervals processed (2559.01MB)
[#] Generating blocks ...
[+] Analysing 91309 BED regions using 60 threads ...
[%] : 100%|█████████████████████████████████████████████████████| 91.3k/91.3k [07:17<00:00, 209it/s]
[+] Made 2,002,637 blocks covering 128,168,768 b (86.21% of BED intervals, 49.04% of genome) (18220.91MB)
[#] Generating output ...
[%] : 100%|███████████████████████████████████████████████████| 2.00M/2.00M [05:07<00:00, 6.50kit/s]^[[D
[>] Created hdf5 store: '/scratch/dlaetsch/gIMble_slim/hmel/gimble.blocks.h5'
[+] Total runtime: 1154.254s

##### Compression
- hmel.multiinter.samples_as_string.only_intergenic.bed : 2.8G
- gimble.blocks.h5 : 2.3G

### bigwig+full

#### 2. [VARIANTS]

./gIMble variants -p hmel/ -s input/hmel.samples.csv -b hmel/gimble.blocks.h5 -g input/hmel.autosomes.genomefile -t 60 -v input/ros10_chi10.DP8MIN2MAC1.vcf.gz
[#] Parsing parameters ...
[+] Read parameters in 0.001s (54.55MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.014s.
[+] Read 747 sequences with total length of 261,379,282 b in 0.021s
[#] Loading blocks ...
[%] : 100%|█████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 2002637/2002637 [18:22<00:00, 1815.63it/s]
[+] Read 2,002,637 blocks covering 128,168,768 b (49.04% of genome) (58116.22MB)
[#] Fetching variants ...
[%] : 100%|███████████████████████████████████████████████████| 2.00M/2.00M [1:01:38<00:00, 541it/s]
[+] VCF parsed in 3876.593s (87158.00MB)
[#] Generating output ...
[>] Created hdf5 store: '/scratch/dlaetsch/gIMble_slim/hmel/gimble.variants.h5'
[#] Analysing variants in blocks...
[%] : 100%|█████████████████████████████████████████████████████| 2.00M/2.00M [46:32<00:00, 717it/s]
[#] Creating dataframe of global mutuple tallies...
[#] Creating dataframe of blocks...
[#] Creating dataframe for samples...

### Sample metrics: hmel/
    population_id    sample_id       blocks      sites    heterozygosity    missingness
--  ---------------  ------------  --------  ---------  ----------------  -------------
 0  chi              chi.CJ564      1603664  102634496         0.0204236     0.00779976
 1  chi              chi.CAM580     1624895  103993280         0.0196914     0.0102578
 2  chi              chi.CAM585     1572093  100613952         0.0199986     0.00949341
 3  chi              chi.CAM582     1712982  109630848         0.0200005     0.00989418
 4  chi              chi.CAM25137   1627096  104134144         0.0201602     0.00966924
 5  chi              chi.CJ560      1600673  102443072         0.0197556     0.00776403
 6  chi              chi.CJ553      1587364  101591296         0.0202169     0.00780736
 7  chi              chi.CAM25091   1566136  100232704         0.0201546     0.00980374
 8  chi              chi.CJ565      1669596  106854144         0.0206402     0.00778122
 9  chi              chi.CAM586     1621294  103762816         0.0198096     0.010188
10  ros              ros.CJ533      1329239   85071296         0.0190412     0.00559111
11  ros              ros.CAM2059    1808683  115755712         0.0198402     0.00657347
12  ros              ros.CAM1841    1678676  107435264         0.0201452     0.00645494
13  ros              ros.CAM2045    1774441  113564224         0.019802      0.00676597
14  ros              ros.CAM2519    1819802  116467328         0.0199502     0.00657092
15  ros              ros.CJ531      1391958   89085312         0.0190928     0.0055863
16  ros              ros.CAM1880    1693292  108370688         0.0200578     0.00640872
17  ros              ros.CAM2552    1452316   92948224         0.0197366     0.00675117
18  ros              ros.CJ2071     1679604  107494656         0.0199186     0.00545453
19  ros              ros.CJ546      1327917   84986688         0.0190452     0.0056218
[#] Creating dataframe for populations...

### Population metrics: hmel/
    population_id         blocks        sites    heterozygosity    missingness
--  ---------------  -----------  -----------  ----------------  -------------
 0  chi              1.61858e+06  1.03589e+08         0.0200851     0.00904587
 1  ros              1.59559e+06  1.02118e+08         0.019663      0.00617789
[#] Creating dataframe for dataset...

### Dataset metrics: hmel/
         blocks        sites         FGVs    pi_chi     pi_ros        dxy       fst
--  -----------  -----------  -----------  --------  ---------  ---------  --------
 0  2.00264e+06  1.28169e+08  0.000588693  0.017665  0.0168134  0.0241658  0.167289
[+] Total runtime: 8664.533s

#### 3. [MODIFY BLOCKS]

/gIMble modify blocks -p hmel/ -s input/hmel.samples.csv -b hmel/gimble.blocks.h5 -g input/hmel.autosomes.genomefile -t 60 -c input/hmel.chrom_coordinates.txt


#### 4. [MODIFY VARIANTS]

./gIMble modify variants -p hmel/ -v hmel/gimble.variants.h5 --max_missing 4 --max_multiallelic 64 -t 60

#### 5. [WINDOWS]

## HMEL CHR18

#### 1. [BLOCKS]

./gIMble blocks -p hmel.chr18/ -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.bed -g input/hmel.chr18.genomefile -t 4
[#] Parsing parameters ...
[+] Read parameters in 0.001s (52.89MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.010s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.013s
[#] Processing BED file ...
[#] Splitting bed file into chunks for downstream processing (this might take a while) ...
[+] Found 1,023,304 BED intervals adding up to 10,627,982 b (63.25% of genome) ...
[+] BED intervals processed (501.07MB)
[#] Generating blocks ...
[+] Analysing 5749 BED regions using 4 threads ...
[%] : 100%|█████████████████████████████████████████████████████| 5.75k/5.75k [00:13<00:00, 426it/s]
[+] Made 145,282 blocks covering 9,298,048 b (87.49% of BED intervals, 55.34% of genome) (1547.89MB)
[#] Generating output ...
[%] : 100%|█████████████████████████████████████████████████████| 145k/145k [00:10<00:00, 14.5kit/s]
[>] Created hdf5 store: '/Users/dlaetsch/git/gIMble/hmel.chr18/gimble.blocks.h5'
[+] Total runtime: 43.649s

#### 2. [VARIANTS]

./gIMble variants -p hmel.chr18/ -s input/hmel.samples.csv -b hmel.chr18/gimble.blocks.h5 -g input/hmel.chr18.genomefile -v input/hmel.chr18.vcf.gz -t 4
[#] Parsing parameters ...
[+] Read parameters in 0.001s (54.45MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.008s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.011s
[#] Loading blocks ...
[%] : 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 145282/145282 [00:30<00:00, 4695.43it/s]
[+] Read 145,282 blocks covering 9,298,048 b (55.34% of genome) (4496.64MB)
[#] Fetching variants ...
[%] : 100%|███████████████████████████████████████████████████████| 145k/145k [05:42<00:00, 425it/s]
[+] VCF parsed in 345.304s (5512.95MB)
[#] Generating output ...
[>] Created hdf5 store: '/Users/dlaetsch/git/gIMble/hmel.chr18/gimble.variants.h5'
[#] Analysing variants in blocks...
[%] : 100%|█████████████████████████████████████████████████████| 145k/145k [01:21<00:00, 1.78kit/s]
[#] Creating dataframe of global mutuple tallies...
[#] Creating dataframe of blocks...
[#] Creating dataframe for samples...

### Sample metrics: hmel.chr18/
    population_id    sample_id       blocks    sites    heterozygosity    missingness
--  ---------------  ------------  --------  -------  ----------------  -------------
 0  chi              chi.CAM25137    120722  7726208         0.0191124     0.0086015
 1  chi              chi.CAM580      120817  7732288         0.0184301     0.00914516
 2  chi              chi.CJ564       119279  7633856         0.019049      0.00707847
 3  chi              chi.CJ560       119312  7635968         0.0190792     0.00698981
 4  chi              chi.CJ553       117805  7539520         0.0179754     0.00702498
 5  chi              chi.CJ565       123764  7920896         0.0197831     0.00710033
 6  chi              chi.CAM582      126575  8100800         0.018985      0.00888752
 7  chi              chi.CAM586      120946  7740544         0.0187298     0.0090566
 8  chi              chi.CAM25091    116117  7431488         0.0186989     0.00874374
 9  chi              chi.CAM585      116865  7479360         0.0188766     0.00857079
10  ros              ros.CAM1841     123821  7924544         0.0185607     0.0057197
11  ros              ros.CAM2552     107775  6897600         0.0178275     0.00589219
12  ros              ros.CJ531       104426  6683264         0.0180502     0.00492364
13  ros              ros.CAM2045     131256  8400384         0.0186686     0.0060077
14  ros              ros.CAM2519     134260  8592640         0.0191108     0.00581998
15  ros              ros.CJ2071      124463  7965632         0.0187912     0.0048885
16  ros              ros.CAM2059     133127  8520128         0.0183754     0.00578184
17  ros              ros.CAM1880     125359  8022976         0.018798      0.00564865
18  ros              ros.CJ533        99415  6362560         0.0178471     0.00493685
19  ros              ros.CJ546        99343  6357952         0.0175486     0.00493602
[#] Creating dataframe for populations...

### Population metrics: hmel.chr18/
    population_id      blocks        sites    heterozygosity    missingness
--  ---------------  --------  -----------  ----------------  -------------
 0  chi                120220  7.69409e+06         0.018872      0.00811989
 1  ros                118324  7.57277e+06         0.0183578     0.00545551
[#] Creating dataframe for dataset...

### Dataset metrics: hmel.chr18/
      blocks        sites         FGVs     pi_chi     pi_ros        dxy     fst
--  --------  -----------  -----------  ---------  ---------  ---------  ------
 0    145282  9.29805e+06  0.000506558  0.0165783  0.0157398  0.0235341  0.1858
[+] Total runtime: 508.928s

#### 3. [MODIFY BLOCKS]

./gIMble modify blocks -p hmel.chr18/ -s input/hmel.samples.csv -c input/hmel.chrom_coordinates.txt -b hmel.chr18/gimble.blocks.h5  -g input/hmel.chr18.genomefile
[#] Parsing parameters ...
[+] Read parameters in 0.001s (52.66MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.012s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.015s
[#] Parse coordinate transformation file ...
[%] : 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 520/520 [00:00<00:00, 567830.79it/s]
[+] Read file in 0.021s (55.14MB)
[#] Loading blocks ...
[%] : 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 145282/145282 [00:31<00:00, 4645.63it/s]
[+] Read 145,282 blocks covering 9,298,048 b (55.34% of genome) (4495.08MB)
[#] Transforming coordinates ...
[%] : 100%|██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 145k/145k [00:00<00:00, 366kit/s]
[+] Transformed coordinates in 0.427s (4495.05MB)
[#] Writing output ...
[#] Generating output...
[%] : 100%|█████████████████████████████████████████████████████| 145k/145k [00:08<00:00, 16.7kit/s]
[#] 145282 of 145282 blocks (100.00%) are being retained ...
[>] Created hdf5 store: '/Users/dlaetsch/git/gIMble/hmel.chr18/gimble.blocks.modified.h5'
[+] Total runtime: 44.807s

#### 4. [MODIFY VARIANTS]

./gIMble modify variants -v hmel.chr18/gimble.variants.h5 -p hmel.chr18/ -m 4 -M 100
[+] Read parameters in 0.000s (52.64MB)
[#] Loading variants ...
[+] Excluded 197391 out of 10730431 mutuples (1.84%)
[>] Created hdf5 store: '/Users/dlaetsch/git/gIMble/hmel.chr18/gimble.variants.modified.h5'
[+] Total runtime: 26.707s

#### 5. [WINDOWS]

./gIMble windows -p hmel.chr18/ -s input/hmel.samples.csv -b hmel.chr18/gimble.blocks.modified.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18/gimble.variants.modified.h5 -w 300 -l 50 -m 1000000 -t 4

blocks: 1139Mb
variants: 3464Mb
modify blocks: 1813Mb
modify variants: 1504Mb
windows: 4196Mb

# Hmel hmel.chr18
./gIMble blocks -o hmel.chr18/ -s input/hmel.chr18.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.sorted.bed -g input/hmel.chr18.genomefile -t 4

./gIMble variants -o hmel.chr18/ -s input/hmel.chr18.samples.csv -b hmel.chr18/gimble.blocks.h5 -g input/hmel.chr18.genomefile -t 4 -v input/hmel.chr18.vcf.gz

./gIMble modify blocks -o hmel.chr18/ -s input/hmel.chr18.samples.csv -c input/hmel.chr18.chrom_coordinates.txt -b hmel.chr18/gimble.blocks.h5 -g input/hmel.chr18.genomefile

./gIMble modify variants -v hmel.chr18/gimble.variants.h5 -o hmel.chr18/ -m 4 -M 100

./gIMble windows -o hmel.chr18/ -s input/hmel.chr18.samples.csv -b hmel.chr18/gimble.blocks.modified.h5 -g input/hmel.chr18.new_coordinates.genomefile -v hmel.chr18/gimble.variants.modified.h5