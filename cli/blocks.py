"""usage: blocktools blocks -s <FILE> -b <FILE> -g <FILE> [-p <STR> -x <INT> -l <INT> -m <INT> -a <INT> -t <INT> -h]
    
    -h, --help
    -x, --block_length <INT>                    Block length in bases [default: 64]
    -l, --max_block_length <INT>                Maximum length of block from start to end [default: 80]
    -m, --min_interval_length <INT>             Minimum length in bases of a BED interval to be considered [default: 0]
    -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
    -b, --bed_file <FILE>                       multiintersectBED file of 
    -g, --genome_file <FILE>                    Genome file (as used in BedTools)
    -t, --threads <INT>                         Number of threads to use [default: 1]
    -p, --prefix <STR>                          Folder for output
"""

'''
[To Do]
- output distribution of distance between blocks
- output barchart of bases blocked per sample
- output bases blocked vs per sequence

-a, --max_interval_distance <INT>           Maximum distance in bases between BED intervals to be considered [default: 20] 
=> make internal only, no need to expose to user
- gaps

./gIMble blocks -s input/hmel.samples.csv -b input/hmel.chr18.multiinter.samples_as_string.only_intergenic.bed -g input/hmel.chr18.genomefile -p hmel.chr18
[#] Parsing parameters ...
[+] Read parameters in 0.001s (52.77MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.008s.
[+] Read 19 sequences with total length of 16,802,090 b in 0.010s
[#] Processing BED file ...
[#] Splitting bed file into chunks for downstream processing (this might take a while) ...
[+] Found 1,023,304 BED intervals adding up to 10,627,982 b (63.25% of genome) ...
[+] BED intervals processed (500.61MB)
[#] Generating blocks ...
[+] Analysing 5749 BED regions using 1 threads ...
[%] : 100%|█████████████████████████████████████████████████████| 5.75k/5.75k [00:21<00:00, 266it/s]
[+] Made 145,282 blocks covering 9,298,048 b (87.49% of BED intervals, 55.34% of genome) (1303.26MB)
[#] Generating output ...
[%] : 100%|█████████████████████████████████████████████████████| 145k/145k [00:08<00:00, 17.0kit/s]
[>] Created: '/Users/dlaetsch/git/gIMble/hmel.chr18.blocks.bed'
[>] Created: '/Users/dlaetsch/git/gIMble/hmel.chr18.blocks.tsv'
[+] Total runtime: 46.311s
'''
from timeit import default_timer as timer
from docopt import docopt
import lib.blocks

def main():
    start_time = timer()
    args = docopt(__doc__)
    blockParameterObj = lib.blocks.task_generate_parameterObj(args)
    entityCollection = lib.blocks.task_generate_entityCollection(blockParameterObj)
    region_dfs = lib.blocks.task_generate_region_dfs(blockParameterObj, entityCollection)
    lib.blocks.task_make_blocks(blockParameterObj, region_dfs, entityCollection)
    lib.blocks.task_write_block_output(blockParameterObj, entityCollection)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass