"""usage: blocktools variants -v <FILE> -b <FILE> -s <FILE> -g <FILE> [-p <STR> -t <INT> -h]
    
    -h, --help
    -v, --vcf <FILE>                            VCF file (tabix'ed)
    -b, --blocks <FILE>                         blocks HDF5 file 
    -p, --prefix <STR>                          Folder for output
    -t, --threads <FILE>                        Threads [default: 1]
    -g, --genome_file <FILE>                    Genome file as used in BedTools
    -s, --sample_file <FILE>                    CSV file ("sample_id,population_id")
"""

'''

 [ To Do ]
 - Variants: filter blocks on 4 missing both for global counts and window counts


- population variants: 
    - collect gt_counts and mean



./gIMble variants -v input/hmel.chr18.vcf.gz -b hmel.chr18/gimble.blocks.bed -B hmel.chr18/gimble.blocks.tsv -s input/hmel.samples.csv -g input/hmel.chr18.genomefile -t 4 -p hmel.chr18/
[#] Parsing parameters ...
[+] Read parameters in 0.001s (54.22MB)
[#] Building entities based on samples and sequences...
[+] Read 20 samples from 2 populations and generated 100 pairs in 0.008s.
[+] Read 19 sequences with total length of 16802090 b in 0.010s
[#] Loading blocks ...
[%] : 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 145282/145282 [00:30<00:00, 4703.12it/s]
[+] Read 145282 blocks in 33.268s (4390.63MB)
[#] Fetching variants ...
[%] : 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 145k/145k [07:14<00:00, 334it/s]
[+] VCF parsed in 438.443s (5029.77MB)
[#] Generating output ...
[%] : 100%|█████████████████████████████████████████████████████| 145k/145k [00:41<00:00, 3.51kit/s]
[>] Created: '/Users/dlaetsch/git/gIMble/hmel.chr18/gimble.block_variants.tsv'
[>] Created: '/Users/dlaetsch/git/gIMble/hmel.chr18/gimble.pair_variants.tsv'
Sample=chi.CAM25091 total=7,583,168 b sampled=(total-miss)=7,518,189 b het=138961 hom=705873 miss=64979 het_rate=(het/sampled)=1.85%
Sample=chi.CAM25137 total=7,880,832 b sampled=(total-miss)=7,814,375 b het=147666 hom=733079 miss=66457 het_rate=(het/sampled)=1.89%
Sample=chi.CAM580 total=7,888,192 b sampled=(total-miss)=7,817,479 b het=142507 hom=729888 miss=70713 het_rate=(het/sampled)=1.82%
Sample=chi.CAM582 total=8,263,808 b sampled=(total-miss)=8,191,812 b het=153794 hom=761760 miss=71996 het_rate=(het/sampled)=1.88%
Sample=chi.CAM585 total=7,632,320 b sampled=(total-miss)=7,568,216 b het=141185 hom=707119 miss=64104 het_rate=(het/sampled)=1.87%
Sample=chi.CAM586 total=7,895,552 b sampled=(total-miss)=7,825,449 b het=144979 hom=730797 miss=70103 het_rate=(het/sampled)=1.85%
Sample=chi.CJ553 total=7,696,832 b sampled=(total-miss)=7,643,867 b het=135526 hom=722270 miss=52965 het_rate=(het/sampled)=1.77%
Sample=chi.CJ560 total=7,795,264 b sampled=(total-miss)=7,741,890 b het=145688 hom=723118 miss=53374 het_rate=(het/sampled)=1.88%
Sample=chi.CJ564 total=7,793,024 b sampled=(total-miss)=7,738,988 b het=145417 hom=723221 miss=54036 het_rate=(het/sampled)=1.88%
Sample=chi.CJ565 total=8,085,184 b sampled=(total-miss)=8,028,943 b het=156700 hom=746348 miss=56241 het_rate=(het/sampled)=1.95%
Sample=ros.CAM1841 total=8,085,184 b sampled=(total-miss)=8,039,858 b het=147085 hom=776104 miss=45326 het_rate=(het/sampled)=1.83%
Sample=ros.CAM1880 total=8,183,872 b sampled=(total-miss)=8,138,553 b het=150816 hom=786377 miss=45319 het_rate=(het/sampled)=1.85%
Sample=ros.CAM2045 total=8,569,984 b sampled=(total-miss)=8,519,517 b het=156823 hom=813344 miss=50467 het_rate=(het/sampled)=1.84%
Sample=ros.CAM2059 total=8,693,312 b sampled=(total-miss)=8,644,050 b het=156561 hom=830916 miss=49262 het_rate=(het/sampled)=1.81%
Sample=ros.CAM2519 total=8,765,760 b sampled=(total-miss)=8,715,751 b het=164212 hom=831931 miss=50009 het_rate=(het/sampled)=1.88%
Sample=ros.CAM2552 total=7,043,904 b sampled=(total-miss)=7,003,262 b het=122967 hom=679481 miss=40642 het_rate=(het/sampled)=1.76%
Sample=ros.CJ2071 total=8,129,600 b sampled=(total-miss)=8,090,660 b het=149684 hom=776972 miss=38940 het_rate=(het/sampled)=1.85%
Sample=ros.CJ531 total=6,827,648 b sampled=(total-miss)=6,794,742 b het=120634 hom=652063 miss=32906 het_rate=(het/sampled)=1.78%
Sample=ros.CJ533 total=6,501,632 b sampled=(total-miss)=6,470,221 b het=113553 hom=620007 miss=31411 het_rate=(het/sampled)=1.76%
Sample=ros.CJ546 total=6,499,008 b sampled=(total-miss)=6,467,625 b het=111573 hom=620669 miss=31383 het_rate=(het/sampled)=1.73%
[+] Total runtime: 570.308s
'''
from docopt import docopt
from timeit import default_timer as timer
import lib.variants

def main():
    start_time = timer()
    args = docopt(__doc__)
    parameterObj = lib.variants.task_parse_parameters(args)
    entityCollection = lib.variants.task_generate_entityCollection(parameterObj)
    #print(entityCollection)
    lib.variants.task_load_blockObjs(parameterObj, entityCollection)
    lib.variants.task_infer_mutypes(parameterObj, entityCollection)
    lib.variants.task_write_variant_output(parameterObj, entityCollection)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass