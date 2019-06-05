"""
Usage: ./gIMble <module> [<args>...]

Modules:
    blocks                Makes blocks (block.bed) from multiinter BED file
    variants              Fetches and analyses variants based on VCF file and block.bed
    modify                 Relabels BED file of blocks according to coordinate file
    windows               Makes windows of blocks based on profile.tsv and    
    plot                  Makes plots
    graph                 Generate a state graph
    probs_sympy           Calculate probabilities

Options:
    -h, --help                         Show this screen.
    -v, --version                      Show version.

Help:
    https://gimble.readme.io/ (TBD)

Dependencies (via conda):
    conda create -n gimblenv -c conda-forge more-itertools tqdm scipy numpy matplotlib sympy giac networkx psutil pandas docopt
"""

'''
[chr18]
grep chr18 ../../blocktools/data/input/Hmel2.chromTransfers.2016_09_20.txt | cut -f4,5,6 > chr18.contigs.bed
grep -wFf <(cut -f1 chr18.contigs.bed) hmel.autosomes.genomefile > hmel.chr18.genomefile
bedtools intersect -b chr18.contigs.bed -a hmel.multiinter.samples_as_string.only_intergenic.bed > hmel.chr18.multiinter.samples_as_string.only_intergenic.bed
bcftools view -R chr18.contigs.bed -O z -o hmel.chr18.vcf.gz ros10_chi10.DP8MIN2MAC1.vcf.gz
tabix -p vcf hmel.chr18.vcf.gz

[ToDo]
- make output realise what a dir is and what the prefix is and based on that create output dir

'''
import sys
from docopt import docopt
from timeit import default_timer as timer

def main():
    try:
        __version__ = '0.3.0'
        start_time = timer()
        args = docopt(__doc__, version=__version__, options_first=True)
        if args['<module>'] == 'blocks':
            import cli.blocks as blocks
            blocks.main()
        elif args['<module>'] == 'variants':
            import cli.variants as variants
            variants.main()
        elif args['<module>'] == 'windows':
            import cli.windows as windows
            windows.main()
        elif args['<module>'] == 'modify':
            import cli.modify as modify
            modify.main()
        elif args['<module>'] == 'fixcoordinates':
            import cli.fixcoordinates as fixcoordinates
            fixcoordinates.main()
        elif args['<module>'] == 'graph':
            import cli.graph as graph
            graph.main()
        elif args['<module>'] == 'probs_sympy':
            import lib.probs_sympy as probs_sympy
            probs_sympy.main()
        else:
            sys.exit("%r is not a gIMble module. See 'gIMble -help'." % args['<module>'])
    except KeyboardInterrupt:
        sys.stderr.write("\n[X] Interrupted by user after %i seconds!\n" % (timer() - start_time))
        sys.exit(-1)
