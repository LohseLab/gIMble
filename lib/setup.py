class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.module = 'setup'
        self.vcf_file = args['--vcf']
        self.bed_file = args['--bed']
        self.genome_file = args['--genome']
        self.sample_file = args['--sample'] 
        self.outprefix = args['--outprefix']

        self.block_length = int(args['--block_length'])
        self.block_span = int(args['--block_span'])
        self.block_gap_run = int(args['--block_gap_run'])

        self.pairedness = int(args['--pairedness'])