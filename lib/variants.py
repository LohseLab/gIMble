from timeit import default_timer as timer
from lib.functions import check_file, memory_usage_psutil, poolcontext, check_path, check_prefix, format_count, format_percentage, format_bases
from lib.classes import EntityCollection
import collections

from tqdm import tqdm
import pysam

'''
[To Do]

'''

NUCLEOTIDES = set(['A', 'G', 'C', 'T'])

MUTYPE_BY_ZYGOSITY = {
    'HOM' : {
        'HOM' : 'fixed', 
        'HET' : 'hetB',
        'MISS': 'missing'
    },
    'HET' : {
        'HOM' : 'hetA',
        'HET' : 'hetAB',
        'MISS': 'missing'
    },
    'MISS': {
        'HOM' : 'missing',
        'HET' : 'missing',
        'MISS': 'missing'
    }
}

class ParameterObj(object):
    def __init__(self, args):
        self.args = args
        self.sample_file = check_file(args.get('--sample_file', None))
        self.genome_file = check_file(args.get('--genome_file', None))
        self.blocks_file = check_file(args.get('--blocks', None))
        self.vcf_file = check_file(args.get('--vcf', None))
        self.path = check_path(args.get('--prefix', None))
        self.prefix = check_prefix(args.get('--prefix', None))
        self.dataset = self.prefix if self.path is None else (self.path / self.prefix) 
        self.threads = int(args.get('--threads'))
        self.block_length = None
        self.max_block_length = None
        self.variant_count = 0

def gt_zygosity(gt):
    gt_set = set(gt)
    if None in gt_set:
        return 'MISS'
    if len(gt_set) == 1:
        return 'HOM'
    return 'HET'

def parse_vcf(parameterObj, entityCollection):
    params = []
    for blockObj in entityCollection.blockObjs:
        params.append((blockObj, entityCollection.pairObjs, parameterObj.vcf_file))
    if parameterObj.threads < 2:
        with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
            for param in params:
                blockObj = infer_mutypes(param)
                block_idx = entityCollection.block_idx_by_block_id[blockObj.id]
                entityCollection.blockObjs[block_idx] = blockObj
                pbar.update()
    else:
        with poolcontext(processes=parameterObj.threads) as pool:
            with tqdm(total=len(params), desc="[%] ", ncols=100, unit_scale=True) as pbar:
                for blockObj in pool.imap_unordered(infer_mutypes, params):    
                    block_idx = entityCollection.block_idx_by_block_id[blockObj.id]
                    entityCollection.blockObjs[block_idx] = blockObj
                    pbar.update()
    return entityCollection

def infer_mutypes(param):
    blockObj, pairObjs, vcf_file = param
    genotypes_by_sample_id = collections.defaultdict(list)
    vcf_reader = pysam.VariantFile(vcf_file, duplicate_filehandle=True)
    # should speed it up. a bit
    if blockObj.length == blockObj.span:
        for record in vcf_reader.fetch(contig=blockObj.sequence_id, start=blockObj.start, stop=blockObj.end):
            #print("\t".join(["%s=%s" % (sample, record.samples[sample]['GT']) for sample in blockObj.sample_ids]))
            try:
                if all([alt in NUCLEOTIDES for alt in record.alts]): # is SNP
                    for sample_id in blockObj.sample_ids:
                        genotypes_by_sample_id[sample_id].append([
                            record.samples[sample_id]['GT'][0], 
                            record.samples[sample_id]['GT'][1]
                        ])
            except TypeError: # monomorphic
                pass
    else:
        for contig_id, start, end in blockObj.bed_tuples:
            #print("# %s %s %s" % (contig_id, start, end))
            for record in vcf_reader.fetch(contig=contig_id, start=start, stop=end):
                #print("\t".join(["%s=%s" % (sample, record.samples[sample]['GT']) for sample in blockObj.sample_ids]))
                try:
                    if all([alt in NUCLEOTIDES for alt in record.alts]): # is SNP
                        for sample_id in blockObj.sample_ids:
                            genotypes_by_sample_id[sample_id].append([
                                record.samples[sample_id]['GT'][0], 
                                record.samples[sample_id]['GT'][1]
                            ])
                except TypeError: # monomorphic
                    pass
    # get zygosity for each sample
    zygosities_by_sample_id = collections.defaultdict(list)
    for sample_id, genotypes in genotypes_by_sample_id.items():
        for genotype in genotypes:
            zygosity = gt_zygosity(genotype)
            zygosities_by_sample_id[sample_id].append(zygosity)
        blockObj.gt_counter_by_sample_id[sample_id] = collections.Counter(zygosities_by_sample_id[sample_id])
        blockObj.gt_counter_by_sample_id[sample_id]['blocks'] += 1
    # get mutype for each pair
    for pair_idx in blockObj.pair_idxs:
        mutypes = []
        #print("## pair %s : %s" % (pair_idx, pairObjs[pair_idx].id))
        sample_id_1, sample_id_2 = pairObjs[pair_idx].id
        for genotype_1, zygosity_1, genotype_2, zygosity_2 in zip(
                genotypes_by_sample_id[sample_id_1], zygosities_by_sample_id[sample_id_1],
                genotypes_by_sample_id[sample_id_2], zygosities_by_sample_id[sample_id_2]):
            #print(sample_id_1, sample_id_2, ":", genotype_1, zygosity_1, genotype_2, zygosity_2)
            genotype_set = set(genotype_1 + genotype_2)
            if len(genotype_set) > 2:
                mutypes.append('multiallelic')
            elif len(genotype_set) == 2:
                mutypes.append(MUTYPE_BY_ZYGOSITY[zygosity_1][zygosity_2])
            else:
                pass
        #print(collections.Counter(mutypes))
        blockObj.mutype_counter_by_pair_idx[pair_idx] = collections.Counter(mutypes)
    return blockObj

# def infer_mutypes_cyvcf(param):
#     from cyvcf2 import VCF
#     blockObj, pairObjs, vcf_file = param
#     genotype_by_sample_id = collections.defaultdict(list)    
#     vcf_reader = VCF(vcf_file, samples=list(blockObj.sample_ids))       
#     for contig_id, start, end in blockObj.bed_tuples:        
#         for record in vcf_reader('%s:%s-%s' % (contig_id, start, end)):            
#             if record.is_snp:                
#                 for sample_id, genotype in zip(vcf_reader.samples, record.genotypes):                    
#                     genotype_by_sample_id[sample_id].append([genotype[0], genotype[1]])
#     mutypes_by_pair_idx = {}
#     for pair_idx in blockObj.pair_idxs:
#         mutypes = []
#         sample_id_1, sample_id_2 = pairObjs[pair_idx].id        
#         for genotype_1, genotype_2 in zip(genotype_by_sample_id[sample_id_1], genotype_by_sample_id[sample_id_2]):            
#             genotype_set = set(genotype_1 + genotype_2)
#             if len(genotype_set) > 2:
#                 mutypes.append('multiallelic')
#             elif -1 in genotype_set:
#                 mutypes.append('missing')
#             else:
#                 zygosity_1 = gt_zygosity(genotype_1)
#                 zygosity_2 = gt_zygosity(genotype_2)
#                 mutypes.append(MUTYPE_BY_ZYGOSITY[zygosity_1][zygosity_2])
#         mutypes_by_pair_idx[pair_idx] = collections.Counter(mutypes)
#     return blockObj.id, mutypes_by_pair_idx

def task_generate_entityCollection(parameterObj):
    start = timer()
    print("[#] Building entities based on samples and sequences...")
    entityCollection = EntityCollection()
    entityCollection.parse_sample_file(parameterObj)
    print("[+] Read %s samples from %s populations and generated %s pairs in %.3fs." % (\
        entityCollection.count('samples'), 
        entityCollection.count('populations'), 
        entityCollection.count('pairs'), 
        timer() - start))
    entityCollection.parse_genome_file(parameterObj)
    print("[+] Read %s sequences with total length of %s in %.3fs" % (
        entityCollection.count('sequences'), 
        format_bases(entityCollection.count('bases')), 
        timer() - start))
    return entityCollection

def task_parse_parameters(args):
    start = timer()
    print("[#] Parsing parameters ...")
    parameterObj = ParameterObj(args)
    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return parameterObj

def task_load_blockObjs(parameterObj, entityCollection):
    print("[#] Loading blocks ...")
    entityCollection.load_blocks(parameterObj)
    block_count = format_count(entityCollection.count('blocks'))
    blocked_bases = entityCollection.count('blocks') * parameterObj.block_length
    total_bases_percentage = format_percentage(blocked_bases / entityCollection.count('bases')) 
    print("[+] Read %s blocks covering %s (%s of genome) (%.2fMB)" % (
        block_count, 
        format_bases(blocked_bases),
        total_bases_percentage,
        memory_usage_psutil()
        ))

def task_infer_mutypes(parameterObj, entityCollection):
    start = timer()
    print("[#] Fetching variants ...")
    entityCollection = parse_vcf(parameterObj, entityCollection)
    print("[+] VCF parsed in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))

def task_write_variant_output(parameterObj, entityCollection):
    print("[#] Generating output ...")
    entityCollection.generate_variant_output(parameterObj)
   