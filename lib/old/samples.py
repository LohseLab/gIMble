"""usage: blocktools samples -s <FILE> [-h]
    
    -h, --help
    -s, --sample_file <FILE>                    Sample file ("sample_id,population_id")
    -o, --outprefix <STR>                       Prefix to use for sample table output                


Remarks:    
    - Sample IDs must be unique
    - Only sample IDs specified in Sample file will be considered 
        when processing BED and VCF files (others will be ignored)
    - Population and sample IDs are sorted lexicographically
"""

from docopt import docopt
from timeit import default_timer as timer
from src.classes import SamplesParameterObj
from src.functions import generate_sample_table

#def task_parse_parameters(args):
#    start = timer()
#    print("[#] Parsing parameters ...")
#    parameterObj = parse_parameters(args)
#    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
#    return parameterObj
#
#def task_load_blockDataObj(parameterObj):
#    start = timer()
#    print("[#] Loading blocks ...")
#    blockDataObj = load_blockDataObj(parameterObj)
#    print("[+] Read %s blocks in %.3fs (%.2fMB)" % (len(blockDataObj), timer() - start, memory_usage_psutil()))
#    return blockDataObj
#
#def task_fetch_variants(parameterObj, blockDataObj):
#    start = timer()
#    print("[#] Fetching variants ...")
#    blockDataObj = parse_vcf(parameterObj, blockDataObj)
#    print("[+] VCF parsed in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
#    return blockDataObj
#
#def task_analyse_variants(parameterObj, blockDataObj):
#    start = timer()
#    print("[#] Analysing variants ...")
#    blockDataObj = analyse_variants(parameterObj, blockDataObj)
#    print("[+] Analysed variants in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
#    return blockDataObj
#
#def task_write_variant_output(parameterObj, blockDataObj):
#    start = timer()
#    print("[#] Writing output ...")
#    fn_profiles_by_block_id_tsv, profileObjs_by_pair_idx = blockDataObj.write_variant_blocks(parameterObj)
#    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_profiles_by_block_id_tsv, timer() - start, memory_usage_psutil()))
#    fn_profiles_summary_tsv = blockDataObj.write_variant_summary(parameterObj, profileObjs_by_pair_idx)
#    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_profiles_summary_tsv, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    samplesParameterObj = SamplesParameterObj(args)
    generate_sample_table(samplesParameterObj)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass