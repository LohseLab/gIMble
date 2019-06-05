"""usage: blocktools plot -y <FILE> [-g <FILE> -t <INT> -p <STR> -1 -2 -3 -4 -h]
    
    -h, --help
    -y, --yaml <FILE>                       YAML file with parameters
    -t, --threads <INT>                     Number of threads to use
    -g, --genome <FILE>                     Alternative genomefile
    -p, --outprefix <STR>                   Alternative outprefix
    -1                                      plot_coverage
    -2                                      plot_window_coverage
    -3                                      plot_variant_pairs
    -4                                      plot_window_variant
"""

from docopt import docopt
from lib.src.code import *
from timeit import default_timer as timer

def task_parse_parameters(args):
    start = timer()
    print("[#] Parsing parameters ...")
    parameterObj = parse_parameters(args)
    print("[+] Read parameters in %.3fs (%.2fMB)" % (timer() - start, memory_usage_psutil()))
    return parameterObj

def task_parse_genome_f(parameterObj):
    start = timer()
    print("[#] Parsing genome file ...")
    sequence_OrdDict = parse_genome_f(parameterObj)
    print("[+] Read %s sequences in %.3fs (%.2fMB)" % (len(sequence_OrdDict), timer() - start, memory_usage_psutil()))
    return sequence_OrdDict

def task_plot_coverage(parameterObj, sequence_OrdDict):
    start = timer()
    print("[#] Parsing %s ..." % parameterObj.block_pairs_f)
    fn_block_span_tsv = plot_coverage(parameterObj, sequence_OrdDict)
    print("[+] Wrote '%s' in %.3fs (%.2fMB)" % (fn_block_span_tsv, timer() - start, memory_usage_psutil()))

def task_plot_window_coverage_tsv(parameterObj, sequence_OrdDict):
    start = timer()
    print("[#] Parsing %s ..." % parameterObj.window_coverage_tsv_f)
    fn_mean_block_density_hist, fn_sample_cov_genome = plot_window_coverage_tsv(parameterObj, sequence_OrdDict)
    print("[+] Wrote \n\t'%s' and \n\t'%s' in %.3fs (%.2fMB)" % (fn_mean_block_density_hist, fn_sample_cov_genome, timer() - start, memory_usage_psutil()))

def task_plot_variant_pairs_tsv(parameterObj, sequence_OrdDict):
    start = timer()
    print("[#] Parsing %s ..." % parameterObj.variant_pairs_tsv_f)
    piA_fn, piB_fn, dxy_fn, fst_fn = plot_variant_pairs_tsv(parameterObj, sequence_OrdDict)
    print("[+] Wrote \n\t'%s' and \n\t'%s' and \n\t'%s' and \n\t'%s' in %.3fs (%.2fMB)" % (piA_fn, piB_fn, dxy_fn, fst_fn, timer() - start, memory_usage_psutil()))

def task_plot_window_variant_tsv(parameterObj, sequence_OrdDict):
    start = timer()
    print("[#] Parsing %s ..." % parameterObj.window_variant_tsv_f)
    dxy_fst_fn, piA_piB_fn, tuple_fn = plot_window_variant_tsv(parameterObj, sequence_OrdDict)
    print("[+] Wrote \n\t'%s' and \n\t'%s' and \n\t'%s' in %.3fs (%.2fMB)" % (dxy_fst_fn, piA_piB_fn, tuple_fn, timer() - start, memory_usage_psutil()))

def main():
    start_time = timer()
    args = docopt(__doc__)
    parameterObj = task_parse_parameters(args)
    sequence_OrdDict = task_parse_genome_f(parameterObj)
    plots = { \
        '-1': task_plot_coverage, \
        '-2': task_plot_window_coverage_tsv, \
        '-3': task_plot_variant_pairs_tsv, \
        '-4': task_plot_window_variant_tsv \
        }
    commands = []
    for key, command in plots.items():
        if key in args:
            commands.append(plots[key])
    if not commands:
        commands = [command for command in plots.values()]
    for command in commands:
        
        command(parameterObj, sequence_OrdDict)
    print("[+] Total runtime: %.3fs" % (timer() - start_time))

if __name__ == "__main__":
    pass