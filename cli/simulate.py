#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""usage: gIMble simulate                   -z FILE
                                            (-m FILE | -w FILE) 
                                            [-c FILE] 
                                            [-b INT] [-r INT]
                                            [-t INT] [-h|--help]
                                            
    Options:
        -h --help                                   show this
        -m, --model_file FILE                       Model file to analyse
        -w, --windows_file FILE                     Output path for tsv file
        -z, --zarr FILE                             Path to zarr store
        -c, --config_file FILE                      Config file with parameters (if not present, empty config file is created)
        -P, --precision INT                         Floating point precision of probabilities [default: 25]
        -b, --blocks INT                            Number of blocks per replicate
        -r, --replicates INT                        Number of replicates per parametercombo
        -t, --threads INT                           Threads [default: 1]
        
"""
import pathlib
import oyaml as yaml
import collections
from timeit import default_timer as timer
from docopt import docopt
import sys, os
from lib.gimble import RunObj
from lib.gimble import Store
import lib.simulate
import numpy as np
import zarr
import pandas as pd

"""
test command:
gIMble simulate -m 
./gIMble model -s A,B -p 2 -n 1,1 -m 'A>B' -j 'A,B' -o output/test
./gIMble model -s A,B -p 2 -n 2,1 -j 'A,B' -o output/test

./gIMble simulate -m output/s_A_B.p2.n_1_1.m_AtoB.j_A_B.model.tsv
./gIMble simulate -m output/test.model.tsv

./gIMble simulate -m output/s_A_B.p2.n_1_1.m_AtoB.j_A_B.model.tsv -c /Users/s1854903/git/gIMble/output/s_A_B.p2.n_1_1.m_AtoB.j_A_B.model.config.yaml -z /Users/s1854903/git/gIMble/test.z
./gIMble simulate -z output/sims.z -m output/test.model.tsv -c output/test.model.config.yaml
"""


class ParameterObj(RunObj):
    """Sanitises command line arguments and stores parameters."""

    def __init__(self, params, args):
        super().__init__(params)
        if args["--windows_file"]:
            self.windows_path = args["--windows_file"]
            self._get_zarr_store(args)
            store = lib.gimble.load_store(self)
            self.return_windows_tsv(store)
            sys.exit("[X] generated windows BED file")
        self.model_file = self._get_path(args["--model_file"])
        self.config_file = self._get_path(args["--config_file"])
        self.threads = self._get_int(args["--threads"])
        self._config = self._get_or_write_config(args["--blocks"], args["--replicates"])
        self.data_type = "simulations"
        self._get_zarr_store(args)
        self.parameter_grid = None

    def _get_zarr_store(self, args):
        z = args['--zarr']
        if z:
            #check if file exists
            if z.endswith('.z'):
                self.path = z
                z_path = pathlib.Path(z)
                if z_path.is_dir() or z_path.is_file():
                    self.zstore = z
                else:
                    self.zstore= z
                    self.data = zarr.open(self.path, mode='w')
            else:
                sys.exit("[X] Specify the path to a zarr file ending in .z .")
        else:
            sys.exit("[X] Specify path for zarr file.")

    def _get_or_write_config(self, blocks, replicates):
        #in case no config file is provided
        if self.config_file is None:
            print("[-] No config file found.")
            print("[+] Generating config file for model %r" % self.model_file)
            """for now we use the following dict until columns are fixed in gimble model"""
            if blocks is None:
                blocks = 1
            if replicates is None:
                replicates = 1
            config = {
                "version": self._VERSION,
                "model": self.model_file.as_posix(),
                # "random_seed": 12345,
                "precision": 25,
                "blocks": int(blocks),
                "blocklength": 1000,
                "replicates": int(replicates),
                #'k_max': collections.defaultdict(dict),
                "parameters": collections.defaultdict(dict),
                "boundaries": collections.defaultdict(list),
                "recombination":collections.defaultdict(dict),
                "grid":collections.defaultdict(dict)
            }
            (pop_configs, columns) = self._parse_model_file()
            config["ploidy"] = int(pop_configs["ploidy"])
            config["parameters"]["sample_size_A"] = pop_configs["A"]
            config["parameters"]["sample_size_B"] = pop_configs["B"]
            config["parameters"]["theta"] = "FLOAT"
            for column in columns:
                if column.startswith("C_") or column.startswith("M_"):
                    config["parameters"][column] = "FLOAT"
            config["parameters"]["T"] = "FLOAT"
            config["recombination"]["rate"] = "FLOAT"
            config["recombination"]["recombination_map"] = "FILE_PATH"
            config["recombination"]["cutoff"] = "INT_0_100"
            config["recombination"]["number_bins"] = "INT"
            config["recombination"]["scale"]="LINEAR/LOG"
            config["grid"]["file"] = "FILE_PATH"
            config["grid"]["parameters"] = "PARAMETERS"
            for parameter in config["parameters"]:
                if parameter not in ["sample_size_A", "sample_size_B"]:
                    config["boundaries"][parameter] = ["MIN", "MAX", "STEPSIZE"]
            config_file = pathlib.Path(self.model_file).with_suffix(".config.yaml")
            yaml.add_representer(
                collections.defaultdict, yaml.representer.Representer.represent_dict
            )
            with open(config_file, "w") as fh:
                yaml.dump(config, fh)
            print("[+] Wrote file %r" % str(config_file))
            sys.exit(
                "[X] Please specify parameters in config file %r" % str(config_file)
            )

        #if a config file is given
        else:
            print("[+] Reading config %r" % self.config_file)
            config_raw = yaml.safe_load(open(self.config_file, "r"))
            config = {}
            for k in ['version', 'model', 'precision']:
                config[k]=config_raw[k]
            for k in ['blocks', 'blocklength', 'replicates', 'ploidy']:
                assert isinstance(config_raw[k], int), f"integer value required for {k}"
                config[k]=int(config_raw[k])
            config["parameters"] = {}

            #grid
            #check whether file_path is path
            p_grid_fpath = config_raw['grid']['file']
            p_grid_names = []
            if os.path.isfile(p_grid_fpath):
                p_grid_names = config_raw['grid']['parameters']
                self.parameter_grid = pd.read_csv(p_grid_fpath, 
                    names=p_grid_names, index=False, header=None, sep='\t')
                
            
            #boundaries and parameters
            for key, value in config_raw['boundaries'].items():
                if any(isinstance(v, str) for v in value):
                    if not key in p_grid_names:
                        assert not isinstance(config_raw['parameters'][key], str), f"value required for parameter {k}"
                        config['parameters'][key] = [config_raw['parameters'][key],]
                else:
                    config['parameters'][key] = value
                    center = config_raw['parameters'][key]
                    if not isinstance(center,str):
                        config['parameters'][key].append(center)
            
            #verify there is no overlap between parameters and gridfile
            assert len(config_raw['boundaries']) == len(p_grid_names)+len(config['parameters']), "The same parameter has been specified both in the gridfile and in the yamlfile."
            
            #sample size
            config["parameters"]["sample_size_A"] = [config_raw["parameters"]["sample_size_A"],]
            config["parameters"]["sample_size_B"] = [config_raw["parameters"]["sample_size_B"],]

            #recombination
            file_path = config_raw['recombination']['recombination_map']
            if file_path != 'FILE_PATH':
                file_path = pathlib.Path(file_path)
                if file_path.is_file():
                    cutoff = config_raw['recombination']['cutoff']
                    rbins = config_raw['recombination']['number_bins']
                    if isinstance(cutoff, str):
                        cutoff = 90
                    if isinstance(rbins, str):
                        rbins = 10
                    if isinstance(config_raw['recombination']['rate'], float):
                        sys.exit("[X] Either provide a recombination map or a rate")
                    scale = config_raw['recombination']['scale']
                    window_bin, to_be_simulated = self._parse_recombination_map(file_path.as_posix(), cutoff, rbins, scale)
                    #at the moment window_bin info is not used but should be
                    #indicates to which bin each window belongs
                    config["parameters"]["recombination"] = to_be_simulated
                else:
                    sys.exit("[X] Incorrect path to recombination map")
            else:
                rec_rate = config_raw['recombination']['rate']
                if isinstance(rec_rate, float):
                    config["parameters"]["recombination"] = [rec_rate,]
                else:
                    config["parameters"]["recombination"] = [0.0,]

            (pop_configs, columns) = self._parse_model_file()
            assert (
                config["parameters"]["sample_size_A"][0] == pop_configs["A"]
                and config["parameters"]["sample_size_B"][0] == pop_configs["B"]
            ), "sample size does not match model sample size"
            return config

    def _parse_model_file(self):
        with open(self.model_file) as fh:
            first_line = fh.readline().rstrip()
            columns = first_line.split()
            second_line = fh.readline().rstrip().split()
            A, B = second_line[3].split(";")
            pop_configs = {}
            for pop, name in zip([A, B], ["A", "B"]):
                pop = pop.lstrip(f"{name}=[").rstrip("]")
                unique_el, count_el = np.unique(pop.split(","), return_counts=True)
                pop_configs[name] = len(unique_el)
            pop_configs["ploidy"] = count_el[0]
            return (pop_configs, columns)

    def _parse_recombination_map(self, v, cutoff, bins, scale):
        #load bedfile
        hapmap = pd.read_csv(v, sep='\t', 
            names=['chrom', 'chromStart', 'chromEnd', 'rec'], header=None)
        #from cM/Mb to rec/bp
        hapmap['rec_scaled'] = hapmap['rec']*1e-8
        return self._make_bins(hapmap['rec_scaled'], scale, cutoff, bins)

    def _make_bins(self, x, scale, cutoff=90, bins=10):
        x=np.array(x)
        clip_value = np.percentile(x,cutoff)
        clip_array = np.clip(x, a_min=None, a_max=clip_value)
        with_zero = np.count_nonzero(x==0.0)
        if scale.lower() == "log":
            start, stop = np.min(clip_array), np.max(clip_array)
            if with_zero:
                zero_index = np.where(clip_array==0.0)
                start=np.min(np.delete(clip_array, zero_index))
            start, stop = np.log10(start), np.log10(stop)
            bin_edges = np.logspace(start, stop, num=bins)
            counts, bin_edges = np.histogram(clip_array, bins=bin_edges)
        
        elif scale.lower() == "linear": 
            counts, bin_edges = np.histogram(clip_array, bins=bins)
        else:
            sys.exit("[-] Scale of recombination values to be simulated should either be LINEAR or LOG")
        bin_with_counts = counts>0
        no_counts = bins-sum(1 for v in bin_with_counts if v)
        print(f"[+] There are {no_counts} recombination bins without counts")
        to_be_simulated  = [(stop + start)/2 for start, stop, c in zip(bin_edges[:-1],bin_edges[1:], bin_with_counts) if c]
        if with_zero:
            to_be_simulated = [0] + list(to_be_simulated)
        #assign each window to a specific bin
        #needs to be adjusted to the actual bins that are simulated
        #if there are windows with rec rate zero these will be simulated as well.
        window_bin = np.digitize(clip_array, bin_edges, right=True)
        return (window_bin, list(to_be_simulated))

    def _check_recombination_map(self, store, df):
        starts = store.data[chrom]['windows/starts'][:]
        ends = store.data[chrom]['windows/ends'][:]
        if set(starts) != set(df['starts']):
            sys.exit("[X] Starts recombination map do not match window coordinates")
        if set(ends) != set(df['ends']):
            sys.exit("[X] Ends recombination map do not match window coordinates")

    def return_windows_tsv(self, store):
        df_list = [] 
        path = self.windows_path
        for chrom in store.data.group_keys():
            if chrom not in ['sims', 'simulations']:
                #chromosome
                starts = store.data[chrom]['windows/starts'][:]
                ends = store.data[chrom]['windows/ends'][:]
                d = [(chrom, int(start), int(end)) for start, end in zip(starts, ends)]
                d = set(d)
                df_small = pd.DataFrame(d, columns=['chrom', 'starts', 'ends'])
                df_small=df_small.sort_values(by=['starts', 'ends'])
                df_list.append(df_small)
        df = pd.concat(df_list)
        df.to_csv(path, index=False , sep='\t')

    def simulate(self):
        replicate = lib.simulate.run_sim(self)
        lib.simulate.get_genotypes(self, replicate)


def main(params):
    try:
        start_time = timer()
        args = docopt(__doc__)
        parameterObj = ParameterObj(params, args)
        lib.simulate.run_sim(parameterObj)

        print("[*] Total runtime: %.3fs" % (timer() - start_time))
    except KeyboardInterrupt:
        print("\n[X] Interrupted by user after %s seconds!\n" % (timer() - start_time))
        exit(-1)
