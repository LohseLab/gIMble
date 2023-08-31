import itertools
from tqdm import tqdm
from tqdm.dask import TqdmCallback
import allel
import demes
import math
import numpy as np
import dask
import pandas as pd
import shutil
import zarr
import numcodecs
import os
import string
import collections
import sys
import warnings
import pathlib
import contextlib
import multiprocessing
import lib.simulate
import fractions
import copy
import tabulate
from timeit import default_timer as timer
import agemo
import lib.optimize
import msprime
import pkg_resources

# np.set_printoptions(threshold=sys.maxsize)
#from lib.GeneratingFunction.gf import togimble

"""
- Coordinate systems:
    {GIMBLE}                    : 0-based       | {GIMBLE}  
    -----------------------------------------------------------
    {GENOME_FILE} : LENGTH      : 0-based   =>  | [0, l)
    {BED_FILE}    : START:END   : 0-based   =>  | [START, END)
    {VCF_FILE}    : POS         : 1-based   =>  | [POS-1)
"""

###
RECOMBINATION_SCALER = 1e-8
MUTYPES = ["m_1", "m_2", "m_3", "m_4"]
SPACING = 16
GRIDSEARCH_DTYPE = np.float32  # -3.4028235e+38 ... 3.4028235e+38
AGEMO_VERSION = pkg_resources.get_distribution("agemo").version
MIN_AGEMO_VERSION = '0.0.2'
# GRIDSEARCH_DTYPE=np.float64 # -1.7976931348623157e+308 ... 1.7976931348623157e+308
###

DFRM = "[──"
DFRT = "├──"
DPPM = "    "
DFRL = "└──"
DPPL = "│   "

SIGNS = {
    "T": "%s%s%s" % ("\u251C", "\u2500", "\u2500"),  # '├──'
    "S": "%s%s%s" % (" ", " ", " "),  # '    '
    "F": "%s%s%s" % ("\u2514", "\u2500", "\u2500"),  # '└──'
    "W": "%s%s%s" % ("\u250C", "\u2500", "\u2500"),  # '└──'
    "P": "%s%s%s" % ("\u2502", " ", " "),  # '│   '
    "B": "%s%s%s" % ("\u2500", "\u2500", "\u2500"),
}  # '───'

@contextlib.contextmanager
def poolcontext(*args, **kwargs):
    pool = multiprocessing.Pool(*args, **kwargs)
    yield pool
    pool.terminate()

def get_fgv_idxs(kmax):
    # DOES NOT WORK FOR ndim=4 (only for ndim=(5 or 6)) !!! needs to be extended!
    # - get indices of FGVs (hetAB>0 & fixed>0) in bsfs based on kmax
    # - assumes kmax := np.array([hetB, hetA, hetAB, fixed])
    # - used for setting FGVs in data bsfs to 0 
    shape = kmax + 2
    hetAB_idxs, fixed_idx = zip(*[
        (hetAB_idx, fixed_idx) for 
            hetAB_idx, fixed_idx in list(itertools.product(np.arange(shape[-2]), np.arange(shape[-1]))) 
                if hetAB_idx > 0 and fixed_idx >0])
    fgv_idxs = (..., np.array(hetAB_idxs), np.array(fixed_idx))
    return fgv_idxs

def get_parameter_dicts_from_user_parameters(**kwargs):
    # order is determined by input : Ne_A, Ne_B, Ne_A_B, T, me
    # replaces parameters_expanded and only returns the combinations of parameters needed for the GDIs
    # NOT the array of all values by parameter name (ie. parameters_expanded) 
    def value_list_to_array(value_list=[]):
        if len(value_list) == 1:
            value_array = np.array(value_list)
            return value_array
        if len(value_list) == 4:
            value_min, value_max, value_num, distr = value_list
            if distr == 'lin':
                value_array = np.linspace(value_min, value_max, num=value_num, endpoint=True, dtype=np.float64)
            if distr == 'log':
                try:
                    value_array = np.geomspace(value_min, value_max, num=value_num, endpoint=True, dtype=np.float64)
                except ValueError:
                    _value_min = 10 ** (np.log10(value_max) - value_num - (1 if value_min == 0 else 0))
                    value_array = np.geomspace(_value_min, value_max, num=value_num, endpoint=True, dtype=np.float64)
                    if value_min == 0:
                        np.insert(value_array, 0, 0)
            return value_array
        return []
    value_arrays_by_parameter = {k: value_list_to_array(v) for k, v in kwargs.items() if v is not None}
    return [dict(zip(value_arrays_by_parameter.keys(), instance)) for instance in itertools.product(*value_arrays_by_parameter.values())]

class ReportObj(object):
    """Report class for making reports"""

    def __init__(self, width=80):
        self.width = width
        self.out = []

    def add_line(self, prefix="[+]", branch="", left="", center="", right="", fill=" "):
        """
        Writes line to Report Object
        Lines is composed of {prefix} {branch} {left} {fill} {center} {fill} {right}
        SIGNS = {'T': '├──', 'F': '└──', 'S': '    ', 'P': '│   ', 'B': '───'}
        """
        start = (
            "%s %s" % (prefix, "".join([SIGNS[b] for b in branch]))
            if branch
            else prefix
        )
        w = self.width - len(left) - len(right) - len(start)
        self.out.append(
            "%s %s %s %s"
            % (
                start,
                left,
                (" %s " % center).center(w, fill) if center else fill * w,
                right,
            )
        )

    def __add__(self, other):
        if isinstance(other, ReportObj):
            if self.width == other.width:
                result = ReportObj(width=self.width)
                result.out = self.out + other.out
                return result
            else:
                raise ValueError("addition only supported for ReportObj of equal width")
        else:
            raise ValueError("addition only supported for ReportObj's")

    def __radd__(self, other):
        return self.__add__(other)

    def __repr__(self):
        return "\n".join(self.out)

def fullprint(*args, **kwargs):
    from pprint import pprint

    opt = np.get_printoptions()
    np.set_printoptions(threshold=np.inf)
    pprint(*args, **kwargs)
    np.set_printoptions(**opt)

def DOL_to_LOD(DOL):
    """
    converts dict of lists to list of dicts
    """
    reshape_DOL = list(zip(*DOL.values()))
    return [{k: v for k, v in zip(DOL.keys(), sublist)} for sublist in reshape_DOL]


def LOD_to_DOL(LOD):
    """
    converts list of dicts to dict of lists
    """
    reshape_LOD = list(zip(*(d.values() for d in LOD)))
    return {
        k: np.array(v, dtype=np.float64) for k, v in zip(LOD[0].keys(), reshape_LOD)
    }

def _return_np_type(x):
    return np.min_scalar_type((1 if np.min(x) >= 0 else -1) * np.max(np.abs(x)))

class GimbleDemographyInstance(object):
    def __init__(self, model=None, Ne_A=None, Ne_B=None, Ne_A_B=None, me=None, T=None, mu=None, ref_pop=None, sync_pops=None, block_length=None, kmax=None):
        self._SUPPORTED_MODELS = ["DIV", "IM_BA", "IM_AB", "MIG_BA", "MIG_AB"]
        self.model = self.validate_model(model)
        self.order_of_parameters = self.get_order_of_parameters()
        
        self.pop_sizes = set(self.order_of_parameters) - set(['me', 'T'])
        self.mu = mu
        self.block_length = block_length
        self.sync_pops = self.validate_sync_pops(sync_pops)
        self.ref_pop = self.validate_ref_pop(ref_pop)
        self.kmax = self.validate_kmax(kmax)
        #
        self.Ne_A = Ne_A
        self.Ne_B = Ne_B
        self.Ne_A_B = Ne_A_B
        self.me = me
        self.T = T

        self.fixed_parameters = self.get_fixed_parameters()

    def __str__(self):
        return "GimbleDemographyInstance(%s)" % (", ".join(["=".join([k, str(v)]) for k, v in self.get_parameter_dict(nones=False).items()]))

    def __repr__(self):
        return "GimbleDemographyInstance(%s)" % (", ".join(["=".join([k, str(v)]) for k, v in self.get_parameter_dict(nones=True).items()]))

    def get_fixed_parameters(self):
        return set([p for p in self.order_of_parameters if getattr(self, p) is not None])

    def scale_parameters(self, values_by_parameter={}):
        SCALING_FACTOR = 2
        values_by_parameter_scaled = {}
        values_by_parameter_unscaled = {}
        parameters_self = set([parameter for parameter in self.order_of_parameters if getattr(self, parameter, None) is not None])
        parameters_args = set(values_by_parameter.keys())
        parameters_missing = set(self.order_of_parameters) - parameters_self - parameters_args - set(self.sync_pops)
        if not parameters_missing:
            ref_pop = self.ref_pop if not self.ref_pop in self.sync_pops else "Ne_s"
            ref_Ne_value = values_by_parameter[ref_pop] if ref_pop in values_by_parameter else getattr(self, ref_pop)
            for parameter in self.order_of_parameters:
                param = 'Ne_s' if parameter in self.sync_pops else parameter
                value = values_by_parameter[param] if param in values_by_parameter else getattr(self, param)
                if param.startswith("Ne"):
                    values_by_parameter_scaled[parameter] = ref_Ne_value / value
                    values_by_parameter_unscaled[parameter] = value
                elif param.startswith("me"):
                    values_by_parameter_scaled[parameter] = SCALING_FACTOR * ref_Ne_value * value
                    values_by_parameter_unscaled[parameter] = value
                elif param.startswith("T"):
                    values_by_parameter_scaled[parameter] = value / (SCALING_FACTOR * ref_Ne_value) 
                    values_by_parameter_unscaled[parameter] = value
                else:
                    raise ValueError('[X] Unknown parameter %r with value %r' % (parameter, unscaled_value))
            values_by_parameter_scaled['theta_branch'] = SCALING_FACTOR * ref_Ne_value * self.mu * self.block_length # 2Ne*mu
        return (values_by_parameter_scaled, values_by_parameter_unscaled)

    def get_agemo_values(self, scaled_values_by_parameter={}, fallback=False):
        '''
        - fallback means falling back to simpler model if me=0
        '''
        if not scaled_values_by_parameter:
            scaled_values_by_parameter, values_by_parameter_unscaled = self.scale_parameters()
        #print("scaled_values_by_parameter", scaled_values_by_parameter, values_by_parameter_unscaled)
        theta_along_branchtypes = np.full(len(self.kmax), scaled_values_by_parameter['theta_branch'], dtype=np.float64) # len(branch_type_object)
        time = scaled_values_by_parameter.get('T', 0)
        if fallback:
            if 'me' in scaled_values_by_parameter and scaled_values_by_parameter['me'] == 0:
                var_values = [scaled_values_by_parameter[parameter] for parameter in self.order_of_parameters if parameter.startswith("Ne")]
                fallback_flag = True
            else:
                var_values = [scaled_values_by_parameter[parameter] for parameter in self.order_of_parameters if parameter.startswith("Ne") or (parameter.startswith("me"))]
                fallback_flag = False
        else:
            var_values = [scaled_values_by_parameter[parameter] for parameter in self.order_of_parameters if parameter.startswith("Ne") or parameter.startswith("me")]
            fallback_flag = False
        var = np.hstack([np.array(var_values, dtype=np.float64), theta_along_branchtypes])
        return (scaled_values_by_parameter['theta_branch'], var, time, fallback_flag)

    def get_parameter_dict(self, nones=False):
        if nones:
            return {key: value for key, value in self.__dict__.items() if not key.startswith("_")}    
        return {key: value for key, value in self.__dict__.items() if not key.startswith("_") and not value is None}

    def validate_sync_pops(self, sync_pops): 
        if sync_pops is None:
            return []
        pops_invalid = [pop for pop in sync_pops if not pop in self.pop_sizes]
        if pops_invalid:
            return sys.exit("[X] Syncing of population sizes in analysis not possible: %r not valid for model %r (%s)" % (", ".join(pops_invalid), self.model, ", ".join(self.pop_sizes)))
        return sync_pops

    def validate_ref_pop(self, ref_pop):
        if ref_pop is None:
            return None
        if ref_pop in self.pop_sizes:
            return ref_pop
        return sys.exit("[X] Reference population %r not valid for model: %r" % (ref_pop, self.model))

    def validate_kmax(self, kmax):
        '''needs some validation once format is decided upon'''
        return kmax

    def validate_model(self, model):
        if model not in set(self._SUPPORTED_MODELS):
            return sys.exit("[X] Model %s is not supported. Supported models are: %s" % (model, ", ".join(self._SUPPORTED_MODELS)))
        return model

    def get_order_of_parameters(self):
        return {
            "DIV": ['Ne_A_B', 'Ne_A', 'Ne_B', 'T'],
            "MIG_AB": ['Ne_A', 'Ne_B', 'me'],
            "MIG_BA": ['Ne_A', 'Ne_B', 'me'],
            "IM_AB": ['Ne_A_B', 'Ne_A', 'Ne_B', 'me', 'T'],
            "IM_BA": ['Ne_A_B', 'Ne_A', 'Ne_B', 'me', 'T'],
            }.get(self.model, None)

    def get_events(self):
        return {
            "DIV": {'coalescence': ['C_A', 'C_B', 'C_A_B'], 'exodus': [(0, 1, 2)]},
            "MIG_AB": {'coalescence': ['C_A', 'C_B'], 'migration': [(0, 1)]},
            "MIG_BA": {'coalescence': ['C_A', 'C_B'], 'migration': [(0, 1)]},
            "IM_AB": {'coalescence': ['C_A', 'C_B', 'C_A_B'], 'migration': [(0, 1)], 'exodus': [(0, 1, 2)]},
            "IM_BA": {'coalescence': ['C_A', 'C_B', 'C_A_B'], 'migration': [(1, 0)], 'exodus': [(0, 1, 2)]},
            }.get(self.model, None)

    def is_fully_parametrized(self):
        if None in set([getattr(self, parameter, None) for parameter in self.order_of_parameters]):
            return False
        return True

    def get_demes_graph(self):
        '''
        in demes, source and destination of migration is specified FW in time
        in gimble, source and destination of migration is specified BW in time
        '''
        if not self.is_fully_parametrized():
            raise Exception("GimbleDemographyInstance is not fully parametrized: %s" % (self.__repr__)) 
        graph = demes.Builder(time_units="generations")
        if self.model == "DIV":
            graph = demes.Builder(time_units="generations")
            graph.add_deme("A_B", epochs=[dict(end_time=self.T, start_size=self.Ne_A_B, end_size=self.Ne_A_B)])
            graph.add_deme("A", ancestors=["A_B"], defaults=dict(epoch=dict(start_size=self.Ne_A, end_size=self.Ne_A)))
            graph.add_deme("B", ancestors=["A_B"], defaults=dict(epoch=dict(start_size=self.Ne_B, end_size=self.Ne_B)))
        elif self.model == "MIG_AB":
            graph = demes.Builder(time_units="generations")
            graph.add_deme("A", ancestors=None, defaults=dict(epoch=dict(start_size=self.Ne_A, end_size=self.Ne_A)))
            graph.add_deme("B", ancestors=None, defaults=dict(epoch=dict(start_size=self.Ne_B, end_size=self.Ne_B)))
            '''Source and destination demes refer to individuals migrating forwards in time.'''
            # graph.add_migration(source="A", dest="B", rate=self.me) # wrong
            graph.add_migration(source="B", dest="A", rate=self.me)
        elif self.model == "MIG_BA":
            graph = demes.Builder(time_units="generations")
            graph.add_deme("A", ancestors=None, defaults=dict(epoch=dict(start_size=self.Ne_A, end_size=self.Ne_A)))
            graph.add_deme("B", ancestors=None, defaults=dict(epoch=dict(start_size=self.Ne_B, end_size=self.Ne_B)))
            '''Source and destination demes refer to individuals migrating forwards in time.'''
            # graph.add_migration(source="B", dest="A", rate=self.me) # wrong
            graph.add_migration(source="A", dest="B", rate=self.me)
        elif self.model == "IM_AB":
            graph = demes.Builder(time_units="generations")
            graph.add_deme("A_B", epochs=[dict(end_time=self.T, start_size=self.Ne_A_B, end_size=self.Ne_A_B)])
            graph.add_deme("A", ancestors=["A_B"], defaults=dict(epoch=dict(start_size=self.Ne_A, end_size=self.Ne_A)))
            graph.add_deme("B", ancestors=["A_B"], defaults=dict(epoch=dict(start_size=self.Ne_B, end_size=self.Ne_B)))
            '''Source and destination demes refer to individuals migrating forwards in time.'''
            # graph.add_migration(source="A", dest="B", rate=self.me) # wrong
            graph.add_migration(source="B", dest="A", rate=self.me)
        elif self.model == "IM_BA":
            graph = demes.Builder(time_units="generations")
            graph.add_deme("A_B", epochs=[dict(end_time=self.T, start_size=self.Ne_A_B, end_size=self.Ne_A_B)])
            graph.add_deme("A", ancestors=["A_B"], defaults=dict(epoch=dict(start_size=self.Ne_A, end_size=self.Ne_A)))
            graph.add_deme("B", ancestors=["A_B"], defaults=dict(epoch=dict(start_size=self.Ne_B, end_size=self.Ne_B)))
            '''Source and destination demes refer to individuals migrating forwards in time.'''
            # graph.add_migration(source="B", dest="A", rate=self.me) # wrong
            graph.add_migration(source="A", dest="B", rate=self.me)
        else:
            pass
        return graph

    def add_optimize_values(self, value_by_nlopt_parameter):
        for parameter in self.order_of_parameters:
            value = value_by_nlopt_parameter.get(parameter, value_by_nlopt_parameter.get('Ne_s', None))
            setattr(self, parameter, value)

    def demes_dump(self, fn=None, fmt='yaml', simplified=True):
        if self.is_fully_parametrized():
            graph = self.get_demes_graph()
            if fn is None:
                return demes.dumps(graph.resolve(), format=fmt, simplified=simplified)
            demes.dump(graph.resolve(), filename=fn, format=fmt, simplified=simplified)
        return ""

    def get_demography(self):
        demes_graph = self.get_demes_graph()
        return msprime.Demography.from_demes(demes_graph.resolve())

def recursive_get_size(path):
    """Gets size in bytes of the given path, recursing into directories."""
    if os.path.isfile(path):
        return os.path.getsize(path)
    if not os.path.isdir(path):
        return 0
    return sum(
        recursive_get_size(os.path.join(path, name)) for name in os.listdir(path)
    )

def parse_csv(csv_f="", dtype=[], usecols=[], sep=",", header=None):
    '''dtypes := "object", "int64", "float64", "bool", "datetime64", "timedelta", "category"'''
    df = pd.read_csv(
        csv_f,
        sep=sep,
        usecols=usecols,
        names=list(dtype.keys()),
        header=header,
        dtype=dtype,
    )
    if df.isnull().values.any():
        sys.exit("[X] Bad file format %r." % csv_f)
    return df

def format_query_meta(meta, ignore_long=False):
    LONG_THRESHOLD = 80
    lines = []
    for key, value in meta.items():
        if isinstance(value, int):
            formatted_value = format_count(value)
        elif isinstance(value, list):
            formatted_value = str(value)
        else:
            formatted_value = value
        if ignore_long:
            if len(formatted_value) < LONG_THRESHOLD:
                lines.append("[+]\t%s: %s" % (key, formatted_value))
        else:
            lines.append("[+]\t%s: %s" % (key, formatted_value))
    return "\n".join(lines)

def format_bases(bases):
    if bases in set(["-", "N/A"]):
        return bases
    return "%s b" % format(bases, ",d")

def format_percentage(fraction, precision=2):
    if fraction in set(["-", "N/A"]):
        return fraction
    return "{:.{}%}".format(fraction, precision)

def format_proportion(fraction, precision=2):
    if fraction in set(["-", "N/A"]):
        return fraction
    return "{:.{}f}".format(fraction, precision)

def format_count(count):
    if count in set(["-", "N/A"]):
        return count
    return "%s" % str(format(count, ",d"))

def format_time(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:0>2}h:{:0>2}m:{:06.3f}s".format(int(hours), int(minutes), seconds)

def format_bytes(size, precision=1):
    power = 2**10
    n = 0
    power_labels = {0: " B", 1: "KB", 2: "MB", 3: "GB", 4: "TB"}
    while size > power:
        size /= power
        n += 1
    return "{:.{}f} {}".format(size, precision, power_labels[n])

def get_n50_from_lengths(lengths):
    length_sorted = sorted(lengths, reverse=True)
    cum_sum = np.cumsum(length_sorted)
    half = int(sum(lengths) / 2)
    cum_sum_2 = min(cum_sum[cum_sum >= half])
    n50_idx = np.where(cum_sum == cum_sum_2)
    return length_sorted[int(n50_idx[0][0])]

def parse_intervals(bed_f, target_sequences, target_samples):
    intervals_df = parse_csv(
        csv_f=bed_f,
        sep="\t",
        usecols=[0, 1, 2, 4],
        dtype={
            "sequence": "category",
            "start": "int64",
            "end": "int64",
            "samples": "category",
        },
        header=None,
    )
    # Subset target_sequences, sort, reindex
    intervals_df = (
        intervals_df[intervals_df["sequence"].isin(target_sequences)]
        .sort_values(["sequence", "start"], ascending=[True, True])
        .reset_index(drop=True)
    )
    # Convert samples column to sample-matrix
    intervals_df = pd.concat(
        [
            intervals_df,
            intervals_df.samples.str.get_dummies(sep=",").filter(target_samples),
        ],
        axis=1,
    ).drop(columns=["samples"])
    # Subset target samples
    samples_in_df = [sample for sample in intervals_df.columns[3:]]
    target_samples_in_df = ordered_intersect(
        a=samples_in_df, b=target_samples, order="a"
    )
    target_samples_not_in_df = set(target_samples).difference(set(target_samples_in_df))
    if target_samples_not_in_df:
        sys.exit(
            "[X] Samples in SAMPLE_FILE not found in BED_FILE: %s"
            % ", ".join(list(target_samples_not_in_df))
        )
    non_target_samples_in_df = set(target_samples).difference(set(samples_in_df))
    intervals_df = intervals_df.drop(non_target_samples_in_df, axis=1)
    # Add length column
    intervals_df["length"] = intervals_df["end"] - intervals_df["start"]
    return intervals_df

def bsfs_to_2d(bsfs):
    """Converts 4D bsfs to 2D array with counts, mutuples.
       Converts 5D bsfs to 2D array with window_idx, counts, mutuples.

    Parameters
    ----------
    bsfs : ndarray, int, ndim (4 or 5)

    Returns
    -------
    out : ndarray, int, ndim (2)
    """
    if not np.any(bsfs):
        return None
    non_zero_idxs = np.nonzero(bsfs)
    if bsfs.ndim == 4:  # blocks
        return np.concatenate(
            [
                bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1),
                np.array(non_zero_idxs).T,
            ],
            axis=1,
        )
    elif bsfs.ndim == 5:  # windows
        non_zero_idxs_array = np.array(non_zero_idxs).T
        first = non_zero_idxs_array[:, 0].reshape(non_zero_idxs[0].shape[0], 1)
        second = bsfs[non_zero_idxs].reshape(non_zero_idxs[0].shape[0], 1)
        third = non_zero_idxs_array[:, 1:]
        return np.concatenate([first, second, third], axis=1)
    else:
        raise ValueError("bsfs_to_2d: bsfs.ndim must be 4 (blocks) or 5 (windows)")

def calculate_blocks_report_metrics(
    tally, sample_sets_count, block_length, intervals_span
):
    BRM = collections.defaultdict(lambda: "-")
    BRM["blocks_total"] = np.sum(tally[:, 0]) if np.any(tally) else 0
    effective_length = block_length * BRM["blocks_total"]
    BRM["interval_coverage"] = (
        effective_length / sample_sets_count
    ) / intervals_span  # / sample_sets_count
    if BRM["blocks_total"]:
        BRM["blocks_invariant"] = tally[0, 0] / BRM["blocks_total"]
        BRM["blocks_fgv"] = (
            np.sum(tally[(tally[:, 3] > 0) & (tally[:, 4] > 0)][:, 0])
            / BRM["blocks_total"]
        )
        hetB = np.sum((tally[:, 0] * tally[:, 1]))
        hetA = np.sum((tally[:, 0] * tally[:, 2]))
        hetAB = np.sum((tally[:, 0] * tally[:, 3]))
        fixed = np.sum((tally[:, 0] * tally[:, 4]))
        # not sure about heterozygosity ...
        BRM["heterozygosity_A"] = (hetA + hetAB) / effective_length
        BRM["heterozygosity_B"] = (hetB + hetAB) / effective_length  # hetB or hetA ?
        BRM["heterozygosity_intra"] = (
            float(fractions.Fraction(1, 2) * (hetA + hetB) + hetAB)
        ) / effective_length
        # BRM['heterozygosity_B'] = BRM['heterozygosity_A'] # this is how it was before...
        BRM["dxy"] = ((hetA + hetB + hetAB) / 2.0 + fixed) / effective_length
        mean_pi = (BRM["heterozygosity_A"] + BRM["heterozygosity_B"]) / 2.0
        BRM["fst"] = (
            (BRM["dxy"] - mean_pi) / (BRM["dxy"] + mean_pi) if (BRM["dxy"] + mean_pi) else "N/A"
        )
        total_segregating = np.sum(tally[:, 0, None] * tally[:, 1:])
        BRM["pi"] = (
            float(
                fractions.Fraction(1, 2) * (hetA + hetB)
                + fractions.Fraction(2, 3) * (hetAB + fixed)
            )
            / effective_length
        )
        BRM["watterson_theta"] = (
            total_segregating / float(harmonic(3)) / effective_length
        )
    return BRM

def write_info_report(version, report, prefix):
    # OUTPUTLIB
    txt = ["# %s" % version, str(report)]
    out_f = "%s.info.txt" % prefix
    print("[+] Writing info file %s ..." % out_f)
    with open(out_f, "w") as out_fh:
        out_fh.write("\n".join(txt) + "\n")
    return out_f

def get_popgen_metrics(array, sites=0):
    """only works for mutypes=4
    Possible inputs:
    A) block-tally  := ndim (2); shape (n, 5)
    B) window-tally := ndim (2); shape (n, 6)
    C) window-bsfs  := ndim (4); shape (maxk_m1, maxk_m2, maxk_m3, maxk_m4)
    D) window-bsfs  := ndim (5); shape (w, maxk_m1, maxk_m2, maxk_m3, maxk_m4)
    sites = (block_length * block_count)
    """
    assert sites > 0, "sites must be positive integer"
    array = array if array.ndim == 2 else bsfs_to_2d(np.array(array))
    if array.shape[1] == 5:  # block tally
        mutype_array = np.array(
            np.sum(array[:, 0, np.newaxis] * array[:, 1:], axis=0)
        ).reshape(
            -1, 4
        )  # must be reshaped to 2D, so that indexing works
    elif array.shape[1] == 6:  # window tally
        mutype_array = np.vstack(
            [
                np.bincount(array[:, 0], weights=(array[:, 1] * array[:, (2 + m_idx)]))
                for m_idx in range(4)
            ]
        ).T
    else:
        raise ValueError(
            "get_popgen_metrics: if tally.ndim is 2, then shape must be (n, 5) (blocks) or (n, 6) (windows)"
        )
    heterozygosity_A = (mutype_array[:, 1] + mutype_array[:, 2]) / sites
    heterozygosity_B = (mutype_array[:, 0] + mutype_array[:, 2]) / sites
    d_xy = (
        (mutype_array[:, 1] + mutype_array[:, 0] + mutype_array[:, 2]) / 2.0
        + mutype_array[:, 3]
    ) / sites
    mean_pi = (heterozygosity_A + heterozygosity_B) / 2.0
    f_st = np.full(mutype_array.shape[0], np.nan)
    np.true_divide(
        (d_xy - mean_pi), (d_xy + mean_pi), out=f_st, where=((d_xy + mean_pi) > 0)
    )
    return np.vstack([heterozygosity_A, heterozygosity_B, d_xy, f_st])

def check_unique_pos(pos_array):
    unique_pos, counts_pos = np.unique(pos_array, return_counts=True)
    duplicates = unique_pos[counts_pos > 1]
    if duplicates.any():
        print(
            "\n[-] %s VCF records with non-unique positions found. Rescuing records by shifting position... (abort if this is not desired)"
            % (len(duplicates))
        )
        pos_array = fix_pos_array(pos_array)
    return pos_array

def fix_pos_array(pos_array):
    """
    De-duplicates array by shifting values forward until there aren't any collisions
    """
    # get boolean array for first and subsequent duplicates (True) (required sorted)
    idxs = np.insert((np.diff(pos_array) == 0).astype(bool), 0, False)
    if np.any(idxs):
        # if there are duplicates, get new values by incrementing by one
        new_values = pos_array[idxs] + 1
        # get non-duplicate values
        uniq_values = pos_array[~idxs]
        # insert new_values in non-duplicated values (required sorted)
        new_idxs = np.searchsorted(uniq_values, new_values)
        # recursive call
        return fix_pos_array(np.sort(np.insert(uniq_values, new_idxs, new_values)))
    # if there are no duplicated values
    return pos_array

def szudzik_pairing(folded_minor_allele_counts):
    # adapted from: https://drhagen.com/blog/superior-pairing-function/
    return np.where(
        (folded_minor_allele_counts[:, 0] >= folded_minor_allele_counts[:, 1]),
        (
            np.square(folded_minor_allele_counts[:, 0])
            + folded_minor_allele_counts[:, 0]
            + folded_minor_allele_counts[:, 1]
        ),
        (
            np.square(folded_minor_allele_counts[:, 1])
            + folded_minor_allele_counts[:, 0]
        ),
    )

def _harmonic(a, b):
    if b - a == 1:
        return fractions.Fraction(1, a)
    m = (a + b) // 2
    return _harmonic(a, m) + _harmonic(m, b)

def harmonic(n):
    """https://fredrik-j.blogspot.com/2009/02/how-not-to-compute-harmonic-numbers.html"""
    return _harmonic(1, n + 1)


def chisq(sample_set_idxs, window_samples_set_idxs):
    if window_samples_set_idxs.size == 0:
        return 0.0
    spacer = np.max(window_samples_set_idxs) + 1
    window_count = window_samples_set_idxs.shape[0]
    window_size = window_samples_set_idxs.shape[1]
    temp_sites = window_samples_set_idxs + (
        spacer * np.arange(window_count, dtype=np.int64).reshape(window_count, 1)
    )
    obs = np.bincount(temp_sites.ravel(), minlength=(window_count * spacer)).reshape(
        -1, spacer
    )[:, sample_set_idxs]
    # print(obs)
    exp = np.full(obs.shape, window_size / sample_set_idxs.shape[0])
    return np.sum((((obs - exp) ** 2) / exp), axis=1)  # / sample_sets.shape[0]


def mse(sample_set_idxs, window_samples_set_idxs):
    """measure of eveness"""
    if window_samples_set_idxs.size == 0:
        return 0.0
    spacer = np.max(window_samples_set_idxs) + 1
    window_count = window_samples_set_idxs.shape[0]
    window_size = window_samples_set_idxs.shape[1]
    temp_sites = window_samples_set_idxs + (
        spacer * np.arange(window_count, dtype=np.int64).reshape(window_count, 1)
    )
    obs = np.bincount(temp_sites.ravel(), minlength=(window_count * spacer)).reshape(
        -1, spacer
    )[:, sample_set_idxs]
    exp = np.full(obs.shape, window_size / sample_set_idxs.shape[0])
    # Gertjan: scale by max
    max_mse_obs = np.zeros(sample_set_idxs.shape[0])
    max_mse_obs[0] = window_size
    max_mse = np.sum(((max_mse_obs - exp) ** 2), axis=1)
    return np.sum((((obs - exp) ** 2)), axis=1) / max_mse


def blocks_to_windows(
    sample_set_idxs,
    block_variation,
    start_array,
    end_array,
    block_sample_set_idxs,
    window_size,
    window_step,
):
    # order of blocks is defined by end_array
    # coordinate_sorted_idx is the order of blocks if one were to sort them by end_array
    coordinate_sorted_idx = np.argsort(end_array)
    # elements in windows are defined by window_idxs -> shape(n, window_size)
    window_idxs = np.arange(coordinate_sorted_idx.shape[0] - window_size + 1)[
        ::window_step, None
    ] + np.arange(window_size)
    # window_idxs_alt = np.arange((block_variation.shape[0] - window_size) + 1)[::window_step, None] + np.arange(window_size)
    # all taking is done with coordinate_sorted_idx and window_idxs
    window_variation = block_variation.take(coordinate_sorted_idx, axis=0).take(
        window_idxs, axis=0
    )
    block_starts = start_array.take(coordinate_sorted_idx, axis=0).take(
        window_idxs, axis=0
    )
    window_starts = np.min(
        start_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0),
        axis=1,
    ).T
    block_ends = end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0)
    window_ends = np.max(
        end_array.take(coordinate_sorted_idx, axis=0).take(window_idxs, axis=0), axis=1
    ).T
    # needs some solution for chisq-calculation by window ...
    window_samples_set_idxs = block_sample_set_idxs.take(
        coordinate_sorted_idx, axis=0
    ).take(window_idxs, axis=0)
    # np.set_printoptions(threshold=sys.maxsize)
    balance = chisq(sample_set_idxs, window_samples_set_idxs)
    mse_sample_set_cov = mse(sample_set_idxs, window_samples_set_idxs)
    block_midpoints = (block_starts / 2) + (block_ends / 2)
    window_pos_mean = np.rint(np.mean(block_midpoints, axis=1).T)
    window_pos_median = np.rint(np.median(block_midpoints, axis=1).T)
    return (
        window_variation,
        window_starts,
        window_ends,
        window_pos_mean,
        window_pos_median,
        balance,
        mse_sample_set_cov,
    )


def block_sites_to_variation_arrays(
    block_sites, cols=np.array([1, 2, 3]), max_type_count=7
):
    temp_sites = block_sites + (
        max_type_count
        * np.arange(block_sites.shape[0], dtype=np.int64).reshape(
            block_sites.shape[0], 1
        )
    )
    multiallelic, missing, monomorphic, variation = np.hsplit(
        np.bincount(
            temp_sites.ravel(), minlength=(block_sites.shape[0] * max_type_count)
        ).reshape(-1, max_type_count),
        cols,
    )
    return (multiallelic, missing, monomorphic, variation)


def gt2fmac(genotype_array):
    """turns scikit-allel/numpy-array (n,s,p) of genotypes into folded-minor-allel-count-np-array (n,s)
    [To Do] multiallelics have to be done differently if we use non-(n,2,2)-GTs ..."""
    np_genotype_array = np.array(genotype_array)
    sa_genotype_array = allel.GenotypeArray(genotype_array)
    n_samples = sa_genotype_array.n_samples
    n_variants = sa_genotype_array.n_variants
    np_allele_count_array = np.ma.masked_equal(
        sa_genotype_array.count_alleles(), 0, copy=True
    )
    missing_genotypes = sa_genotype_array.is_missing()
    if np.all(missing_genotypes):
        folded_minor_allele_counts = (
            np.ones((n_variants, n_samples), dtype=np.int8) * -1
        )  # -1, -1 for missing => -1
    else:
        idx_max_global_allele_count = np.nanargmax(np_allele_count_array, axis=1)
        idx_min_global_allele_count = np.nanargmin(np_allele_count_array, axis=1)
        has_major_allele = idx_max_global_allele_count != idx_min_global_allele_count
        idx_min_prime_allele = np.amin(np_genotype_array[:, 0], axis=1)
        idx_min_global_allele = np.amin(np.amin(np_genotype_array, axis=1), axis=1)
        idx_max_global_allele = np.amax(np.amax(np_genotype_array, axis=1), axis=1)
        idx_major_allele = np.where(
            has_major_allele, idx_max_global_allele_count, idx_min_prime_allele
        )
        idx_minor_allele = np.where(
            has_major_allele,
            idx_min_global_allele_count,
            np.where(
                (idx_min_global_allele == idx_min_prime_allele),
                np.max((idx_min_global_allele, idx_max_global_allele), axis=0),
                np.min((idx_min_global_allele, idx_max_global_allele), axis=0),
            ),
        )
        allele_map = np.ones((np_allele_count_array.shape), dtype=np.int64) * np.arange(
            np_allele_count_array.shape[-1], dtype=np.int8
        )  # changed to np.int8
        # for each genotype (np.arange(allele_map.shape[0])), set minor allele to 1 (1st do minor, so that overwritten if monomorphic)
        allele_map[np.arange(allele_map.shape[0]), idx_minor_allele] = 1
        # for each genotype (np.arange(allele_map.shape[0])), set major allele to 0
        allele_map[np.arange(allele_map.shape[0]), idx_major_allele] = 0
        folded_minor_allele_counts = sa_genotype_array.map_alleles(allele_map).to_n_alt(
            fill=-1
        )
        folded_minor_allele_counts[np.any(missing_genotypes, axis=1)] = (
            np.ones(n_samples) * -1
        )  # -1, -1 for missing => -1
        # multiallelics take precedence over missing
        folded_minor_allele_counts[(np_allele_count_array.count(axis=1) > 2)] = np.ones(
            2
        ) * (
            -1,
            -2,
        )  # -1, -2 for multiallelic => -2
    return folded_minor_allele_counts


def intervals_to_sites(intervals):
    """starts = np.array([0, 5, 8, 11, 15])
       ends   = np.array([2, 8, 9, 13, 18])
    clens : array([2, 5, 6, 8, 11])
    _sites: array([1, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1])
    _sites: array([0, 1, 1, 1, 1, 1,  1, 1, 1, 1, 1])
    _sites: array([1, 1, 4, 1, 1, 1,  3, 1, 3, 1, 1])
    sites : array([0, 1, 5, 6, 7, 8, 11,12,15,16,17])"""
    starts, ends = intervals
    if starts is None or ends is None:
        return None
    clens = np.cumsum(ends - starts)
    if np.any(clens):
        sites = np.ones(clens[-1], dtype=np.int64)
        sites[0] = starts[0]
        sites[clens[:-1]] = starts[1:] - ends[:-1] + 1
        sites = sites.cumsum()
        return sites
    return None


def sites_to_blocks(sites, block_length, block_span, sample_set, debug=False):
    # block_sites are 0-based, but they are SITES (numbering the bases) as opposed to coordinates (numbering between the bases)
    # 0 1 2 3 4 5 6 BED
    # A T G C A G
    # 1 2 3 4 5 6  VCF
    # 0 1 2 3 4 5  Scikit allel, sites
    if sites is None:
        return None
    max_gap = block_span - block_length
    if max_gap == 0:  # no gaps within blocks
        block_sites = sites[: block_length * (sites.shape[0] // block_length)].reshape(
            -1, block_length
        )
    else:  # gaps within blocks are allowed
        block_sites = np.concatenate(
            [
                x[: block_length * (x.shape[0] // block_length)].reshape(
                    -1, block_length
                )
                for x in np.split(sites, np.where(np.diff(sites) > max_gap)[0] + 1)
            ]
        )
    ## no splitting
    # block_sites = sites[:block_length * (sites.shape[0] // block_length)].reshape(-1, block_length)
    ## yes splitting
    # block_sites = np.concatenate([
    #    x[:block_length * (x.shape[0] // block_length)].reshape(-1, block_length)
    #        for x in np.split(sites, np.where(np.diff(sites) > max_gap)[0] + 1)])
    block_sites_valid_mask = (
        (block_sites[:, -1] - block_sites[:, 0] + 1)
    ) <= block_span  # +1 is needed because block sites are sites (not coordinates)
    if debug:
        print("[+] sample_set", sample_set)
        print("[+] sites", sites)
        success = (
            (block_sites[block_sites_valid_mask].shape[0] * block_length)
            / sites.shape[0]
            if np.count_nonzero(block_sites[block_sites_valid_mask])
            else 0
        )
        print("[+] blocks", block_sites[block_sites_valid_mask])
        print("[+] block_sites_valid_mask", block_sites_valid_mask)
        print("[+] block_sites_valid", (block_sites[:, -1] - block_sites[:, 0] + 1))
        print(
            "[+] sites=%s : blocks=%s success=%.2f"
            % (sites.shape[0], block_sites[block_sites_valid_mask].shape[0], success)
        )
    if np.any(block_sites_valid_mask):
        return block_sites[block_sites_valid_mask]
    return None


def subset_gt_matrix(meta_seqs, sample_set, indices, gt_matrix):
    if gt_matrix is not None:
        sample_set_vcf_idxs = np.array(
            [meta_seqs["variants_idx_by_sample"][sample] for sample in sample_set]
        )
        return gt_matrix.subset(indices, sample_set_vcf_idxs)
    return None


def blocks_to_arrays(blocks, gts, pos):
    starts = np.array(blocks[:, 0], dtype=np.int64)
    ends = np.array(blocks[:, -1] + 1, dtype=np.int64)  # BED: end+1
    pos_in_block_sites = np.isin(pos, blocks, assume_unique=True)
    if np.any(pos_in_block_sites):  # if variants in blocks
        folded_minor_allele_counts = gt2fmac(gts)
        block_sites_in_pos = np.isin(blocks, pos, assume_unique=True)
        blocks[block_sites_in_pos] = (
            szudzik_pairing(folded_minor_allele_counts) + 2
        )  # add 2 so that not negative for bincount
        blocks[
            ~block_sites_in_pos
        ] = 2  # monomorphic = 2 (0 = multiallelic, 1 = missing)
        temp_sites = blocks + (
            7 * np.arange(blocks.shape[0], dtype=np.int64).reshape(blocks.shape[0], 1)
        )
        multiallelic, missing, monomorphic, variation = np.hsplit(
            np.bincount(temp_sites.ravel(), minlength=(blocks.shape[0] * 7)).reshape(
                -1, 7
            ),
            np.array([1, 2, 3]),
        )
    else:
        multiallelic = np.zeros((blocks.shape[0], 1), dtype=np.int64)
        missing = np.zeros((blocks.shape[0], 1), dtype=np.int64)
        monomorphic = np.full((blocks.shape[0], 1), blocks.shape[1], dtype=np.int64)
        variation = np.zeros((blocks.shape[0], 4), dtype=np.int64)
    return (starts, ends, multiallelic, missing, monomorphic, variation)

def tally_to_bsfs(tally, max_k, data="blocks"):
    if data == "blocks":
        counts = tally[:, 0]
        dtype = _return_np_type(counts)
        out = np.zeros((np.array(max_k, dtype="int") + 2), dtype)
        out[tuple(tally[:, 1:].T)] = counts
    elif data == "windows":
        counts = tally[:, 1]
        window_idxs = tally[:, 0][np.newaxis]
        num_windows = np.max(window_idxs) + 1
        out_shape = np.insert((np.array(max_k, dtype="int") + 2), 0, num_windows)
        dtype = _return_np_type(counts)
        out = np.zeros(out_shape, dtype)
        idx_and_counts = np.hstack((window_idxs.T, tally[:, 2:]))
        out[tuple(idx_and_counts.T)] = counts

    else:
        raise ValueError(f"2d_to_bsfs not implemtened for {data}.")
    return out


def tally_variation(variation, form="bsfs", max_k=None):
    """
    Parameters
    ----------
    variation : ndarray, np.uint64, ndim (2 or 3)

    form : 'bsfs' for bsfs tensors
           'tally' for 2D tally array

    max_k : ndarray, int, ndim (1) for capping mutypes
            defaults to np.array([8,8,8,8]) if (max_k == None & form == 'bsfs')
    Returns
    -------
    out : ndarray, int
          ndim (4 or 5) if form == 'bsfs'
          ndim (2 or 3) if form == 'tally'
    """
    if max_k is None:
        max_k = np.array([8, 8, 8, 8]) if form == "bsfs" else None
    else:
        max_k = max_k + 1  # for clipping
    if variation.ndim == 2:
        mutuples = np.clip(variation, 0, max_k)
    elif variation.ndim == 3:
        index = np.repeat(np.arange(variation.shape[0]), variation.shape[1]).reshape(
            -1, 1
        )
        mutuples = np.concatenate(
            (index, np.clip(variation.reshape(-1, variation.shape[-1]), 0, max_k)),
            axis=-1,
        )
    else:
        raise ValueError(
            "variation.ndim is %r, should either be 2 (blocks) or 3 (windows)"
            % variation.ndim
        )
    try:
        #print("mutuples", mutuples)
        mutuples_unique, counts = np.unique(mutuples, return_counts=True, axis=0)
        #print("mutuples_unique", mutuples_unique.shape)
        #print("counts", counts.shape)
        #print("form", form)
        #print("variation", variation.shape, variation)
        dtype = _return_np_type(counts)
        if form == "bsfs":
            out = np.zeros(tuple(max_k + 1), dtype) if variation.ndim == 2 else np.zeros(tuple([variation.shape[0]]) + tuple(max_k + 1), dtype)
            #print("out", out.shape)
            out[tuple(mutuples_unique.T)] = counts
        elif form == "tally":
            out = np.concatenate(
                (counts.reshape(counts.shape[0], 1), mutuples_unique), axis=1
            )
            if variation.ndim == 3:
                out[:, [0, 1]] = out[
                    :, [1, 0]
                ]  # for window variation, order must be [idx, count, mutuple]
            out = out.astype(dtype)
        else:
            raise ValueError("form must be %r or %r, was %r" % ("bsfs", "tally", form))
    except MemoryError as e:
        sys.exit(
            "[X] tally_variation() ran out of memory. Try specifying lower k-max values. %s."
            % str(e)
        )
    return out

def calculate_marginality_of_variation(data, max_k=None):
    """to be run on variation arrays of ndim = 2 or 3"""
    assert data.shape[-1] == 4 and (
        data.ndim == 3 or data.ndim == 2
    ), "[X] data.ndim must be 2 or 3, data.shape[-1] must be 4"
    if max_k is None:
        return 0.0
    is_above_max_k = np.any((data - max_k) > 0, axis=-1)
    return np.sum(is_above_max_k) / is_above_max_k.flatten().shape[0]

def ordered_intersect(a=[], b=[], order="a"):
    A, B = a, b
    if order == "b":
        B, A = a, b
    return [_a for _a in A if _a in set(B)]

def gridsearch_dask(tally=None, grid=None, num_cores=1, chunksize=500, desc=None):
    """returns 2d array of likelihoods of shape (windows, grid)"""
    if grid is None or tally is None:
        return None
    grid_log = np.zeros(grid.shape)
    np.log(grid, where=grid > 0, out=grid_log)
    tally_expansion_axis = (0, 1) if tally.ndim == 4 else 1
    tally = np.expand_dims(tally, axis=tally_expansion_axis)
    tally_chunks = tuple((min(x, chunksize) if i==0 else x) for i, x in enumerate(tally.shape))
    grid_chunks = tuple((min(x, chunksize) if i==0 else x) for i, x in enumerate(grid_log.shape))
    scheduler = "processes" if num_cores > 1 else "single-threaded"
    with warnings.catch_warnings():
        with dask.config.set(scheduler=scheduler, n_workers=num_cores):
            warnings.simplefilter("ignore")
            tally = dask.array.from_array(tally, chunks=tally_chunks)
            grid_log = dask.array.from_array(grid_log, chunks=grid_chunks)
            product = dask.array.multiply(tally, grid_log)
            # product:
            # - keep first two dimensions (product.shape[0], product.shape[1])
            # - reshape to the prod of remaining dimensions: np.prod((4,4,4,4)) = 256 
            product = product.reshape((product.shape[0], product.shape[1], np.prod(product.shape[2:])))
            #print("# product.ndim", product.ndim, product.shape, product.dtype, product)
            result = dask.array.sum(product, axis=-1)
            with TqdmCallback(desc=desc, ncols=100, position=0):
                out = result.compute()
            #print(np.sum(out))
            return out
          
class Store(object):
    def __init__(self, prefix=None, path=None, create=False, overwrite=False):
        self.prefix = (
            prefix if not prefix is None else str(pathlib.Path(path).resolve().stem)
        )
        self.path = path if not path is None else "%s.z" % prefix
        self.data = self._init_store(create, overwrite)

    def tree(self):
        print(self.data.tree())

    def log_action(self, module, command):
        self.data.attrs[module] = command

    def get_stage(self, stage):
        return self.data.attrs[stage]

    def has_stage(self, stage):
        return stage in self.data.attrs

    def list_keys(self, category=None):
        category = category if not category in set(['all', None]) else None
        available_keys_by_category = collections.defaultdict(set)
        max_depth_by_key = {
            "blocks": 1,
            "windows": 1,
            "tally": 2,
            "makegrid": 2,
            "optimize": 3,
            "gridsearch": 3,
            "simulate": 2,
        }
        category_by_module = {
            "blocks": "parse",
            "windows": "parse",
            "tally": "tally",
            "makegrid": "makegrid",
            "simulate": "simulate",
            "optimize": "optimize",
            "gridsearch": "gridsearch",
        }
        def data_key_finder(path):
            key_list = str(path).split("/")
            module = key_list[0]
            if len(key_list) == max_depth_by_key.get(module, None):
                available_keys_by_category[category_by_module[module]].add(
                    "/".join(key_list[0 : max_depth_by_key.get(key_list[0], 0)])
                )
        self.data.visit(data_key_finder)
        if available_keys_by_category:
            if category and category in available_keys_by_category:
                print("[#] %s" % category)
                print("\t- %s" % "\n\t- ".join(sorted(available_keys_by_category[category])))
            else:
                for category, available_keys in available_keys_by_category.items():
                    if category == "tally":
                        available_keys = [key for key in available_keys if key not in set(['tally/blocks_raw', 'tally/windows_raw'])]
                    print("[#] %s" % category)
                    print("\t- %s" % "\n\t- ".join(sorted(available_keys)))
        else:
            print("[X] Gimble store %r seems to be empty." % self.path)
        sys.exit(1)

    def validate_key(self, key, category=None):
        meta = self._get_meta(key)
        if meta:
            return key
        print("[X] Label %r not found in store. Available labels:" % str(key))
        self.list_keys(key)

    def measure(self, genome_f=None, sample_f=None, bed_f=None, vcf_f=None):
        measure_key = "seqs/"
        self._set_meta(measure_key)
        print("[#] Processing GENOME_FILE %r." % genome_f)
        self._read_sequences(measure_key, genome_f)
        print("[#] Processing SAMPLE_FILE %r." % sample_f)
        self._read_samples(measure_key, sample_f)
        print("[#] Processing BED_FILE %r." % bed_f)
        self._read_intervals(measure_key, bed_f)
        print("[#] Processing VCF_FILE %r." % vcf_f)
        self._read_variants(measure_key, vcf_f)
        # print(self.data.tree())

    def _read_sequences(self, measure_key, genome_f):
        sequences_df = parse_csv(
            csv_f=genome_f,
            sep="\t",
            usecols=[0, 1],
            dtype={"sequence_id": "category", "sequence_length": "int64"},
            header=None,
        )
        # meta = self._get_meta(measure_key)
        meta = self._get_meta("seqs")
        meta["seq_names"] = sequences_df["sequence_id"].to_list()
        meta["seq_lengths"] = sequences_df["sequence_length"].to_list()
        meta["seq_n50"] = get_n50_from_lengths(meta["seq_lengths"])
        meta["genome_f"] = genome_f

    def _read_samples(self, measure_key, sample_f):
        try:
            samples_df = parse_csv(
                csv_f=sample_f,
                sep=",",
                usecols=[0, 1],
                dtype={"samples": "object", "populations": "category"},
                header=None,
            )
        except pd.errors.ParserError:
            sys.exit("[X] --sample_f %r can't be parsed. Make sure it is a CSV file composed of two columns (sample ID and population ID)." % sample_f)
        # meta = self._get_meta(measure_key)
        meta = self._get_meta("seqs")
        meta["samples"] = [
            sample.rstrip() for sample in samples_df["samples"].to_list()
        ]
        meta["populations"] = [
            population.rstrip() for population in samples_df["populations"].to_list()
        ]
        meta["population_ids"] = sorted(set(meta["populations"]))
        meta["population_by_letter"] = {
            letter: population_id
            for population_id, letter in zip(
                meta["population_ids"], string.ascii_uppercase
            )
        }
        meta["population_by_sample"] = {
            sample: population
            for sample, population in zip(meta["samples"], meta["populations"])
        }
        meta["sample_sets"] = [
            tuple(
                sorted(
                    x,
                    key=(
                        meta["population_by_sample"].get
                        if meta["population_by_sample"][x[0]]
                        != meta["population_by_sample"][x[1]]
                        else None
                    ),
                )
            )
            for x in itertools.combinations(meta["population_by_sample"].keys(), 2)
        ]
        meta["sample_sets_inter"] = [
            False
            if len(set([meta["population_by_sample"][sample] for sample in sample_set]))
            == 1
            else True
            for sample_set in meta["sample_sets"]
        ]
        meta["sample_sets_intra_A"] = [
            all(
                [
                    meta["population_by_sample"][sample] == meta["population_ids"][0]
                    for sample in sample_set
                ]
            )
            for sample_set in meta["sample_sets"]
        ]
        meta["sample_sets_intra_B"] = [
            all(
                [
                    meta["population_by_sample"][sample] == meta["population_ids"][1]
                    for sample in sample_set
                ]
            )
            for sample_set in meta["sample_sets"]
        ]
        meta["sample_f"] = sample_f
        # ---> reports
        # longest_sample_string = max([len(", ".join(sample_set)) for sample_set in meta['sample_sets']]) + 2
        # meta['spacing'] = longest_sample_string if longest_sample_string > meta['spacing'] else meta['spacing']

    def _read_variants(self, measure_key, vcf_f):
        meta = self._get_meta(measure_key)
        seq_names = meta["seq_names"]
        samples = meta["samples"]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            gt_key, pos_key, sample_key = "calldata/GT", "variants/POS", "samples"
            samples_gt_order = allel.read_vcf(vcf_f, fields=[sample_key])[sample_key]
            query_samples = ordered_intersect(a=samples_gt_order, b=samples, order="a")
            # Check if all samples were found
            if set(query_samples) != set(samples):
                sys.exit(
                    "[X] The following samples in SAMPLE_FILE were not found in VCF_FILE: %s"
                    % (", ".join(list(set(samples).difference(set(query_samples)))))
                )
            # Set up counts arrays
            count_shape = (len(meta["seq_names"]), len(query_samples))
            count_records = np.zeros(count_shape[0], dtype=np.int64)
            count_called = np.zeros(count_shape, dtype=np.int64)
            count_hom_ref = np.zeros(count_shape, dtype=np.int64)
            count_hom_alt = np.zeros(count_shape, dtype=np.int64)
            count_het = np.zeros(count_shape, dtype=np.int64)
            count_missing = np.zeros(count_shape, dtype=np.int64)
            for idx, seq_name in tqdm(
                enumerate(seq_names),
                total=len(seq_names),
                desc="[%] Reading variants",
                ncols=100,
            ):
                vcf_data = allel.read_vcf(
                    vcf_f,
                    region=seq_name,
                    samples=query_samples,
                    fields=[gt_key, pos_key],
                )
                if vcf_data:
                    # genotypes
                    gt_matrix_raw = vcf_data[gt_key]
                    # counts
                    intervals = self._get_interval_coordinates(seq_name=seq_name)
                    sites = intervals_to_sites(intervals)
                    if sites is not None:
                        # positions in VCF
                        pos_array_raw = check_unique_pos(
                            (vcf_data[pos_key] - 1)
                        )  # port to BED (0-based) coordinates
                        # intersection of VCF and BED intervals
                        interval_mask = np.isin(
                            pos_array_raw, sites, assume_unique=True
                        )
                        gt_matrix = gt_matrix_raw[interval_mask]
                        count_records[idx] = gt_matrix.shape[0]
                        pos_array = pos_array_raw[interval_mask]
                        sa_genotype_matrix = allel.GenotypeArray(gt_matrix)
                        count_called[idx, :] = sa_genotype_matrix.count_called(axis=0)
                        count_hom_ref[idx, :] = sa_genotype_matrix.count_hom_ref(axis=0)
                        count_hom_alt[idx, :] = sa_genotype_matrix.count_hom_alt(axis=0)
                        count_het[idx, :] = sa_genotype_matrix.count_het(axis=0)
                        count_missing[idx, :] = sa_genotype_matrix.count_missing(axis=0)
                        self._save_variants(seq_name, pos_array, gt_matrix)
            meta["variants_idx_by_sample"] = {
                query_sample: idx for idx, query_sample in enumerate(query_samples)
            }
        meta["vcf_f"] = vcf_f
        meta["variants_counts"] = int(
            np.sum(count_records)
        )  # ZARR JSON encoder does not like numpy dtypes
        meta["variants_counts_called"] = [
            int(x) for x in np.sum(count_called, axis=0)
        ]  # ZARR JSON encoder does not like numpy dtypes
        meta["variants_counts_hom_ref"] = [
            int(x) for x in np.sum(count_hom_ref, axis=0)
        ]  # ZARR JSON encoder does not like numpy dtypes
        meta["variants_counts_hom_alt"] = [
            int(x) for x in np.sum(count_hom_alt, axis=0)
        ]  # ZARR JSON encoder does not like numpy dtypes
        meta["variants_counts_het"] = [
            int(x) for x in np.sum(count_het, axis=0)
        ]  # ZARR JSON encoder does not like numpy dtypes
        meta["variants_counts_missing"] = [
            int(x) for x in np.sum(count_missing, axis=0)
        ]  # ZARR JSON encoder does not like numpy dtypes
        # QC plots

    def _save_variants(self, sequence, pos_array, gt_matrix):
        self.data.create_dataset(
            "seqs/%s/variants/pos" % sequence, data=pos_array, dtype=np.int64
        )
        self.data.create_dataset(
            "seqs/%s/variants/matrix" % sequence, data=gt_matrix, dtype=np.int64
        )

    def _save_variants_meta(self):
        pass

    def _read_intervals(self, measure_key, bed_f):
        meta = self._get_meta(measure_key)
        target_sequences, target_samples = set(meta["seq_names"]), list(meta["samples"])
        intervals_df = parse_intervals(bed_f, target_sequences, target_samples)
        valid_sequences = intervals_df["sequence"].unique()
        intervals_idx_by_sample = {
            sample: idx for idx, sample in enumerate(intervals_df.columns[3:-1])
        }
        intervals_count = len(intervals_df.index)
        intervals_span = int(intervals_df["length"].sum())
        count_bases_samples = np.zeros(
            (len(valid_sequences), len(target_samples)), dtype=np.int64
        )
        for idx, (sequence, _df) in tqdm(
            enumerate(intervals_df.groupby(["sequence"], observed=True)),
            total=len(valid_sequences),
            desc="[%] Reading intervals",
            ncols=100,
        ):
            interval_matrix = _df[target_samples].to_numpy()
            starts = _df["start"].to_numpy()
            ends = _df["end"].to_numpy()
            self._set_intervals(sequence, interval_matrix, starts, ends)
            length_matrix = interval_matrix * _df["length"].to_numpy().reshape(-1, 1)
            count_bases_samples[idx, :] = np.sum(length_matrix, axis=0)
        self._set_intervals_meta(
            bed_f,
            intervals_idx_by_sample,
            count_bases_samples,
            intervals_count,
            intervals_span,
        )

    def _set_intervals(self, sequence, interval_matrix, starts, ends):
        self.data.create_dataset(
            "seqs/%s/intervals/matrix" % sequence, data=interval_matrix
        )
        self.data.create_dataset("seqs/%s/intervals/starts" % sequence, data=starts)
        self.data.create_dataset("seqs/%s/intervals/ends" % sequence, data=ends)

    def _set_intervals_meta(
        self,
        bed_f,
        intervals_idx_by_sample,
        count_bases_samples,
        intervals_count,
        intervals_span,
    ):
        meta_intervals = self._get_meta("seqs")
        meta_intervals["bed_f"] = bed_f
        meta_intervals["intervals_idx_by_sample"] = intervals_idx_by_sample
        meta_intervals["intervals_span_sample"] = [
            int(x) for x in np.sum(count_bases_samples, axis=0)
        ]  # JSON encoder does not like numpy dtypes
        meta_intervals["intervals_count"] = intervals_count
        meta_intervals["intervals_span"] = intervals_span

    def _preflight_blocks(
        self,
        block_length,
        block_span,
        block_max_multiallelic,
        block_max_missing,
        overwrite,
    ):
        """
        [TODO]: allow for multiple block datasets
        """
        if not self.has_stage("measure") and not self.has_stage("parse"):
            sys.exit(
                "[X] Gimble store %r has no data to block. Please run 'gimble parse'." % self.path
            )
        if self.has_stage("blocks"):
            if not overwrite:
                sys.exit(
                    "[X] Gimble store %r already contains blocks.\n[X] These blocks => %r\n[X] Please specify '--force' to overwrite."
                    % (self.path, self.get_stage("blocks"))
                )
            print(
                "[-] Gimble store %r already contains blocks. But these will be overwritten..."
                % (self.path)
            )
            # wipe bsfs, windows, AND meta, since new blocks...
            blocks_meta = self._get_meta("blocks")
            if blocks_meta:
                self._del_data_and_meta("blocks")
                if "blocks_raw_tally_key" in blocks_meta:
                    self._del_data_and_meta(blocks_meta["blocks_raw_tally_key"])
            windows_meta = self._get_meta("windows")
            if windows_meta:
                self._del_data_and_meta("windows")
                if "windows_raw_tally_key" in windows_meta:
                    self._del_data_and_meta(windows_meta["windows_raw_tally_key"])
                    #self._del_data_and_meta(windows_meta["windowsum_raw_tally_key"])
        config = {
            "blocks_key": "blocks",
            "block_length": block_length,
            "block_span": block_span,
            "block_max_multiallelic": block_max_multiallelic,
            "block_max_missing": block_max_missing,
            "blocks_raw_by_sample_set_idx": collections.Counter(),  # all possible blocks
            "blocks_by_sample_set_idx": collections.Counter(),  # all valid blocks => only these get saved to store
            "blocks_by_sequence": collections.Counter(),  # all valid blocks
            "overwrite": overwrite,
        }
        return config

    def blocks(
        self,
        block_length=64,
        block_span=128,
        block_max_multiallelic=3,
        block_max_missing=3,
        overwrite=False,
    ):
        config = self._preflight_blocks(
            block_length,
            block_span,
            block_max_multiallelic,
            block_max_missing,
            overwrite,
        )
        config = self._make_blocks(config)
        if config["count_total"] == 0:
            sys.exit("[X] No blocks could be generated from data given the parameters.")
        meta = {
            "length": config["block_length"],
            "span": config["block_span"],
            "max_missing": config["block_max_missing"],
            "max_multiallelic": config["block_max_multiallelic"],
            "count_by_sample_set_idx": dict(config["blocks_by_sample_set_idx"]),  # keys are strings
            "count_by_sequence": dict(config["blocks_by_sequence"]),  # keys are strings
            "count_raw_by_sample_set_idx": dict(config["blocks_raw_by_sample_set_idx"]),  # keys are strings
            "count_total": config["count_total"],
            "count_total_raw": config["count_total_raw"],
        }
        self._set_meta(config["blocks_key"], meta=meta)
        # after making blocks, make tally 'raw' without kmax. This is needed for heterozygosity, d_xy, F_sts metrics, etc
        meta["blocks_raw_tally_key"] = self.tally(
            "blocks", "blocks_raw", None, "X", None, None, overwrite=True, verbose=False
        )
        self._set_meta(config["blocks_key"], meta=meta)
        meta_seqs = self._get_meta(config["blocks_key"])

    def windows(self, window_size=500, window_step=100, overwrite=False):
        config = self._preflight_windows(window_size, window_step, overwrite)
        config = self._make_windows(config)
        if config["window_count"] == 0:
            sys.exit("[X] No windows could be made.")
        else:
            meta = {
                "window_size": config["window_size"],
                "window_step": config["window_step"],
                "sample_sets": config["sample_sets"],
                "windows_key": config["windows_key"],
                "window_count": config["window_count"],
            }
            self._set_meta(config["windows_key"], meta=meta)
            meta["windows_raw_tally_key"] = self.tally(
                "windows",
                "windows_raw",
                None,
                "X",
                None,
                None,
                overwrite=True,
                verbose=False,
                tally_form='tally'
            )
            # meta['windowsum_raw_tally_key'] = self.tally('windows', 'windowsum_raw', None, 'X', None, None, overwrite=True, verbose=False)
            self._set_meta(config["windows_key"], meta=meta)

    def _preflight_simulate(self, kwargs):
        kwargs["simulate_key"] = "simulate/%s" % kwargs["simulate_label"]
        if self._has_key(kwargs["simulate_key"]):
            if not kwargs["overwrite"]: 
                sys.exit(
                    "[X] Simulated results with label %r already exist. Use '-f' to overwrite."
                    % (kwargs["simulate_key"])
                )
            self._del_data_and_meta(kwargs["simulate_key"])
        if kwargs["gridsearch_key"]:
            gridsearch_key = self.validate_key(kwargs["gridsearch_key"], 'gridsearch')
            meta_gridsearch = self._get_meta(gridsearch_key)
            if meta_gridsearch['data_source'] == 'sims':
                sys.exit("[X] Gridsearch results stored under %r are based on simulated data %r and can't be used as input for 'simulate'." % (
                    gridsearch_key, meta_gridsearch['data_key']))
            kwargs['demographies'] = self.get_demographies_from_gridsearch(gridsearch_key, kwargs["constraint"])
            if not kwargs['windows'] == len(kwargs['demographies']):
                sys.exit("[X] Specified '--windows' (%s) MUST match windows in gridsearch results (%s intervals)." % (kwargs['windows'], len(kwargs['demographies'])))
        else:
            kwargs['demographies'] = [GimbleDemographyInstance(
                model=kwargs['model'], 
                Ne_A=kwargs['Ne_A'],
                Ne_B=kwargs['Ne_B'],
                Ne_A_B=kwargs['Ne_A_B'],
                me=kwargs['me'], 
                T=kwargs['T'], 
                mu=kwargs['mu'])] * kwargs['windows']
        if kwargs['rec_map']:
            # format is "chr", "start", "end", "rec" (rec in cM/Mb !!!)
            recombination_df = pd.read_csv(kwargs['rec_map'], sep="\t", names=["chr", "start", "end", "rec"], header=None)
            # from cM/Mb to rec/bp (simulator needs rec/b) 
            recombination_df["rec_scaled"] = recombination_df["rec"] * RECOMBINATION_SCALER
            kwargs['recombination_rate'] = tuple(recombination_df['rec_scaled'])
            if not kwargs['windows'] == len(kwargs['recombination_rate']):
                sys.exit("[X] Specified '--windows' (%s) MUST match values in '--rec_map' (%s intervals)." % (kwargs['windows'], len(kwargs['recombination_rate'])))
        else:
            kwargs['recombination_rate'] = [kwargs['rec_rate'] * RECOMBINATION_SCALER] * kwargs['windows']            
        #print("recombination_rates", len(kwargs['recombination_rate']))
                
        print("[+] Simulating %s replicates of %s window(s) of %s blocks" % 
            (kwargs['replicates'],
            kwargs['windows'],
            kwargs['blocks'] 
            ))
        rng = np.random.default_rng(kwargs['seed'])
        kwargs['ancestry_seeds_by_replicate'] = {replicate_idx: rng.integers(1, 2**32, kwargs['windows']) for replicate_idx in range(kwargs['replicates'])}
        kwargs['mutation_seeds_by_replicate'] = {replicate_idx: rng.integers(1, 2**32, kwargs['windows']) for replicate_idx in range(kwargs['replicates'])}
        return kwargs

    def get_demographies_from_gridsearch(self, gridsearch_label, gridsearch_constraint={}):
        print("[+] Gridsearch constraint: %s" % (",".join(
            ["%s=%s" % (k, v) for k, v in gridsearch_constraint.items()]) if gridsearch_constraint else None))
        meta_gridsearch = self._get_meta(gridsearch_label)
        meta_makegrid = self._get_meta(meta_gridsearch["makegrid_key"])
        model = meta_makegrid['model']
        modelObjs = []
        parameter_names = list(meta_gridsearch["grid_dict"].keys())
        parameter_array = np.array([np.array(v, dtype=np.float64) for k, v in meta_gridsearch["grid_dict"].items()]).T
        print("[+] Generating input for 'simulate' from gridsearch results...")
        for gridsearch_key in meta_gridsearch["gridsearch_keys"]: # this loop is because data tried to emulate sims (which was stupid) 
            lncls = np.array(self._get_data(gridsearch_key))  # (w, gridpoints)
            if gridsearch_constraint: 
                fixed_param_index = self._get_fixed_param_index(gridsearch_constraint, parameter_names, parameter_array, meta_makegrid) 
                lncls_fixed = lncls[:,fixed_param_index]
                lncls_fixed_max_idx = np.argmax(lncls_fixed, axis=1)
                lncls_max = lncls_fixed[np.arange(lncls_fixed.shape[0]), lncls_fixed_max_idx]
                lncls_max_parameters = parameter_array[fixed_param_index][lncls_fixed_max_idx]
            else:
                lncls_max_idx = np.argmax(lncls, axis=1)
                lncls_max = lncls[np.arange(lncls.shape[0]), lncls_max_idx]
                lncls_max_parameters = parameter_array[lncls_max_idx]
            for idx in tqdm(range(lncls_max_parameters.shape[0]), desc="[%] Generating", ncols=100):
                model_dict = {parameter_name: parameter_value for parameter_name, parameter_value in zip(parameter_names, lncls_max_parameters[idx,:])}
                modelObj = GimbleDemographyInstance(model=meta_makegrid['model'], **model_dict) # model comes from makegrid
                modelObjs.append(modelObj)
        if len(modelObjs) == 0:
            sys.exit("[X] get_demographies_from_gridsearch() was unsuccessful...")
        return modelObjs

    def simulate(self, **kwargs):
        kwargs = self._preflight_simulate(kwargs)
        simulate_jobs_by_replicate_idx = lib.simulate.get_sim_args_by_replicate_idx(kwargs)
        # create empty arrays in zarr store
        print("[+] Running simulations...")
        with tqdm(total=(kwargs['windows'] * kwargs['replicates']), desc="[%] Simulating", ncols=100) as pbar:
            for replicate_idx in range(kwargs['replicates']):
                replicate_key = "%s/%s" % (kwargs['simulate_key'], replicate_idx)
                replicate_shape = tuple([kwargs['windows']] + list(kwargs['kmax'] + 2))
                #self.data.create_dataset(key, data=zarr.zeros(pileup_shape, chunks=(chunk_size, len(fields)), dtype='i8'), overwrite=True)
                self.data.create_dataset(replicate_key, data=zarr.zeros(replicate_shape), dtype='i8', overwrite=True)
                if kwargs['processes'] <= 1:
                    for simulate_job in simulate_jobs_by_replicate_idx[replicate_idx]:
                        simulate_window = lib.simulate.simulate_call(simulate_job)
                        self.data[replicate_key][simulate_window['window_idx']] = simulate_window['bsfs']
                        pbar.update(1)
                else:
                    with poolcontext(processes=kwargs['processes']) as pool:
                        for simulate_window in pool.imap_unordered(lib.simulate.simulate_call, simulate_jobs_by_replicate_idx[replicate_idx]):
                            self.data[replicate_key][simulate_window['window_idx']] = simulate_window['bsfs']
                            pbar.update(1)
                kwargs['idx'] = replicate_idx
                replicate_key = "%s/%s" % (kwargs['simulate_key'], replicate_idx)
                replicate_meta = {
                    "idx": kwargs['idx'],
                    "ancestry_seeds": tuple([int(s) for s in kwargs["ancestry_seeds_by_replicate"][kwargs["idx"]]]),
                    "mutation_seeds": tuple([int(s) for s in kwargs["mutation_seeds_by_replicate"][kwargs["idx"]]])
                }
                self.data[replicate_key].attrs.put(replicate_meta)
        simulate_meta = {
            'simulate_key': kwargs['simulate_key'],
            'simulate_label': kwargs['simulate_label'],
            'data_source': 'sims',
            'samples_A': kwargs['samples_A'],
            'samples_B': kwargs['samples_A'],
            'replicates': kwargs['replicates'],
            'windows': kwargs['windows'],
            'blocks': kwargs['blocks'],
            'block_length': kwargs['block_length'],
            'continuous_genome': kwargs['continuous_genome'],
            'max_k': list([int(k) for k in kwargs['kmax']]),
            'mu': kwargs['mu'],
            'model': kwargs['model'],
            'gridsearch_key': kwargs['gridsearch_key'],
            'constraint': kwargs['constraint'],
            'rec_rate': kwargs['rec_rate'],
            'rec_map': kwargs['rec_map'],
            'recombination_rate': kwargs['recombination_rate'],
            'Ne_A': kwargs['Ne_A'],
            'Ne_B': kwargs['Ne_B'],
            'Ne_A_B': kwargs['Ne_A_B'],
            'T': kwargs['T'],
            'me': kwargs['me'],
            'processes': kwargs['processes'],
            'seed': kwargs['seed'],
            'overwrite': kwargs['overwrite'],
            'grid_dict': {k: list(v) for k, v in LOD_to_DOL(get_parameter_dicts_from_user_parameters(Ne_A=[kwargs['Ne_A']], Ne_B=[kwargs['Ne_B']], Ne_A_B=[kwargs['Ne_A_B']], T=[kwargs['T']], me=[kwargs['me']])).items()},
        }
        self.data[kwargs['simulate_key']].attrs.put(simulate_meta)
        print("[+] Tally of simulation can be accessed with %r" % kwargs['simulate_key'])

    def _preflight_query(self, version, data_key, extended, fixed_param, sliced_param, diss):
        data_key = self.validate_key(data_key)
        config = {"data_key": data_key, "version": version}
        if not self._has_key(data_key):
            sys.exit(
                "[X] ZARR store %s has no data under the key %r."
                % (self.path, data_key)
            )
        config["data_type"] = data_key.split("/")[0]
        config["extended"] = extended
        config["fixed_param"] = fixed_param
        config["sliced_param"] = sliced_param
        return config

    def _format_grid(self, meta_makegrid):
        grid_dict = meta_makegrid['grid_dict']
        parameters = {k:meta_makegrid[k] for k in grid_dict.keys()}
        grid_points = meta_makegrid['parameters_grid_points']
        output = [
            "[#]\t Makegrid : %s" % meta_makegrid['makegrid_key'], 
            "[#]\t Gridpoints: %s" % grid_points]
        output.append("[#]\t Parameters:")
        for parameter, string in parameters.items():
            output.append("[#]\t\t[%s]" % (parameter))
            output.append("[#]\t\t- config := %s" % ", ".join([str(s) for s in string]))
            output.append("[#]\t\t- values := %s" % (sorted(set(grid_dict[parameter]))))
        return "\n".join(output)

    def query(self, version, data_key, extended, fixed_param, sliced_param, diss):
        config = self._preflight_query(version, data_key, extended, fixed_param, sliced_param, diss)
        if config["data_type"] == "tally":
            self._write_tally_tsv(config)
        elif config["data_type"] == "optimize":
            optimize_meta = self._get_meta(config["data_key"])
            print(format_query_meta(optimize_meta))
            self._write_optimize_tsv(config)
            if not optimize_meta['deme'] == "N/A":
                self._write_demes_yaml(config)
        elif config["data_type"] == "windows":
            self._write_window_bed(config)
        elif config["data_type"] == "gridsearch":
            self._write_gridsearch_results(config)
            #self._write_gridsearch_results_old(config)
        elif config["data_type"] == "makegrid":
            self._write_makegrid(config)
        elif config["data_type"] == "simulate":
            self._write_tally_tsv(config)
        elif config["data_type"] == "diss":
            self._write_diss_tsv(config)
        else:
            sys.exit("[X] Not implemented.")

    def _write_makegrid(self, config):
        '''should create folder with TSVs of each gridpoint'''
        def make_grid_2d(grid):
            idxs = np.nonzero(grid>-10) # should capture all values..
            idxs_array = np.array(idxs).T
            rows = idxs_array.shape[0] # number of rows in array (based on gridpoints and kmax)
            first = idxs_array[:,0].reshape(rows, 1) # gridpoint idx
            second = grid[idxs].reshape(rows, 1) # float in grid
            third = idxs_array[:,1:] # mutuple
            return np.concatenate([first, second, third], axis=1)
        meta_makegrid = self._get_meta(config["data_key"])
        print(self._format_grid(meta_makegrid))
        grid = np.array(self._get_data(meta_makegrid['makegrid_key']))
        grid_2d = make_grid_2d(grid)
        parameter_names = list(meta_makegrid["grid_dict"].keys())
        parameter_array = np.array([np.array(v, dtype=np.float64) for k, v in meta_makegrid["grid_dict"].items()]).T
        dtypes = {"grid_idx": "int64","P": "float64","het_b": "int64","het_a": "int64","het_ab": "int64","fixed_diff": "int64"}
        columns = ["grid_idx", "P", "het_b", "het_a", "het_ab", "fixed_diff"]
        for idx in range(parameter_array.shape[0]):
            fn = "%s.%s.tsv" % (meta_makegrid['makegrid_label'], ".".join(["%s=%s" % (name, float(value)) for name, value in zip(parameter_names, parameter_array[idx])]))
            pd.DataFrame(data=grid_2d[grid_2d[:,0]==idx], columns=columns).astype(dtype=dtypes).to_csv(fn, header=True, index=False, sep="\t")
            print("[#] Wrote file %r." % fn)

    def _write_demes_yaml(self, config):
        optimize_meta = dict(self._get_meta(config["data_key"]))
        fn_yaml = "%s.%s.optimize.demes.yaml" % (
                    self.prefix,
                    config["data_key"].replace("/", "_"),
                )
        demes.dump(demes.loads(optimize_meta['deme']), filename=fn_yaml)
        print("[#] Wrote file %r." % fn_yaml)

    def _write_optimize_tsv(self, config):
        optimize_meta = self._get_meta(config["data_key"])
        single_file_flag = len(optimize_meta["optimize_key"]) == 1
        for idx in range(optimize_meta["nlopt_runs"]):
            instance_key = "%s/%s" % (optimize_meta["optimize_key"], idx)
            instance_meta = dict(self._get_meta(instance_key))
            df = pd.DataFrame(data=instance_meta["optimize_results"])
            if single_file_flag:
                fn_tsv = "%s.%s.optimize.tsv" % (
                    self.prefix,
                    config["data_key"].replace("/", "_"),
                )
            else:
                if optimize_meta['data_source'] == 'sims':
                    df.insert(0, 'replicate_idx', idx)
                fn_tsv = "%s.%s.optimize.%s.tsv" % (
                    self.prefix,
                    config["data_key"].replace("/", "_"),
                    idx,
                )
            with open(fn_tsv, "w") as fh_tsv:
                header = ["# %s" % config["version"]]
                fh_tsv.write("\n".join(header) + "\n")
            if single_file_flag:
                df.to_csv(
                    fn_tsv, index=True, mode='a', index_label="window_idx", sep="\t"
                )
            else:
                if len(df) > 1:
                    df.insert(1, 'window_idx', np.arange(len(df)))
                df.to_csv(
                    fn_tsv, index=False, mode='a', sep="\t"
                )
            print("[#] Wrote file %r." % fn_tsv)

    def _write_diss_tsv(self, config):
        meta_seq = self._get_meta("seqs")
        sample_sets = ['X', 'A', 'B']
        # option A
        blocks_fgv_by_sample_set = collections.Counter()
        blocks_by_sample_set = collections.Counter()
        with tqdm(
            total=(len(meta_seq["seq_names"]) * 3),
            desc="[%] Preparing data",
            ncols=100,
            unit_scale=True,
            ) as pbar:
            for sample_set in sample_sets:
                fn = "gimble.blocks.diss.%s.tsv" % sample_set
                #fn_T = "gimble.blocks.diss.tally.%s.tsv" % sample_set
                diss_dfs = []
                tally_dfs = []
                for seq_name in meta_seq['seq_names']:
                    tally = tally_variation(self._get_variation(data_type="blocks", sequences=[seq_name], sample_sets=sample_set), form="tally")
                    #tally_df = pd.DataFrame({'count': tally[:,0], 'm1': tally[:,1], 'm2': tally[:,2], 'm3': tally[:,3], 'm4': tally[:,4]})
                    #tally_dfs.append(tally_df)
                    blocks_fgv_by_sample_set[sample_set] += np.sum(tally[(tally[:, 3] > 0) & (tally[:, 4] > 0)][:,0])
                    blocks_by_sample_set[sample_set] += np.sum(tally[:,0])
                    counts = tally[:,0]
                    mutations = np.sum(tally[:,1:], axis=1)
                    diss_df = pd.DataFrame({'count': counts, 'mutations': mutations}).groupby(['mutations']).sum().reset_index()
                    diss_df['seq'] = seq_name
                    diss_df = diss_df[['seq', 'count', 'mutations']]
                    diss_dfs.append(diss_df)
                    pbar.update()
                pd.concat(diss_dfs).to_csv(fn, index=False, sep="\t")
                #pd.concat(tally_dfs).to_csv(fn_T, index=False, sep="\t")
        print("[+] FGV blocks for each sample set:")
        for sample_set in sample_sets:
            print("[+]\t%s: %s / %s = %s" % (
                sample_set, 
                int(blocks_fgv_by_sample_set[sample_set]), 
                int(blocks_by_sample_set[sample_set]), 
                (blocks_fgv_by_sample_set[sample_set] / blocks_by_sample_set[sample_set])))
        
        # option B
        #for sample_set in ['X', 'A', 'B']:
        #    fn = "gimble.blocks.diss.option_B.%s.tsv" % sample_set
        #    diss_dfs = []
        #    for seq_name in meta_seq['seq_names']:
        #        tally = tally_variation(self._get_variation(data_type="blocks", sequences=[seq_name], sample_sets="X"), form="tally")
        #        counts = tally[:,0]
        #        mutations = np.sum(tally[:,1:], axis=1)
        #        diss_df = pd.DataFrame({'seq': seq_name, 'count': counts, 'mutations': mutations})
        #        diss_dfs.append(diss_df)
        #    pd.concat(diss_dfs).groupby(['seq', 'mutations']).sum().reset_index().to_csv(fn, index=False, sep="\t")
        #    print("[+] Wrote file %r." % fn)
        
        

    def _write_tally_tsv(self, config):
        data = self._get_data(config["data_key"])
        if not isinstance(data, zarr.core.Array):
            """intermediate fix until we figure out how to retrieve simulate tallies properly..."""
            sys.exit("[X] %r does not point to a tally." % config["data_key"])
        config["data"] = np.array(data)
        # print('config["data"]', type(config["data"]), config["data"].shape, config["data"])
        config["meta"] = dict(self._get_meta(config["data_key"]))
        print("[+] Tally ...")
        print(format_query_meta(config["meta"]))
        config["header"] = ["count", "het_b", "het_a", "het_ab", "fixed_diff"]
        if config["data"].ndim == 5:
            config["header"] = ["window_idx"] + config["header"]
            #if config["data_type"] == "tally":
            #    config["header"] = ["window_idx"] + config["header"]
            #else:
            #    config["header"] = ["window_idx"] + config["header"]
        config["filename"] = "%s.%s.tsv" % (
            self.prefix,
            config["data_key"].replace("/", "_"),
        )
        pd.DataFrame(
            data=bsfs_to_2d(config["data"]), columns=config["header"], dtype="int64"
        ).to_csv(config["filename"], index=False, sep="\t")
        print("[#] Wrote file %r." % config["filename"])

    def _write_window_bed(self, config):
        meta_blocks = self._get_meta("blocks")
        # print(dict(meta_blocks))
        meta_windows = self._get_meta("windows")
        query_meta_blocks = format_query_meta(meta_blocks, ignore_long=True)
        # print(query_meta_blocks)
        query_meta_windows = format_query_meta(meta_windows, ignore_long=True)
        # print(query_meta_windows)
        # print("[+] Getting data for BED ...")
        (
            sequences,
            starts,
            ends,
            index,
            pos_mean,
            pos_median,
            balance,
            mse_sample_set_cov,
        ) = self._get_window_bed()
        # print("[+] Calculating popgen metrics ...")
        tally = self._get_data(meta_windows["windows_raw_tally_key"])
        meta_tally = self._get_meta(meta_windows["windows_raw_tally_key"])

        pop_metrics = get_popgen_metrics(
            tally, sites=(meta_tally["block_length"] * meta_windows["window_size"])
        )
        int_bed = np.vstack(
            [
                starts,
                ends,
                index,
                pos_mean,
                pos_median,
                pop_metrics,
                balance,
                mse_sample_set_cov,
            ]
        ).T
        columns = [
            "sequence",
            "start",
            "end",
            "index",
            "pos_mean",
            "pos_median",
            "heterozygosity_A",
            "heterozygosity_B",
            "d_xy",
            "f_st",
            "balance",
            "mse_sample_set_cov",
        ]
        header = ["# %s" % config["version"], "# %s" % "\t".join(columns)]
        out_f = "%s.windows.popgen.bed" % self.prefix
        print("[+] Writing BED file %r ..." % out_f)
        with open(out_f, "w") as out_fh:
            out_fh.write("\n".join(header) + "\n")
        dtypes = {
            "start": "int64",
            "end": "int64",
            "index": "int64",
            "pos_mean": "int64",
            "pos_median": "int64",
            "heterozygosity_A": "float64",
            "heterozygosity_B": "float64",
            "d_xy": "float64",
            "f_st": "float64",
            "balance": "float64",
        }
        bed_df = pd.DataFrame(data=int_bed, columns=columns[1:]).astype(dtype=dtypes)
        bed_df["sequence"] = sequences
        # bed_df.sort_values(['sequence', 'start'], ascending=[True, True]).to_csv(out_f, na_rep='NA', mode='a', sep='\t', index=False, header=False, columns=columns, float_format='%.5f')
        bed_df.to_csv(
            out_f,
            na_rep="NA",
            mode="a",
            sep="\t",
            index=False,
            header=False,
            columns=columns,
            float_format="%.5f",
        )
        return out_f

    def _get_window_bed_columns(self):
        meta_seqs = self._get_meta("seqs")
        meta_windows = self._get_meta("windows")
        MAX_SEQNAME_LENGTH = max([len(seq_name) for seq_name in meta_seqs["seq_names"]])
        sequences = np.zeros(meta_windows["window_count"], dtype="<U%s" % MAX_SEQNAME_LENGTH)
        starts = np.zeros(meta_windows["window_count"], dtype=np.int64)
        ends = np.zeros(meta_windows["window_count"], dtype=np.int64)
        offset = 0
        for seq_name in tqdm(
            meta_seqs["seq_names"],
            total=len(meta_seqs["seq_names"]),
            desc="[%] Preparing output",
            ncols=100,
            unit_scale=True,
        ):
            start_key = "windows/%s/starts" % (seq_name)
            end_key = "windows/%s/ends" % (seq_name)
            if start_key in self.data:
                start_array = np.array(self.data[start_key])
                window_count = start_array.shape[0]
                starts[offset : offset + window_count] = start_array
                ends[offset : offset + window_count] = np.array(self.data[end_key])
                sequences[offset : offset + window_count] = np.full_like(
                    window_count, seq_name, dtype="<U%s" % MAX_SEQNAME_LENGTH
                )
                offset += window_count
        index = np.arange(meta_windows["window_count"], dtype=np.int64)
        return (sequences, starts, ends, index)

    def _get_indices_by_sliced_param_value(self, config, parameter_names, parameter_array):
        try:
            sliced_param_idx = parameter_names.index(config["sliced_param"])
        except ValueError:
            sys.exit(
                "[X] Parameter %r is not part of the grid %r\n\tAvailable parameters are: %s"
                % (
                    config["sliced_param"],
                    config["data_key"],
                    ", ".join(parameter_names),
                )
            ) 
        def get_lists_of_indices(sliced_param_array):
            idx_sort = np.argsort(sliced_param_array)
            sorted_sliced_param_array = sliced_param_array[idx_sort]
            vals, idx_start, count = np.unique(
                sorted_sliced_param_array,
                return_counts=True,
                return_index=True,
            )
            res = np.split(idx_sort, idx_start[1:])
            return {"%s=%s" % (config["sliced_param"], value): indices for value, indices in dict(zip(vals, res)).items()}

        indices_by_sliced_param_value = get_lists_of_indices(parameter_array[:, sliced_param_idx])
        return indices_by_sliced_param_value

    def _get_fixed_param_index(self, fixed_params_dict, parameter_names, parameter_array, meta_makegrid):
        fixed_param_mask = np.zeros(parameter_array.shape, dtype=bool)
        for param, value in fixed_params_dict.items(): 
            try:
                fixed_param_idx = parameter_names.index(param)
                fixed_param_mask[:,fixed_param_idx] = (parameter_array[:,fixed_param_idx]==value)
            except ValueError:
                print("[X] --fixed_param %r not part of grid %r." % (param, meta_makegrid['makegrid_label']))
        fixed_param_index = np.zeros(parameter_array.shape[0], dtype=bool)
        # only those rows which have as many "True" values as parameters were specified... 
        fixed_param_index = np.sum(fixed_param_mask,axis=1)==len(fixed_params_dict)
        if not np.any(fixed_param_index):
            sys.exit(
                "[X] Specified constraint %r not found in grid:\n%s" % (
                    ", ".join(["%s=%s" % (param, value) for param, value in fixed_params_dict.items()]), 
                    self._format_grid(meta_makegrid))
            )
        return fixed_param_index

    def _write_sliced_gridsearch_bed(self, config, header, indices_by_sliced_param_value, base_fn, lncls, bed_columns, sequence_array=None, start_array=None, end_array=None, index_array=None):
        sliced_columns = list(indices_by_sliced_param_value.keys())
        sliced_lncls = np.zeros(
            (lncls.shape[0], len(indices_by_sliced_param_value.keys()))
        )
        for sliced_param_idx, (
            sliced_param_value,
            sliced_param_indices,
        ) in enumerate(indices_by_sliced_param_value.items()):
            sliced_lncls[:, sliced_param_idx] = np.max(
                    lncls[:, sliced_param_indices], axis=-1
            )
        if len(bed_columns) == 1:
            sliced_df = pd.DataFrame(
                data=np.vstack([
                    base_fn,
                    sliced_lncls.T]).T, columns=bed_columns + sliced_columns)
        else:
            sliced_df = pd.DataFrame(
                data=np.vstack([
                    sequence_array, 
                    start_array, 
                    end_array, 
                    index_array, 
                    sliced_lncls.T]).T, columns=bed_columns + sliced_columns)
        sliced_out_f = "%s.sliced_param.%s.bed" % (
            base_fn, config["sliced_param"],
        )
        with open(sliced_out_f, "w") as sliced_out_fh:
            sliced_out_header = header + ["# %s" % "\t".join(bed_columns + sliced_columns)]
            sliced_out_fh.write("\n".join(sliced_out_header) + "\n")
        sliced_df.to_csv(
            sliced_out_f,
            na_rep="NA",
            mode="a",
            sep="\t",
            index=False,
            header=False
        )
        print("[+] \tWrote %s ..." % sliced_out_f)

    def _write_fixed_gridsearch_bed(self, config, header, parameter_columns, base_fn, lncls, columns, fixed_param_index, parameter_array, sequence_array=None, start_array=None, end_array=None, index_array=None):
        # subset only those that coincide with fixed param ... 
        lncls_fixed = lncls[:,fixed_param_index]
        lncls_fixed_max_idx = np.argmax(lncls_fixed, axis=1)
        lncls_fixed_max = lncls_fixed[np.arange(lncls_fixed.shape[0]), lncls_fixed_max_idx]
        lncls_fixed_max_parameters = parameter_array[fixed_param_index][lncls_fixed_max_idx]
        if len(columns) == 1: 
            fixed_param_df = pd.DataFrame(
                data=np.vstack([
                    base_fn,
                    lncls_fixed_max_parameters.T, 
                    lncls_fixed_max]).T, columns=columns + parameter_columns)
        else:
            fixed_param_df = pd.DataFrame(
                data=np.vstack([
                    sequence_array, 
                    start_array, 
                    end_array, 
                    index_array, 
                    lncls_fixed_max_parameters.T, 
                    lncls_fixed_max]).T, columns=columns + parameter_columns)
        fixed_param_out_f = "%s.fixed_param.%s.bed" % (base_fn, ".".join(["%s_%s" % (param, str(value).replace(".", "_")) for param, value in config["fixed_param"].items()]))
        with open(fixed_param_out_f, "w") as fixed_param_out_fh:
            fixed_param_out_header = header + ["# --fixed_param %s" % ",".join(["%s=%s" % (param, value) for param, value in config["fixed_param"].items()])] + ["# %s" % "\t".join(columns + parameter_columns)]
            fixed_param_out_fh.write("\n".join(fixed_param_out_header) + "\n")
        fixed_param_df.to_csv(
            fixed_param_out_f,
            na_rep="NA",
            mode="a",
            sep="\t",
            index=False,
            header=False,
        )
        print("[+] \tWrote %s ..." % fixed_param_out_f)

    def _write_sliced_gridsearch_sims(self, header, lncls, indices_by_sliced_param_value, columns, base_fn, config, replicate_idx_array=None, window_idx_array=None, recombination_rate_array=None):
        sliced_columns = list(indices_by_sliced_param_value.keys())
        sliced_lncls = np.zeros(
            (lncls.shape[0], len(indices_by_sliced_param_value.keys()))
        )
        for sliced_param_idx, (
            sliced_param_value,
            sliced_param_indices,
        ) in enumerate(indices_by_sliced_param_value.items()):
            sliced_lncls[:, sliced_param_idx] = np.max(
                    lncls[:, sliced_param_indices], axis=-1
            )
        if len(columns) == 1:
            sliced_df = pd.DataFrame(
            data=np.vstack([
                base_fn, 
                sliced_lncls.T]).T, columns=columns + sliced_columns)
        else:
            sliced_df = pd.DataFrame(
                data=np.vstack([
                    replicate_idx_array, window_idx_array, recombination_rate_array, 
                    sliced_lncls.T]).T, columns=columns + sliced_columns)
        sliced_out_f = "%s.sliced_param.%s.bed" % (
            base_fn, config["sliced_param"],
        )
        with open(sliced_out_f, "w") as sliced_out_fh:
            sliced_out_header = header + ["# %s" % "\t".join(columns + sliced_columns)]
            sliced_out_fh.write("\n".join(sliced_out_header) + "\n")
        sliced_df.to_csv(
            sliced_out_f,
            na_rep="NA",
            mode="a",
            sep="\t",
            index=False,
            header=False
        )
        print("[+] \tWrote %s ..." % sliced_out_f)

    def _write_fixed_gridsearch_sims(self, config, header, parameter_columns, base_fn, lncls, columns, fixed_param_index, parameter_array, replicate_idx_array=None, window_idx_array=None, recombination_rate_array=None):
        # subset only those that coincide with fixed param ... 
        lncls_fixed = lncls[:,fixed_param_index]
        lncls_fixed_max_idx = np.argmax(lncls_fixed, axis=1)
        lncls_fixed_max = lncls_fixed[np.arange(lncls_fixed.shape[0]), lncls_fixed_max_idx]
        lncls_fixed_max_parameters = parameter_array[fixed_param_index][lncls_fixed_max_idx]
        if len(columns) == 1:
            fixed_param_df = pd.DataFrame(
                data=np.vstack([
                    base_fn,
                    lncls_fixed_max_parameters.T, 
                    lncls_fixed_max]).T, columns=columns + parameter_columns)
        else:
            fixed_param_df = pd.DataFrame(
                data=np.vstack([
                    replicate_idx_array,  
                    window_idx_array, 
                    recombination_rate_array, 
                    lncls_fixed_max_parameters.T, 
                    lncls_fixed_max]).T, columns=columns + parameter_columns)
        fixed_param_out_f = "%s.fixed_param.%s.bed" % (base_fn, ".".join(["%s_%s" % (param, str(value).replace(".", "_")) for param, value in config["fixed_param"].items()]))
        with open(fixed_param_out_f, "w") as fixed_param_out_fh:
            fixed_param_out_header = header + ["# --fixed_param %s" % ",".join(["%s=%s" % (param, value) for param, value in config["fixed_param"].items()])] + ["# %s" % "\t".join(columns + parameter_columns)]
            fixed_param_out_fh.write("\n".join(fixed_param_out_header) + "\n")
        fixed_param_df.to_csv(
            fixed_param_out_f,
            na_rep="NA",
            mode="a",
            sep="\t",
            index=False,
            header=False,
        )
        print("[+] \tWrote %s ..." % fixed_param_out_f)

    def _write_gridsearch_results_sims(self, config):
        '''TEMPORARY'''
        # windows = replicates
        # needs recombination rate added
        # 
        #print("_write_gridsearch_results_sims")
        meta_gridsearch = self._get_meta(config["data_key"])
        #print("\nmeta_gridsearch", dict(meta_gridsearch))
        meta_makegrid = self._get_meta(meta_gridsearch["makegrid_key"])
        #print("\nmeta_makegrid", dict(meta_makegrid))
        meta_tally = self._get_meta(meta_gridsearch["data_key"])
        #print("\nmeta_tally", dict(meta_tally))
        data_ndims = meta_gridsearch["data_ndims"]
        #print("data_ndims", data_ndims)
        parameter_names = list(meta_gridsearch["grid_dict"].keys())
        parameter_array = np.array(
            [
                np.array(v, dtype=np.float64) for k, v in meta_gridsearch["grid_dict"].items()
            ]
        ).T
        grid_points = meta_makegrid.get("parameters_grid_points", parameter_array.shape[0])
        columns = []
        fixed_param_columns = []
        header = ["# %s" % config["version"]]
        parameter_columns = parameter_names + ["lnCl"]
        if config["sliced_param"]:
            indices_by_sliced_param_value = self._get_indices_by_sliced_param_value(config, parameter_names, parameter_array)     
        if config["fixed_param"]:
            fixed_param_index = self._get_fixed_param_index(config["fixed_param"], parameter_names, parameter_array, meta_makegrid)
        if data_ndims == 5:
            replicates = meta_tally['replicates']
            windows = meta_tally['windows']
            window_idx_array = np.arange(windows)
            recombination_rate_array = np.array(meta_tally['recombination_rate'])
            columns = ["replicate_idx", "window_idx", "recombination_rate"]
            for replicate_idx, gridsearch_key in enumerate(meta_gridsearch["gridsearch_keys"]):
                replicate_idx_array = np.repeat(replicate_idx, windows)
                base_fn = "%s.%s" % (self.prefix, "_".join(gridsearch_key.split("/")))
                lncls = np.array(self._get_data(gridsearch_key))  # (w, gridpoints)
                if config["sliced_param"]:
                    self._write_sliced_gridsearch_sims(header, lncls, indices_by_sliced_param_value, columns, base_fn, config, replicate_idx_array, window_idx_array, recombination_rate_array)
                if config["fixed_param"]:
                    self._write_fixed_gridsearch_sims(config, header, parameter_columns, base_fn, lncls, columns, fixed_param_index, parameter_array, replicate_idx_array, window_idx_array, recombination_rate_array)

                #print("lncls", lncls.shape, lncls)
                lncls_max_idx = np.argmax(lncls, axis=1)
                lncls_max = lncls[np.arange(lncls.shape[0]), lncls_max_idx]
                lncls_max_parameters = parameter_array[lncls_max_idx]
                lncls_max_df = pd.DataFrame(
                        data=np.vstack([
                            replicate_idx_array,
                            window_idx_array, 
                            recombination_rate_array,
                            lncls_max_parameters.T, 
                            lncls_max]).T, columns=columns + parameter_columns)
                lncls_max_df['replicate_idx'].astype('int64')
                lncls_max_df['window_idx'].astype('int64')
                lncls_max_out_f = "%s.best.bed" % base_fn
                with open(lncls_max_out_f, "w") as lncls_max_out_fh:
                    lncls_max_out_header = header + ["# %s" % "\t".join(columns + parameter_columns)]
                    lncls_max_out_fh.write("\n".join(lncls_max_out_header) + "\n")
                lncls_max_df.to_csv(
                    lncls_max_out_f,
                    na_rep="NA",
                    mode="a",
                    sep="\t",
                    index=False,
                    header=False,
                )
                print("[+] \tWrote %s ..." % lncls_max_out_f)
        if data_ndims == 4:
            replicates = meta_tally['replicates']
            columns = ['dataset']
            for replicate_idx, gridsearch_key in enumerate(meta_gridsearch["gridsearch_keys"]):
                base_fn = "%s.%s" % (self.prefix, "_".join(gridsearch_key.split("/")))
                lncls = np.array(self._get_data(gridsearch_key))  # (w, gridpoints)
                if config["sliced_param"]:
                    self._write_sliced_gridsearch_sims(header, lncls, indices_by_sliced_param_value, columns, base_fn, config)
                if config["fixed_param"]:
                    self._write_fixed_gridsearch_sims(config, header, parameter_columns, base_fn, lncls, columns, fixed_param_index, parameter_array)
                #print("lncls", lncls.shape, lncls)
                lncls_max_idx = np.argmax(lncls, axis=1)
                lncls_max = lncls[np.arange(lncls.shape[0]), lncls_max_idx]
                lncls_max_parameters = parameter_array[lncls_max_idx]
                lncls_max_df = pd.DataFrame(
                        data=np.vstack([
                            base_fn,
                            lncls_max_parameters.T, 
                            lncls_max]).T, columns=columns + parameter_columns)
                lncls_max_out_f = "%s.best.bed" % base_fn
                with open(lncls_max_out_f, "w") as lncls_max_out_fh:
                    lncls_max_out_header = header + ["# %s" % "\t".join(columns + parameter_columns)]
                    lncls_max_out_fh.write("\n".join(lncls_max_out_header) + "\n")
                lncls_max_df.to_csv(
                    lncls_max_out_f,
                    na_rep="NA",
                    mode="a",
                    sep="\t",
                    index=False,
                    header=False,
                )
                print("[+] \tWrote %s ..." % lncls_max_out_f)

    def _write_gridsearch_results(self, config):
        tally_key = config["data_key"] 
        meta_gridsearch = self._get_meta(config["data_key"])
        meta_makegrid = self._get_meta(meta_gridsearch["makegrid_key"])
        meta_tally = self._get_meta(meta_gridsearch["data_key"])
        data_ndims = meta_gridsearch.get("data_ndims", None)
        tally_ndims = meta_tally.get("data_ndims", None)
        data_source = meta_gridsearch.get("data_source", None)
        if not all([data_ndims, data_source]):
            # check whether these fields are set
            sys.exit("[X] The gridsearch result was made with an old version of gIMble. Please re-do the gridsearch.")
        parameter_names = list(meta_gridsearch["grid_dict"].keys())
        parameter_array = np.array([np.array(v, dtype=np.float64) for k, v in meta_gridsearch["grid_dict"].items()]).T
        parameter_columns = parameter_names + ["lnCl"]
        columns = []
        header = ["# %s" % config["version"]]
        if config["sliced_param"]:
            indices_by_sliced_param_value = self._get_indices_by_sliced_param_value(config, parameter_names, parameter_array)     
        if config["fixed_param"]:
            fixed_param_index = self._get_fixed_param_index(config["fixed_param"], parameter_names, parameter_array, meta_makegrid)
        # prepare base_df (on which the remaining part of query gets concatenated onto)
        if data_ndims == 4:
            if data_source == 'meas':
                if tally_ndims == 5: # windows.windowsum
                    columns = ["tally_key", "total_windows_count", "total_blocks_count", 'block_length']
                    base_df = pd.DataFrame({
                        'tally_key': [tally_key],
                        'total_windows_count': [meta_tally['windows']],
                        'total_blocks_count': [meta_tally['blocks']],
                        'block_length': [meta_tally['block_length']],
                        })
                if tally_ndims == 4: # blocks
                    columns = ["tally_key", 'total_blocks_count', 'block_length'] 
                    base_df = pd.DataFrame({
                        'tally_key': [tally_key],
                        'total_blocks_count': [meta_tally['blocks']],
                        'block_length': [meta_tally['block_length']],
                        })
            if data_source == 'sims': # simulate.windowsum
                windows_count = meta_tally['windows']
                row_count = 1
                columns = ["replicate_idx", "total_windows_count", 'total_blocks_count', 'block_length'] 
                replicate_count = meta_tally['replicates']
                base_df = pd.DataFrame({
                        'replicate_idx': np.zeros(row_count),
                        'total_windows_count': np.repeat(meta_tally['windows'], row_count),
                        'total_blocks_count': np.repeat(meta_tally['blocks'], row_count),
                        'block_length': np.repeat(meta_tally['block_length'], row_count)
                        })
        if data_ndims == 5:
            windows_count = meta_tally['windows']
            row_count = windows_count
            window_idx_array = np.arange(windows_count)
            if data_source == 'meas': # windows
                columns = ["sequence", "start", "end", "window_idx"]
                sequence_array, start_array, end_array, index_array = self._get_window_bed_columns()
                base_df = pd.DataFrame({
                    'sequence': sequence_array,
                    'start': start_array,
                    'end': end_array,
                    'window_idx': index_array})
            if data_source == 'sims': # simulate
                columns = ["window_idx", "recombination_rate"]
                recombination_rate_array = np.array(meta_tally['recombination_rate'])
                base_df = pd.DataFrame({
                    'replicate_idx': np.zeros(row_count), 
                    'window_idx': window_idx_array,
                    'recombination_rate': recombination_rate_array})
        # sort out query
        for replicate_idx, gridsearch_key in enumerate(meta_gridsearch["gridsearch_keys"]):
            base_fn = "%s.%s" % (self.prefix, "_".join(gridsearch_key.split("/")))
            lncls = np.array(self._get_data(gridsearch_key))  # (w, gridpoints)
            if data_source == 'sims':
                base_df['replicate_idx'] = np.repeat(replicate_idx, row_count)
            if config["sliced_param"]:
                def _get_sliced_df(lncls, indices_by_sliced_param_value):
                    sliced_columns = list(indices_by_sliced_param_value.keys())
                    sliced_lncls = np.zeros((lncls.shape[0], len(sliced_columns)))
                    for sliced_param_idx, (sliced_param_value, sliced_param_indices) in enumerate(indices_by_sliced_param_value.items()):
                        sliced_lncls[:, sliced_param_idx] = np.max(lncls[:, sliced_param_indices], axis=-1)
                    sliced_df = pd.DataFrame(data=np.vstack([sliced_lncls.T]).T, columns=sliced_columns)
                    return sliced_df
                sliced_fn = "%s.sliced_param.%s.bed" % (base_fn, config["sliced_param"])
                sliced_df = pd.concat([base_df, _get_sliced_df(lncls, indices_by_sliced_param_value)], axis=1)
                with open(sliced_fn, "w") as sliced_fh:
                    sliced_header = header + ["# %s" % "\t".join(list(sliced_df.columns))]
                    sliced_fh.write("\n".join(sliced_header) + "\n")
                sliced_df.to_csv(sliced_fn, na_rep="NA", mode="a", sep="\t", index=False, header=False)
                print("[+] \tWrote %r" % sliced_fn)
            if config["fixed_param"]:
                def _get_fixed_df(lncls, parameter_columns, fixed_param_index):
                    lncls_fixed = lncls[:,fixed_param_index]
                    lncls_fixed_max_idx = np.argmax(lncls_fixed, axis=1)
                    lncls_fixed_max = lncls_fixed[np.arange(lncls_fixed.shape[0]), lncls_fixed_max_idx]
                    lncls_fixed_max_parameters = parameter_array[fixed_param_index][lncls_fixed_max_idx]
                    fixed_df = pd.DataFrame(
                        data=np.vstack([
                            lncls_fixed_max_parameters.T, 
                            lncls_fixed_max]).T, columns=parameter_columns)
                    return fixed_df
                fixed_fn = "%s.fixed_param.%s.bed" % (base_fn, ".".join(["%s_%s" % (param, str(value).replace(".", "_")) for param, value in config["fixed_param"].items()]))
                fixed_df = pd.concat([base_df, _get_fixed_df(lncls, parameter_columns, fixed_param_index)], axis=1)
                with open(fixed_fn, "w") as fixed_fh:
                    fixed_header = header + ["# --fixed_param %s" % ",".join(["%s=%s" % (param, value) for param, value in config["fixed_param"].items()])] + ["# %s" % "\t".join(list(fixed_df.columns))]
                    fixed_fh.write("\n".join(fixed_header) + "\n")
                fixed_df.to_csv(
                    fixed_fn,
                    na_rep="NA",
                    mode="a",
                    sep="\t",
                    index=False,
                    header=False,
                )
                print("[+] \tWrote %r" % fixed_fn)
            def _get_lncls_max_df(lncls, parameter_array):  
                lncls_max_idx = np.argmax(lncls, axis=1)
                lncls_max = lncls[np.arange(lncls.shape[0]), lncls_max_idx]
                lncls_max_parameters = parameter_array[lncls_max_idx]
                lncls_max_df = pd.DataFrame(
                    data=np.vstack([
                        lncls_max_parameters.T, 
                        lncls_max]).T, columns=parameter_columns)
                return lncls_max_df
            lncls_max_fn = "%s.best.bed" % base_fn
            lncls_max_df = pd.concat([base_df, _get_lncls_max_df(lncls, parameter_array)], axis=1)
            with open(lncls_max_fn, "w") as lncls_max_fh:
                lncls_max_header = header + ["# %s" % "\t".join(list(lncls_max_df.columns))]
                lncls_max_fh.write("\n".join(lncls_max_header) + "\n")
            lncls_max_df.to_csv(
                lncls_max_fn,
                na_rep="NA",
                mode="a",
                sep="\t",
                index=False,
                header=False,
            )
            print("[+] \tWrote %r" % lncls_max_fn)

    def _write_gridsearch_results_old(self, config):
        meta_gridsearch = self._get_meta(config["data_key"])
        #print('dict(meta_gridsearch)', dict(meta_gridsearch))
        meta_makegrid = self._get_meta(meta_gridsearch["makegrid_key"])
        #print('dict(meta_makegrid)', dict(meta_makegrid))
        meta_tally = self._get_meta(meta_gridsearch["data_key"])
        #print('dict(meta_tally)', dict(meta_tally))
        data_ndims = meta_gridsearch.get("data_ndims", None)
        if meta_gridsearch['data_source'] == 'sims':
            self._write_gridsearch_results_sims(config)
        parameter_names = list(meta_gridsearch["grid_dict"].keys())
        parameter_array = np.array(
            [
                np.array(v, dtype=np.float64) for k, v in meta_gridsearch["grid_dict"].items()
            ]
        ).T
        grid_points = meta_makegrid.get("parameters_grid_points", parameter_array.shape[0])
        columns = []
        fixed_param_columns = []
        header = ["# %s" % config["version"]]
        if config["sliced_param"]:
            indices_by_sliced_param_value = self._get_indices_by_sliced_param_value(config, parameter_names, parameter_array)     
        if config["fixed_param"]:
            # uses same columns as overall 'best' ...
            fixed_param_index = self._get_fixed_param_index(config["fixed_param"], parameter_names, parameter_array, meta_makegrid)
        parameter_columns = parameter_names + ["lnCl"]

        if data_ndims == 5:
            sequence_array, start_array, end_array, index_array = self._get_window_bed_columns()
            bed_columns = ["sequence", "start", "end", "index"]
            for gridsearch_key in meta_gridsearch["gridsearch_keys"]:
                base_fn = "%s.%s" % (self.prefix, "_".join(gridsearch_key.split("/")))
                lncls = np.array(self._get_data(gridsearch_key))  # (w, gridpoints)
                #print("lncls.shape", lncls.shape)
                #print("lncls", lncls.shape)
                #print("lncls", lncls.shape, lncls)
                if config["sliced_param"]:
                    self._write_sliced_gridsearch_bed(config, header, indices_by_sliced_param_value, base_fn, lncls, bed_columns, sequence_array, start_array, end_array, index_array)
                if config["fixed_param"]:
                    self._write_fixed_gridsearch_bed(config, header, parameter_columns, base_fn, lncls, bed_columns, fixed_param_index, parameter_array, sequence_array, start_array, end_array, index_array)
                lncls_max_idx = np.argmax(lncls, axis=1)
                lncls_max = lncls[np.arange(lncls.shape[0]), lncls_max_idx]
                lncls_max_parameters = parameter_array[lncls_max_idx]
                lncls_max_df = pd.DataFrame(
                        data=np.vstack([
                            sequence_array, 
                            start_array, 
                            end_array, 
                            index_array, 
                            lncls_max_parameters.T, 
                            lncls_max]).T, columns=bed_columns + parameter_columns)
                lncls_max_out_f = "%s.best.bed" % base_fn
                with open(lncls_max_out_f, "w") as lncls_max_out_fh:
                    lncls_max_out_header = header + ["# %s" % "\t".join(bed_columns + parameter_columns)]
                    lncls_max_out_fh.write("\n".join(lncls_max_out_header) + "\n")
                lncls_max_df.to_csv(
                    lncls_max_out_f,
                    na_rep="NA",
                    mode="a",
                    sep="\t",
                    index=False,
                    header=False,
                )
                print("[+] \tWrote %s ..." % lncls_max_out_f)
        
        if data_ndims == 4:
            for gridsearch_key in meta_gridsearch["gridsearch_keys"]:
                base_fn = "%s.%s" % (self.prefix, "_".join(gridsearch_key.split("/")))
                columns = [base_fn]
                lncls = np.array(self._get_data(gridsearch_key))  # (w, gridpoints)
                #print("## lncls", lncls.shape)
                if config["sliced_param"]:
                    self._write_sliced_gridsearch_bed(config, header, indices_by_sliced_param_value, base_fn, lncls, columns)
                if config["fixed_param"]:
                    self._write_fixed_gridsearch_bed(config, header, parameter_columns, base_fn, lncls, columns, fixed_param_index, parameter_array)
                lncls_max_idx = np.argmax(lncls, axis=1)
                #print('lncls_max_idx', lncls_max_idx)
                lncls_max = lncls[np.arange(lncls.shape[0]), lncls_max_idx]
                lncls_max_parameters = parameter_array[lncls_max_idx]
                lncls_max_df = pd.DataFrame(
                        data=np.vstack([
                            base_fn,
                            lncls_max_parameters.T, 
                            lncls_max]).T, columns=columns + parameter_columns)
                lncls_max_out_f = "%s.best.bed" % base_fn
                with open(lncls_max_out_f, "w") as lncls_max_out_fh:
                    lncls_max_out_header = header + ["# %s" % "\t".join(columns + parameter_columns)]
                    lncls_max_out_fh.write("\n".join(lncls_max_out_header) + "\n")
                lncls_max_df.to_csv(
                    lncls_max_out_f,
                    na_rep="NA",
                    mode="a",
                    sep="\t",
                    index=False,
                    header=False,
                )
                print("[+] \tWrote %s ..." % lncls_max_out_f)
        
    def gridsearch_preflight(self, data_key, grid_key, windowsum, overwrite):
        config = {'gridsearch_time': "None"}
        if not self._has_key(grid_key):
            sys.exit("[X] No grid found under key %r." % grid_key)
        config["makegrid_key"] = grid_key
        grid_meta = self._get_meta(config["makegrid_key"])
        grid = np.array(
            self._get_data(config["makegrid_key"]), dtype=np.float64
        )  # grid is likelihoods
        config["grid_dict"] = grid_meta["grid_dict"]  # NEEDED FOR SIMULATE
        config["makegrid_label"] = grid_meta["makegrid_label"]
        config["makegrid_block_length"] = grid_meta["block_length"]
        config["makegrid_kmax"] = grid_meta["max_k"]
        config["makegrid_model"] = grid_meta["model"]
        config["parameters_grid_points"] = grid_meta.get("parameters_grid_points", grid.shape[0])
        # Error if no data
        if not self._has_key(data_key):
            sys.exit("[X] No data found with key %r." % data_key)
        data_meta = self._get_meta(data_key)
        # Error if windowsum on blocks (should eventually be changed to work on ndim's)
        if not data_meta['data_source'] == 'sims' and data_meta['data_type'] == 'blocks' and windowsum:
            sys.exit("[X] '--windowsum' is only applicable to data in windows. Tally in %r consists of blocks." % data_key)
        config["data_key"] = data_key
        config["data_source"] = data_meta['data_source']
        config["data_label"] = "%s%s" % (data_key.split("/")[1], (".windowsum" if windowsum else ""))
        config['data_kmax'] = np.array(data_meta['max_k'])
        config['data_block_length'] = data_meta['block_length']
        config["batch_sites"] = data_meta['block_length'] * data_meta['blocks']
        config['data_replicates'] = data_meta.get('replicates', 1)
        config['data_windows'] = data_meta['windows']
        config['data_blocks'] = data_meta['blocks']
        # optimize
        if config["data_source"] == "sims":
            data = [(replicate_idx, tally) if not windowsum else (replicate_idx, np.sum(tally, axis=0)) for replicate_idx, tally in self.data[config["data_key"]].arrays()]
        else:
            # set FGVs to 0 (needs to be fixed for blocks, not important since FGVs don't come out of agemo unless anomaly)
            # fgv_idxs = get_fgv_idxs(np.array(data_meta['max_k']))
            # tally = self.data[config["data_key"]] if not windowsum else np.sum(self.data[config["data_key"]], axis=0)
            # tally[fgv_idxs] = 0
            # data = [(0, tally)]
            data = [(0, self.data[config["data_key"]]) if not windowsum else (0, np.sum(self.data[config["data_key"]], axis=0))]
        #if sim_key:
        #    sim_meta = self._get_meta(sim_key)
        #    #print('sim_meta', dict(sim_meta))
        #    config["data_key"] = sim_key
        #    config["data_source"] = "sims"
        #    config["data_label"] = sim_meta['simulate_label'] if not windowsum else "%s.windowsum" % sim_meta['simulate_label']
        #    config['data_kmax'] = np.array(sim_meta['max_k'])
        #    config['data_block_length'] = sim_meta['block_length']
        #    config["batch_sites"] = sim_meta['block_length'] * sim_meta['blocks'] # ?
        #    config['data_replicates'] = sim_meta['replicates']
        #    config['data_windows'] = sim_meta['windows']
        #    config['data_blocks'] = sim_meta['blocks']
        #    if windowsum:
        #        data = [(idx, np.sum(tally, axis=0, dtype=np.int64)) for idx, tally in self.data[sim_key].arrays()]
        #    else:
        #        data = [(idx, tally) for idx, tally in self.data[sim_key].arrays()]
        #if tally_key:
        #    tally_meta = self._get_meta(tally_key)
        #    #print('tally_meta', dict(tally_meta))
        #    config["data_source"] = "meas"
        #    config["data_label"] = "/".join(tally_key.split("/")[1:])
        #    config["data_key"] = tally_key
        #    config['data_kmax'] = np.array(tally_meta['max_k'])
        #    config['data_block_length'] = tally_meta['block_length']
        #    config['data_replicates'] = 1
        #    config['data_windows'] = tally_meta['windows']
        #    config['data_blocks'] = tally_meta['blocks']
        #    data = [(0, self._get_data(tally_key)) for _ in (0,)]
        config['data_ndims'] = data[0][1].ndim
        config["batch_sites"] = config['data_block_length'] * config['data_blocks'] # ?
        config["gridsearch_key"] = "gridsearch/%s/%s" % (config["data_label"], config["makegrid_label"])
        config["gridsearch_keys"] = ["gridsearch/%s/%s/%s" % (config["data_label"], config['makegrid_label'], idx)
                for idx in range(config["data_replicates"])]
        if self._has_key(config["gridsearch_key"]):
            if not overwrite: 
                sys.exit(
                    "[X] Gridsearch results with grid label %r on data %r already exist. Use '-f' to overwrite."
                    % (config["makegrid_label"], config["data_label"])
                )
            self._del_data_and_meta(config["gridsearch_key"])
        if not config["data_block_length"] == config["makegrid_block_length"]:
            sys.exit(
                "[X] Block lengths in data %r (%s) and grid %r (%s) are not compatible."
                % (
                    config["data_label"], 
                    config["data_block_length"], 
                    config["makegrid_label"], 
                    config["makegrid_block_length"]))
        if not np.array_equal(config["data_kmax"], config["makegrid_kmax"]):
            sys.exit(
                "[X] k-max in data %r (%s) and grid %r (%s) are not compatible."
                % (
                    config["data_key"], 
                    ",".join(map(str, config["data_kmax"])),
                    config["makegrid_key"], 
                    ",".join(map(str, config["makegrid_kmax"]))))
        return (config, data, grid)

    def gridsearch(
        #self, tally_key, sim_key, grid_key, windowsum, num_cores, chunksize, overwrite
        self, data_key, grid_key, windowsum, num_cores, chunksize, overwrite
    ):
        print("[+] Gathering data for gridsearch ...")
        config, data, grid = self.gridsearch_preflight(
            #tally_key, sim_key, grid_key, windowsum, overwrite
            data_key, grid_key, windowsum, overwrite
        )
        print(
            "[+] Searching with grid %r along tally %r ..."
            % (config["makegrid_label"], config["data_label"])
        )
        position = 1
        global_desc = "[%] Gridsearch"
        is_single_replicate = (len(data)==1)
        with tqdm(total=len(data), position=position, desc=global_desc, ncols=100, leave=True, disable=is_single_replicate) as pbar:
            for i in range(len(data)):
                #print("tally", type(tally), tally.shape)
                idx, tally = data[i]
                iteration_desc = (global_desc if is_single_replicate else "[%%] Replicate %s" % str(i))
                gridsearch_dask_result = gridsearch_dask(
                    tally=tally, grid=grid, num_cores=num_cores, chunksize=chunksize, desc=iteration_desc)
                position+=1
                self._set_data(config["gridsearch_keys"][i], gridsearch_dask_result)
                pbar.update()
            pbar.clear(nolock=False) # prevents double progress at 100%

        gridsearch_meta = {
            'makegrid_key': config["makegrid_key"],
            'grid_points': config["parameters_grid_points"],
            'makegrid_label': config["makegrid_label"],
            'batch_sites': config["batch_sites"],
            'data_key': config["data_key"],
            'gridsearch_key': config["gridsearch_key"],
            'grid_dict': config["grid_dict"],
            'block_length': config["data_block_length"],
            'data_label': config["data_label"],
            'data_source': config["data_source"],
            'gridsearch_keys': config["gridsearch_keys"],
            'data_ndims': config["data_ndims"],
        }
        self._set_meta(config["gridsearch_key"], gridsearch_meta)
        print("[+] Gridsearch can be accessed with %r" % config["gridsearch_key"])

    def _preflight_optimize(self, kwargs):
        '''
        '''
        kwargs['data_key'] = kwargs['data_key']
        if not self._has_key(kwargs['data_key']):
            sys.exit("[X] No data found under key %r." % kwargs['data_key'])
        data_meta = self._get_meta(kwargs['data_key'])
        kwargs['kmax'] = np.array(data_meta['max_k'])
        kwargs['data_source'] = data_meta['data_source']
        #print(dict(data_meta))
        if kwargs['data_source'] == 'meas' and data_meta['data_type'] == 'blocks' and kwargs['windowsum']:
            sys.exit("[X] '--windowsum' not a valid option for data under %r" % kwargs['data_key'])
        kwargs['data_label'] = "%s%s" % (kwargs['data_key'].split("/")[1], "" if not kwargs['windowsum'] else ".windowsum")
        kwargs['optimize_key'] = "optimize/%s/%s" % (kwargs['data_label'], kwargs["optimize_label"])
        if self._has_key(kwargs['optimize_key']):
            if not kwargs['overwrite']: 
                sys.exit(
                    "[X] Analysis with label %r on data %r already exist. Specify '-f' to overwrite."
                    % (kwargs["optimize_label"],  kwargs['data_key']))
            self._del_data_and_meta(kwargs['optimize_key'])
        
        # [ToDo] data should be partitioned into jobs immediately
        # lib.optimize.get_agemo_nlopt_args()
        kwargs["block_length"] = data_meta["block_length"]
        if kwargs["data_source"] == "meas":
            # set FGVs to 0 (needs to be fixed for blocks, not important since FGVs don't come out of agemo unless anomaly)
            # fgv_idxs = get_fgv_idxs(np.array(data_meta['max_k']))
            # tally = self.data[kwargs["data_key"]] if not kwargs['windowsum'] else np.sum(self.data[kwargs["data_key"]], axis=0)
            # tally[fgv_idxs] = 0
            # data = [(0, tally)]
            data = [(0, self.data[kwargs["data_key"]]) if not kwargs['windowsum'] else (0, np.sum(self.data[kwargs["data_key"]], axis=0))]
            kwargs["nlopt_chains"] = data_meta['windows'] if (data_meta['data_ndims'] == 5 and not kwargs['windowsum']) else 1 # blocks/windowsum => 1, windows => n
            kwargs["nlopt_runs"] = 1
        else:
            data = [(replicate_idx, tally) if not kwargs['windowsum'] else (replicate_idx, np.sum(tally, axis=0)) for replicate_idx, tally in self.data[kwargs["data_key"]].arrays()]
            kwargs["nlopt_chains"] = data_meta['windows'] if not kwargs['windowsum'] else 1
            kwargs["nlopt_runs"] = data_meta['replicates'] 
        #print("len(data)", len(data), data[0][1].shape)
        # demography
        def float_or_none(_dict, _parameter):
            return None if _dict[_parameter] is None or len(_dict[_parameter]) == 2 else _dict[_parameter][0]
        gimbleDemographyInstance = GimbleDemographyInstance(
                model=kwargs['model'], 
                mu=kwargs['mu'], 
                ref_pop=kwargs['ref_pop'], 
                block_length=kwargs['block_length'], 
                sync_pops=kwargs['sync_pops'],
                Ne_A=float_or_none(kwargs, 'Ne_A'),
                Ne_B=float_or_none(kwargs, 'Ne_B'),
                Ne_A_B=float_or_none(kwargs, 'Ne_A_B'),
                me=float_or_none(kwargs, 'me'),
                T=float_or_none(kwargs, 'T'), 
                kmax=kwargs['kmax'])
        #print('gimbleDemographyInstance', gimbleDemographyInstance)
        kwargs['gimbleDemographyInstance'] = gimbleDemographyInstance
        # NLOPT parameters
        kwargs['nlopt_parameters'] = []
        kwargs['nlopt_parameters_fixed'] = list(gimbleDemographyInstance.fixed_parameters)
        kwargs['nlopt_parameters_bound'] = []
        kwargs['nlopt_lower_bound'] = []
        kwargs['nlopt_upper_bound'] = []
        kwargs['optimize_result_keys'] = []
        sync_done = False
        for parameter in gimbleDemographyInstance.order_of_parameters:
            if len(kwargs[parameter]) == 2:
                kwargs['nlopt_parameters_bound'].append(parameter)
                if parameter in gimbleDemographyInstance.sync_pops:
                    if not sync_done:
                        kwargs['nlopt_parameters'].append('Ne_s')
                        kwargs['nlopt_lower_bound'].append(min(kwargs[parameter]))
                        kwargs['nlopt_upper_bound'].append(max(kwargs[parameter]))
                        sync_done = True
                else:
                    kwargs['nlopt_parameters'].append(parameter)
                    kwargs['nlopt_lower_bound'].append(min(kwargs[parameter]))
                    kwargs['nlopt_upper_bound'].append(max(kwargs[parameter]))
        startpoints = {
            'midpoint' : np.mean(np.vstack((kwargs['nlopt_lower_bound'], kwargs['nlopt_upper_bound'])), axis=0),
            'random': np.random.uniform(low=kwargs['nlopt_lower_bound'], high=kwargs['nlopt_upper_bound'])
        }
        kwargs['nlopt_start_point'] = startpoints[kwargs['start_point_method']]
        kwargs['nlopt_seed'] = kwargs['seed']
        kwargs['nlopt_log_prefix'] = "%s.%s" % (self.prefix, kwargs['optimize_key'].replace("/", "_"))
        return (data, kwargs)

    def optimize(self, **kwargs):
        data, config = self._preflight_optimize(kwargs)
        optimize_analysis_start_time = timer()  # used for saving elapsed time in meta
        print("[+] Building agemo evaluators ...")
        evaluator_agemo = self.get_agemo_evaluator(
            model=config['model'], 
            kmax=config["kmax"],
            seed=config["seed"])
        # fallback_evaluator gets created if IM
        fallback_evaluator_agemo=self.get_agemo_evaluator(
            model='DIV', 
            kmax=config["kmax"],
            seed=config["seed"]) if config['model'].startswith('IM') else None
        print("[+] Searching parameter space ...")
        for replicate_idx, dataset in data:
            optimize_instance_start_time = timer()
            optimize_result = lib.optimize.agemo_optimize(
                evaluator_agemo, replicate_idx, dataset, config, fallback_evaluator=fallback_evaluator_agemo
            )
    
            optimize_time = format_time(timer() - optimize_instance_start_time)
            optimize_result_key = self.save_optimize_instance(
                replicate_idx, config, optimize_result, optimize_time, True
            )
            config['optimize_result_keys'].append(optimize_result_key)
        if config['data_source'] == 'meas' and optimize_result['dataset_count'] == 1: # save demes for blocks/windowsum only...
            print("[+] Saving demes graph with inferred values ...")
            config['gimbleDemographyInstance'].add_optimize_values(optimize_result['nlopt_values_by_windows_idx'][0])
            config['deme'] = config['gimbleDemographyInstance'].demes_dump()
        config["optimize_time"] = format_time(timer() - optimize_analysis_start_time)
        optimize_meta = {
            'optimize_label': config['optimize_label'], # : 'optimize_cli', 
            'data_key': config['data_key'], # : data_key
            'windowsum': config['windowsum'], # : False, 
            'Ne_A': config['Ne_A'], # : [10000.0, 100000.0], 
            'Ne_B': config['Ne_B'], # : [10000.0, 100000.0], 
            'Ne_A_B': config['Ne_A_B'], # : [20000.0, 100000.0], 
            'T': config['T'], # : [40000.0, 400000.0], 
            'me': config['me'], # : [0.0, 1e-07], 
            'model': config['model'], # : 'IM_BA', 
            'sync_pops': config['sync_pops'], # : ['Ne_B', 'Ne_A'], 
            'ref_pop': config['ref_pop'], # : 'Ne_A', 
            'mu': config['mu'], # : 2e-09, 
            'max_k': list([int(k) for k in config["kmax"]]), # : array([2, 2, 2, 2]), 
            'processes': config['processes'], # : 1, 
            'seed': config['seed'], # : 20, 
            'overwrite': config['overwrite'], # : True, 
            'start_point_method': config['start_point_method'], # : 'midpoint', 
            'nlopt_maxeval': config['nlopt_maxeval'], # : 100, 
            'nlopt_xtol_rel': config['nlopt_xtol_rel'], # : -1.0, 
            'nlopt_ftol_rel': config['nlopt_ftol_rel'], # : -1.0, 
            'nlopt_algorithm': config['nlopt_algorithm'], # : 'CRS2', 
            'data_key': config['data_key'], # : 'tally/blocks', 
            'data_source': config['data_source'], # : 'meas', 
            'data_label': config['data_label'], # : 'blocks', 
            'optimize_key': config['optimize_key'], # : 'optimize/blocks/optimize_cli', 
            'block_length': config['block_length'], # : 64, 
            'nlopt_chains': config['nlopt_chains'], # : 1, 
            'nlopt_runs': config['nlopt_runs'], # : 1, 
            'nlopt_parameters': config['nlopt_parameters'], # : ['Ne_A_B', 'Ne_s', 'me', 'T'], 
            'nlopt_parameters_fixed': config['nlopt_parameters_fixed'], # : [], 
            'nlopt_parameters_bound': config['nlopt_parameters_bound'], # : ['Ne_A_B', 'Ne_A', 'Ne_B', 'me', 'T'], 
            'nlopt_lower_bound': config['nlopt_lower_bound'], # : [20000.0, 10000.0, 0.0, 40000.0], 
            'nlopt_upper_bound': config['nlopt_upper_bound'], # : [100000.0, 100000.0, 1e-07, 400000.0], 
            'optimize_result_keys': config['optimize_result_keys'], # : ['optimize/blocks/optimize_cli/0'], 
            'nlopt_start_point': list([float(k) for k in config["nlopt_start_point"]]), # : array([ 60000.        ,  55000.        ,      0.00000005, 220000.        ]), 
            'optimize_time': config['optimize_time'], # : '00h:00m:42.680s'}
            'deme': config.get('deme', 'N/A'),
        }
        self._set_meta(config["optimize_key"], meta=optimize_meta)
        print("[+] Optimization results saved under label %r" % config["optimize_key"])

    def save_optimize_instance(
        self, data_idx, config, optimize_result, optimize_time, overwrite
    ):
        optimize_result_key = "%s/%s" % (config["optimize_key"], data_idx)
        #print("optimize_key instance", optimize_result_key)
        optimize_meta = {}
        optimize_meta["nlopt_log_iteration_header"] = optimize_result[
            "nlopt_log_iteration_header"
        ]
        optimize_meta["optimize_results"] = []
        for windows_idx in range(optimize_result["dataset_count"]):
            result = optimize_result["nlopt_values_by_windows_idx"][windows_idx]
            result["nlopt_iterations"] = optimize_result["nlopt_evals_by_windows_idx"][windows_idx]
            result["nlopt_exit_code"] = optimize_result["nlopt_status_by_windows_idx"][windows_idx]
            result["lnCL"] = optimize_result["nlopt_optimum_by_windows_idx"][windows_idx]
            result["anomaly_count"] = optimize_result["nlopt_anomalies_by_windows_idx"][windows_idx]
            result["nlopt_time"] = optimize_time
            optimize_meta["optimize_results"].append(result)
        self._set_meta_and_data(
            optimize_result_key, optimize_meta, optimize_result["nlopt_log_iteration_array"]
        )
        return optimize_result_key

    def _preflight_tally(
        self,
        data_type,
        data_label,
        max_k,
        sample_sets,
        sequence_ids,
        genome_file,
        overwrite,
    ):
        config = {
            "data_ndims": 0,
            "data_key": data_type,
            "data_source": data_type,
            "data_type": data_type,
            "data_label": data_label,
            "tally_key": 'tally/%s' % data_label,
            "max_k": max_k,
            "sample_sets": sample_sets,
            "sequences": sequence_ids,
            "genome_file": genome_file,
            "blocks": 0,
            "windows": 0,
            "marginalty": "0.0%",
            "block_length": 0,
        }
        # check data is there
        if not self._has_key(config["data_key"]):
            sys.exit("[X] gimbleStore has no %r." % config["data_key"])
        config["block_length"] = self._get_meta("blocks")[
            "length"
        ]  # needs fixing if multiple block-datasets
        # check tally key
        if self._has_key(config["tally_key"]):
            if not overwrite:
                sys.exit(
                    "[X] Tally with label %r already exist. Specify '-f' to overwrite."
                    % (data_label)
                )
            self._del_data_and_meta(config["tally_key"])
        # sort out sequences (this could be prettier)
        seq_names = self._get_meta("seqs")["seq_names"]
        if config["sequences"] is None:
            config["sequences"] = seq_names
        else:
            sequences_user = set([_ for _ in config["sequences"].split(",") if _])
            config["sequences"] = [
                seq_name for seq_name in seq_names if seq_name in sequences_user
            ]
        return config

    def tally(
        self,
        data_type,
        data_label,
        max_k,
        sample_sets,
        sequence_ids,
        genome_file,
        overwrite,
        verbose=True,
        tally_form="bsfs",
    ):
        config = self._preflight_tally(
            data_type,
            data_label,
            max_k,
            sample_sets,
            sequence_ids,
            genome_file,
            overwrite,
        )
        variation = self._get_variation(
            data_type=config["data_type"],
            sample_sets=config["sample_sets"],
            sequences=config["sequences"],
            progress=verbose,
        )
        config["windows"] = variation.shape[0] if variation.ndim == 3 else 0
        config["blocks"] = (
            (variation.shape[0] * variation.shape[1])
            if variation.ndim == 3
            else variation.shape[0]
        )
        config["marginality"] = format_percentage(
            calculate_marginality_of_variation(variation, max_k=config["max_k"])
        )
        if verbose:
            if config["data_type"] == "blocks":
                print(
                    "[+] Found %s blocks on %s sequence(s)."
                    % (
                        format_count(config["blocks"]),
                        format_count(len(config["sequences"])),
                    )
                )
            else:
                print(
                    "[+] Found %s blocks across %s (sliding) windows (%s blocks per window) on %s sequence(s)."
                    % (
                        format_count(config["blocks"]),
                        format_count(config["windows"]),
                        format_count(variation.shape[1]),
                        format_count(len(config["sequences"]))))
            print(
                "[+] Percentage of blocks treated as marginals (w/ kmax = %s) = %s"
                % (config["max_k"], config["marginality"])
            )
            print("[+] Tally'ing variation data ... ")
        #if config["data_source"] == "windowsum":  # data_source, NOT data_type
        #    variation = variation.reshape(-1, variation.shape[-1])
        variation_tally = tally_variation(variation, form=tally_form, max_k=config["max_k"])
        config["data_ndims"] = variation_tally.ndim
        self.save_tally(config, variation_tally, verbose)
        return config["tally_key"]

    def save_tally(self, config, tally, verbose=True):
        tally_key = config["tally_key"]
        tally_meta = {
            "data_ndims": config["data_ndims"],
            "data_key": config["data_key"],  # where to find the data that went into tally
            "data_source": 'meas',  # meas | sims ?
            "data_type": config["data_type"],
            "tally_key": config["tally_key"],  # where to find the tally
            "max_k": (None if config["max_k"] is None else tuple([int(v) for v in config["max_k"]])),
            "sample_sets": config.get("sample_sets", "NA"),
            "sequences": config.get("sequences", []),
            "genome_file": config.get("genome_file", None),
            "blocks": config.get("blocks", 0),
            "windows": config.get("windows", 0),
            "marginality": config.get("marginality", "NA"),
            "block_length": config["block_length"],
        }
        self._set_meta_and_data(tally_key, tally_meta, tally)

    def _has_key(self, key):
        return (key in self.data) if key else False

    def _set_data(self, key, array, overwrite=True, compression=False):
        #print("array", array.shape, array.nbytes)
        compressor = numcodecs.Blosc(cname='zstd', clevel=5, shuffle=numcodecs.Blosc.BITSHUFFLE) if compression else None
        self.data.create_dataset(key, data=array, overwrite=overwrite, compressor=compressor)
        #print(self.data[key].info)

    def _set_meta_and_data(self, key, meta, array, overwrite=True):
        self.data.create_dataset(key, data=array, overwrite=overwrite)
        self.data[key].attrs.put(meta)

    def delete_key(self, data_key):
        if self._has_key(data_key):
            size = format_bytes(recursive_get_size("%s/%s" % (self.path, data_key)))
            while True:
                sys.stdout.write("[?] Delete %r of size %s ? [y/N] " % (data_key, size))
                choice = input().lower()
                if choice in set(["", "n", "no"]):
                    print("[+] Deletion cancelled.")
                    return 1
                elif choice in set(["y", "yes", "ye"]):
                    del self.data[data_key]
                    print("[+] %r deleted." % data_key)
                    return 1
                else:
                    sys.stdout.write("Please respond with 'y' or 'n'\n")            
        print("[+] %r not found." % data_key)

    def _del_data_and_meta(self, key):
        if self._has_key(key):
            del self.data[key]

    def _get_data(self, key, dtype=None):
        if self._has_key(key):
            data = self.data[key]
            if dtype is None:
                return self.data[key]
            else:
                if np.can_cast(data.dtype, dtype):
                    return self.data[key]
                else:
                    raise ValueError("Can't cast %r to %r" % (data.dtype, dtype))
        return None

    def _get_meta(self, key):
        if self._has_key(key):
            return self.data[key].attrs
        return None

    def _set_meta(self, key, meta={}):
        self.data.require_group(key)
        self.data[key].attrs.put(meta)

    def makegrid(self, Ne_A, Ne_B, Ne_A_B, T, me, makegrid_label, model, block_length, ref_pop, mu, kmax, processes, seed, overwrite):
        makegrid_key = "makegrid/%s" % (makegrid_label)
        if self._has_key(makegrid_key):
            if not overwrite:
                sys.exit("[X] Grid with label %r already exist. Specify '-f' to overwrite." % makegrid_label)
            self._del_data_and_meta(makegrid_key)
        parameter_dicts = get_parameter_dicts_from_user_parameters(Ne_A=Ne_A, Ne_B=Ne_B, Ne_A_B=Ne_A_B, T=T, me=me)
        def get_agemo_parameters(parameter_dicts, model, block_length, ref_pop, mu, kmax):
            agemo_parameters = []
            for parameter_dict in parameter_dicts:
                model_instance = GimbleDemographyInstance(
                    model=model, 
                    Ne_A=parameter_dict.get('Ne_A', None), 
                    Ne_B=parameter_dict.get('Ne_B', None), 
                    Ne_A_B=parameter_dict.get('Ne_A_B', None), 
                    me=parameter_dict.get('me', None), 
                    T=parameter_dict.get('T', None),
                    mu=mu, 
                    ref_pop=ref_pop, 
                    block_length=block_length, 
                    kmax=kmax) 
                agemo_parameters.append(model_instance.get_agemo_values(fallback=True))
            return agemo_parameters
        agemo_parameters = get_agemo_parameters(parameter_dicts, model, block_length, ref_pop, mu, kmax)
        print(
            "[+] Grid of %s parameter-grid-points will be prepared..."
            % len(agemo_parameters)
        )
        print("[+] Using 'agemo v%s'" % (AGEMO_VERSION))
        try:
            evaluator_agemo = self.get_agemo_evaluator(
                model=model, 
                kmax=kmax,
                seed=seed)
            # fallback_evaluator gets created if IM
            fallback_evaluator_agemo=self.get_agemo_evaluator(
                model='DIV', 
                kmax=kmax,
                seed=seed) if model.startswith('IM') else None
        except TypeError:
            sys.exit("[X] Agemo version is outdated. Please upgrade to 'agemo v%s'\n\tpip install agemo -U" % MIN_AGEMO_VERSION)
        grid = self.evaluate_grid(evaluator_agemo, agemo_parameters, parameter_dicts, processes=processes, fallback_evaluator=fallback_evaluator_agemo)
        
        # ToDo:
        # remove grid_dict: is still required in gridsearch and query ... needs to be removed there
        # rename parameters_grid_points to grid_points
        grid_meta = {
            "makegrid_key": makegrid_key,
            "makegrid_label": makegrid_label,
            "block_length": block_length,
            "mu": mu,
            "model": model,
            "reference_pop_id": ref_pop,
            "max_k": list([int(k) for k in kmax]),
            "Ne_A": Ne_A,
            "Ne_B": Ne_B,
            "Ne_A_B": Ne_A_B,
            "T": T,
            "me": me,
            "processes": processes,
            "seed": seed,
            "overwrite": overwrite,
            "grid_dict": {k: list(v) for k,v in LOD_to_DOL(get_parameter_dicts_from_user_parameters(Ne_A=Ne_A, Ne_B=Ne_B, Ne_A_B=Ne_A_B, T=T, me=me)).items()},
            "parameters_grid_points": len(agemo_parameters),
        } 

        #print(grid_meta)
        self._set_meta_and_data(makegrid_key, grid_meta, grid)
        print(
            "[+] Grid can be accessed with '--grid_key %s'"
            % grid_meta["makegrid_key"]
        )

    def save_grid(self, config, grid):
        grid_meta = config_to_meta(config, "makegrid")
        self._set_meta_and_data(config["key"], grid_meta, grid)
        print(
            "[+] Grid can be accessed with '--grid_label %s'"
            % grid_meta["makegrid_key"]
        )

    def evaluate_grid(self, evaluator_agemo, agemo_parameters, parameter_dicts, processes=1, fallback_evaluator=None, verbose=False, agemo_anomaly_tol=1e-5):
        if verbose:
            for parameter in agemo_parameters:
                print(float(parameter[0]), [float(x) for x in parameter[1]], float(parameter[2]), parameter[3])
        all_ETPs = []
        print("[+] Calculating probabilities of mutation configurations under %s parameter combinations" % len(agemo_parameters))
        if processes==1:
            for agemo_parameter in tqdm(agemo_parameters, desc="[%] Progress", ncols=100):
                theta_branch, var, time, fallback_flag = agemo_parameter
                evaluator = evaluator_agemo if not fallback_flag else fallback_evaluator
                result = evaluator.evaluate(theta_branch, var, time=time)
                all_ETPs.append(result)
        else:
            global EVALUATOR
            global FALLBACK_EVALUATOR
            EVALUATOR = evaluator_agemo
            FALLBACK_EVALUATOR = fallback_evaluator
            with multiprocessing.Pool(processes=processes) as pool:
                for ETP in pool.starmap(self.multi_eval, tqdm(agemo_parameters, ncols=100, desc="[%] Progress")):
                    all_ETPs.append(ETP)
        print("[+] Checking probabilities for anomalies ...")
        anomaly_count = 0
        for idx, etps in enumerate(all_ETPs):
            if (not np.isclose(np.sum(etps), 1, rtol=agemo_anomaly_tol)):
                anomaly_count += 1
                print("[-] Anomaly found: %s" % (", ".join(["%s=%s" % (k, v) for k, v in parameter_dicts[idx].items()])))
        if anomaly_count:
            print("[X] %s anomalies detected. Please consider changing the parameter space as this grid might behave strangely." % anomaly_count)
        else:
            print("[+] No anomalies found. All is well.")
        return np.array(all_ETPs, dtype=np.float64)

    def multi_eval(self, theta_branch, var, time, fallback_flag):
        #theta_branch, var, time, fallback_flag = agemo_parameter
        evaluator = EVALUATOR if not fallback_flag else FALLBACK_EVALUATOR
        return evaluator.evaluate(theta_branch, var, time=time) 

    def get_agemo_evaluator(self, model=None, kmax=None, seed=None):
        mutation_shape = tuple(kmax + 2)
        if model == "DIV":
            sample_configuration = [(), ('a', 'a'), ('b', 'b')]
            events = [agemo.PopulationSplitEvent(len(sample_configuration), 0, 1, 2)]
        if model == "MIG_AB":
            #sample_configuration = [(), ('a', 'a'), ('b', 'b')]
            sample_configuration = [('a', 'a'), ('b', 'b')]
            #events = [agemo.MigrationEvent(len(sample_configuration), 1, 2)]
            events = [agemo.MigrationEvent(len(sample_configuration), 0, 1)]
        if model == "MIG_BA":
            #sample_configuration = [(), ('a', 'a'), ('b', 'b')]
            sample_configuration = [('a', 'a'), ('b', 'b')]
            #events = [agemo.MigrationEvent(len(sample_configuration), 2, 1)]
            events = [agemo.MigrationEvent(len(sample_configuration), 1, 0)]
        if model == "IM_AB":
            sample_configuration = [(), ('a', 'a'), ('b', 'b')]
            events = [
                agemo.MigrationEvent(len(sample_configuration), 1, 2),
                agemo.PopulationSplitEvent(len(sample_configuration) + 1, 0, 1, 2)]
        if model == "IM_BA":
            sample_configuration = [(), ('a', 'a'), ('b', 'b')]
            events = [
                agemo.MigrationEvent(len(sample_configuration), 2, 1),
                agemo.PopulationSplitEvent(len(sample_configuration) + 1, 0, 1, 2)]
        # Default order of branches is different and therefor needs to be adjusted:
        # - agemo : [hetA, hetB, fixed, hetAB]
        #       branchtype_dict = {'a': 0, 'abb': 0, 'b': 1, 'aab': 1, 'aa': 2, 'bb': 2, 'ab': 3} 
        # - gimble : [hetB, hetA, hetAB, fixed]
        #       branchtype_dict = {'a': 1, 'abb': 1, 'b': 0, 'aab': 0, 'aa': 3, 'bb': 3, 'ab': 2}
        branchtype_dict = {'a': 1, 'abb': 1, 'b': 0, 'aab': 0, 'aa': 3, 'bb': 3, 'ab': 2}
        branch_type_object = agemo.BranchTypeCounter(sample_configuration, branchtype_dict=branchtype_dict)
        #num_branchtypes = len(branch_type_object)
        mutation_type_object = agemo.MutationTypeCounter(branch_type_object, mutation_shape)
        gf = agemo.GfMatrixObject(branch_type_object, events)
        evaluator = agemo.BSFSEvaluator(gf, mutation_type_object, seed=seed)
        return evaluator

    def _validate_seq_names(self, sequences=None):
        """Returns valid seq_names in sequences or exits."""
        meta = self._get_meta("seqs")
        if sequences is None:
            return meta["seq_names"]
        if set(sequences).issubset(set(meta["seq_names"])):
            return sequences
        sys.exit(
            "[X] Sequence(s) %r not a subset of sequence(s) %r in ZARR store"
            % (", ".join(sequences), ", ".join(meta["seq_names"]))
        )

    def _get_sample_set_idxs(self, query="X"):
        """Returns list of sample_set_idxs.

        Parameters
        ----------
        query : string or None
                'X' - inter-population sample_sets
                'A' - intra-population sample_sets of population A
                'B' - intra-population sample_sets of population B
                None - all sample_sets

        Returns
        -------
        out : list of strings
            sample_set_idxs that can be used to access data in gimble store
        """
        meta = self._get_meta("seqs")
        if query is None:
            return [str(idx) for idx in range(len(meta["sample_sets"]))]
        elif query == "X":
            return [
                str(idx)
                for (idx, is_cartesian) in enumerate(meta["sample_sets_inter"])
                if is_cartesian
            ]
        elif query == "A":
            return [
                str(idx)
                for (idx, is_intra_A) in enumerate(meta["sample_sets_intra_A"])
                if is_intra_A
            ]
        elif query == "B":
            return [
                str(idx)
                for (idx, is_intra_B) in enumerate(meta["sample_sets_intra_B"])
                if is_intra_B
            ]
        else:
            raise ValueError("'query' must be 'X', 'A', 'B', or None")

    def _get_variation(
        self,
        data_type=None,
        sample_sets="X",
        sequences=None,
        population_by_letter=None,
        progress=False,
    ):
        """Returns variation array of 2 (blocks) or 3 (windows) dimensions.

        Parameters
        ----------
        data_type : 'blocks' or 'windows'
        sample_sets : only needed for data_type 'blocks'. String or None
                None - all sample_sets
                'X' - inter-population sample_sets (default)
                'A' - intra-population sample_sets of population A
                'B' - intra-population sample_sets of population B
            If supplied, array is based only on variation in those sample_sets
        sequences : list of strings or None
            If supplied, array is based only on those sequences
        population_by_letter : dict (string -> string) or None
            Mapping of population IDs to population letter in model (from INI file)

        Returns
        -------
        out : ndarray, int, ndim (mutypes)
        """
        meta = self._get_meta("seqs")
        sequences = self._validate_seq_names(sequences)
        if population_by_letter:
            assert set(population_by_letter.values()) == set(
                meta["population_by_letter"].values()
            ), (
                "population_by_letter %r does not equal populations in ZARR store (%r)"
                % (population_by_letter, meta["population_by_letter"])
            )
        keys = []
        if data_type == "blocks":
            sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
            keys = [
                "blocks/%s/%s/variation" % (seq_name, sample_set_idx)
                for seq_name, sample_set_idx in list(
                    itertools.product(sequences, sample_set_idxs)
                )
            ]
            # print('variation keys', keys)
        elif data_type == "windows":
            keys = ["windows/%s/variation" % (seq_name) for seq_name in sequences]
        else:
            raise ValueError("Invalid datatype: %s" % data_type)
        variations = []
        for key in tqdm(
            keys,
            total=len(keys),
            desc="[%] Preparing data",
            ncols=100,
            unit_scale=True,
            disable=(not progress),
        ):
            # variations.append(np.array(self.data[key], dtype=np.int64))
            if self._has_key(key):
                variations.append(self.data[key])
        if not variations:
            sys.exit("[X] Not enough blocks in Gimble datastore.")
        variation = np.concatenate(variations, axis=0)
        polarise_true = (
            (
                (population_by_letter["A"] == meta["population_by_letter"]["B"])
                and (population_by_letter["B"] == meta["population_by_letter"]["A"])
            )
            if population_by_letter
            else False
        )
        if polarise_true:
            variation[..., [0, 1]] = variation[..., [1, 0]]
        return variation
      
    def _init_store(self, create, overwrite):
        if create:
            if os.path.isdir(self.path):
                print("[-] Gimble datastore %r already exists." % self.path)
                if not overwrite:
                    print("[X] Please specify '-f' to overwrite.")
                    sys.exit(1)
                print("[+] Deleting existing Gimble datastore %r" % self.path)
                shutil.rmtree(self.path)
            print("[+] Creating Gimble datastore in %r" % self.path)
            return zarr.open(str(self.path), mode="w")
        # print("[+] Loading Gimble store from %r" % self.path)
        return zarr.open(str(self.path), mode="r+")

    def _get_window_bed(self):
        meta_seqs = self._get_meta("seqs")
        meta_windows = self._get_meta("windows")
        MAX_SEQNAME_LENGTH = max(
            [len(seq_name) + 1 for seq_name in meta_seqs["seq_names"]]
        )
        window_count = meta_windows.get(
            "window_count", meta_windows.get("count", meta_windows.get("windows", 0))
        )
        index = np.arange(window_count)
        sequences = np.zeros(window_count, dtype="<U%s" % MAX_SEQNAME_LENGTH)
        starts = np.zeros(window_count, dtype=np.int64)
        ends = np.zeros(window_count, dtype=np.int64)
        pos_mean = np.zeros(window_count, dtype=np.float64)
        pos_median = np.zeros(window_count, dtype=np.float64)
        balance = np.zeros(window_count, dtype=np.float64)
        mse_sample_set_cov = np.zeros(window_count, dtype=np.float64)
        offset = 0
        for seq_name in tqdm(
            meta_seqs["seq_names"],
            total=len(meta_seqs["seq_names"]),
            desc="[%] Preparing data",
            ncols=100,
            unit_scale=True,
        ):
            start_key = "windows/%s/starts" % (seq_name)
            end_key = "windows/%s/ends" % (seq_name)
            pos_mean_key = "windows/%s/pos_mean" % (seq_name)
            pos_median_key = "windows/%s/pos_median" % (seq_name)
            balance_key = "windows/%s/balance" % (seq_name)
            mse_sample_set_cov_key = "windows/%s/mse_sample_set_cov" % (seq_name)
            if start_key in self.data:
                start_array = np.array(self.data[start_key])
                _window_count = start_array.shape[0]
                starts[offset : offset + _window_count] = start_array
                ends[offset : offset + _window_count] = np.array(self.data[end_key])
                pos_mean[offset : offset + _window_count] = np.array(
                    self.data[pos_mean_key]
                )
                pos_median[offset : offset + _window_count] = np.array(
                    self.data[pos_median_key]
                )
                balance[offset : offset + _window_count] = np.array(
                    self.data[balance_key]
                )
                mse_sample_set_cov[offset : offset + _window_count] = np.array(
                    self.data[mse_sample_set_cov_key]
                )
                sequences[offset : offset + _window_count] = np.full_like(
                    _window_count, seq_name, dtype="<U%s" % MAX_SEQNAME_LENGTH
                )
                offset += _window_count
        return (
            sequences,
            starts,
            ends,
            index,
            pos_mean,
            pos_median,
            balance,
            mse_sample_set_cov,
        )

    def _preflight_windows(self, window_size, window_step, overwrite=False):
        config = {
            "window_size": window_size,
            "window_step": window_step,
            "sample_sets": "X",
        }
        if not self.has_stage("blocks"):
            sys.exit(
                "[X] Gimble store %r has no blocks. Please run 'gimble blocks'." % self.path
            )
        if self.has_stage("windows"):
            if not overwrite:
                sys.exit(
                    "[X] Gimble store %r already contains windows.\n[X] These windows => %r\n[X] Please specify '--force' to overwrite."
                    % (self.path, self.get_stage("windows"))
                )
            print(
                "[-] Gimble store %r already contains windows. But these will be overwritten..."
                % (self.path)
            )
            """only deletes windows (NOT windows tallies)"""
            self._del_data_and_meta("windows")
        config["windows_key"] = "windows"  # one could allow for multiple windows sets by using labels here
        config["window_count"] = 0
        return config

    def _get_interval_coordinates(self, seq_name=None, sample_set=None):
        if seq_name is None:
            raise ValueError("_get_interval_coordinates: needs seq_name")
        matrix_key = "seqs/%s/intervals/matrix" % seq_name
        if matrix_key in self.data:
            start_key = "seqs/%s/intervals/starts" % seq_name
            end_key = "seqs/%s/intervals/ends" % seq_name
            if sample_set is None:
                return (np.array(self.data[start_key]), np.array(self.data[end_key]))
            meta_seqs = self._get_meta("seqs")
            try:
                sample_set_key = np.array(
                    [
                        meta_seqs["intervals_idx_by_sample"][sample]
                        for sample in sample_set
                    ]
                )
            except KeyError:
                sys.exit(
                    "_get_interval_coordinates: sample_set %s not found in store. Existing samples: %s"
                    % (sample_set, list(meta_seqs["intervals_idx_by_sample"].keys()))
                )
            mask = np.all(np.array(self.data[matrix_key])[:, sample_set_key], axis=1)
            return (
                np.array(self.data[start_key])[mask],
                np.array(self.data[end_key])[mask],
            )
        return (None, None)

    def _get_variants(self, seq_name):
        pos_key = "seqs/%s/variants/pos" % (seq_name)
        gt_key = "seqs/%s/variants/matrix" % (seq_name)
        pos = (
            np.array(self.data[pos_key], dtype=np.int64)
            if pos_key in self.data
            else None
        )
        gt_matrix = (
            allel.GenotypeArray(self.data[gt_key].view(read_only=True))
            if gt_key in self.data
            else None
        )
        if pos is not None and gt_matrix is not None:
            assert (
                pos.shape[0] == gt_matrix.shape[0]
            )  # check whether they are the same length ...
        return (pos, gt_matrix)

    def _make_blocks(self, config):
        meta_seqs = self._get_meta("seqs")
        with tqdm(
            total=(len(meta_seqs["seq_names"]) * len(meta_seqs["sample_sets"])),
            desc="[%] Making pair-blocks",
            ncols=100,
            unit_scale=True,
        ) as pbar:
            for seq_name in meta_seqs["seq_names"]:
                pos, gt_matrix = self._get_variants(seq_name)  # arrays or None
                for sample_set_idx, sample_set in enumerate(meta_seqs["sample_sets"]):
                    # get BED starts/ends from store
                    intervals = self._get_interval_coordinates(
                        seq_name=seq_name, sample_set=sample_set
                    )
                    # turn BED starts/ends into sites-array
                    sites = intervals_to_sites(intervals)
                    # turn sites-array into 2D np.array with block sites (or None)
                    block_sites = sites_to_blocks(
                        sites, config["block_length"], config["block_span"], sample_set
                    )
                    if block_sites is not None:
                        # subset gts of sample_set from gt_matrix (or None)
                        gts = subset_gt_matrix(
                            meta_seqs,
                            sample_set,
                            np.isin(pos, block_sites, assume_unique=True),
                            gt_matrix,
                        )
                        # get block arrays
                        (
                            starts,
                            ends,
                            multiallelic,
                            missing,
                            monomorphic,
                            variation,
                        ) = blocks_to_arrays(block_sites, gts, pos)
                        # save block arrays
                        blocks_raw, blocks_valid = self._set_blocks(
                            seq_name,
                            sample_set_idx,
                            starts,
                            ends,
                            multiallelic,
                            missing,
                            monomorphic,
                            variation,
                            config["block_max_missing"],
                            config["block_max_multiallelic"],
                        )
                        # record counts
                        config["blocks_raw_by_sample_set_idx"][
                            sample_set_idx
                        ] += blocks_raw
                        config["blocks_by_sample_set_idx"][
                            sample_set_idx
                        ] += blocks_valid
                        config["blocks_by_sequence"][seq_name] += blocks_valid
                    pbar.update(1)
        config["count_total"] = sum(
            [count for count in config["blocks_by_sample_set_idx"].values()]
        )
        config["count_total_raw"] = sum(
            [count for count in config["blocks_raw_by_sample_set_idx"].values()]
        )
        return config

    def _set_blocks(
        self,
        seq_name,
        sample_set_idx,
        starts,
        ends,
        multiallelic,
        missing,
        monomorphic,
        variation,
        block_max_missing,
        block_max_multiallelic,
    ):
        valid = (
            np.less_equal(missing, block_max_missing)
            & np.less_equal(multiallelic, block_max_multiallelic)
        ).flatten()
        blocks_starts_key = "blocks/%s/%s/starts" % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_starts_key, data=starts[valid], overwrite=True)
        blocks_ends_key = "blocks/%s/%s/ends" % (seq_name, sample_set_idx)
        self.data.create_dataset(blocks_ends_key, data=ends[valid], overwrite=True)
        blocks_variation_key = "blocks/%s/%s/variation" % (seq_name, sample_set_idx)
        self.data.create_dataset(
            blocks_variation_key, data=variation[valid], overwrite=True
        )
        blocks_missing_key = "blocks/%s/%s/missing" % (seq_name, sample_set_idx)
        self.data.create_dataset(
            blocks_missing_key, data=missing[valid], overwrite=True
        )
        blocks_multiallelic_key = "blocks/%s/%s/multiallelic" % (
            seq_name,
            sample_set_idx,
        )
        self.data.create_dataset(
            blocks_multiallelic_key, data=multiallelic[valid], overwrite=True
        )
        return (valid.shape[0], valid[valid == True].shape[0])

    # def _set_blocks_meta(self, block_length, block_span, block_max_missing, block_max_multiallelic,
    #        blocks_raw_by_sample_set_idx, blocks_by_sample_set_idx, blocks_by_sequence):
    #    meta_blocks = self._get_meta('blocks')
    #    print('meta_blocks', type(meta_blocks), dict(meta_blocks))
    #    if meta_blocks is None:
    #        sys.exit("[X] No blocks could be generated from data given the parameters.")
    #    meta_blocks['length'] = block_length
    #    meta_blocks['span'] = block_span
    #    meta_blocks['max_missing'] = block_max_missing
    #    meta_blocks['max_multiallelic'] = block_max_multiallelic
    #    meta_blocks['count_by_sample_set_idx'] = dict(blocks_by_sample_set_idx) # keys are strings
    #    meta_blocks['count_by_sequence'] = dict(blocks_by_sequence) # keys are strings
    #    meta_blocks['count_raw_by_sample_set_idx'] = dict(blocks_raw_by_sample_set_idx) # keys are strings
    #    meta_blocks['count_total'] = sum([count for count in blocks_by_sample_set_idx.values()])
    #    meta_blocks['count_total_raw'] = sum([count for count in blocks_raw_by_sample_set_idx.values()])
    #    print('meta_blocks', type(meta_blocks), dict(meta_blocks))

    def _get_block_coordinates(self, sample_sets=None, sequences=None):
        sequences = self._validate_seq_names(sequences)
        sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
        keys_start_end = [
            (
                "blocks/%s/%s/starts" % (seq_name, sample_set_idx),
                "blocks/%s/%s/ends" % (seq_name, sample_set_idx),
            )
            for seq_name, sample_set_idx in list(
                itertools.product(sequences, sample_set_idxs)
            )
        ]
        block_starts, block_ends = [], []
        for start_key, end_key in keys_start_end:
            block_starts.append(np.array(self.data[start_key], dtype=np.int64))
            block_ends.append(np.array(self.data[end_key], dtype=np.int64))
        block_start = np.concatenate(block_starts, axis=0)
        block_end = np.concatenate(block_ends, axis=0)
        return (block_start, block_end)

    def _get_block_sample_set_idxs(self, sample_sets=None, sequences=None):
        sequences = self._validate_seq_names(sequences)
        sample_set_idxs = self._get_sample_set_idxs(query=sample_sets)
        key_by_sample_set_idx = {
            sample_set_idx: "blocks/%s/%s/starts" % (seq_name, sample_set_idx)
            for seq_name, sample_set_idx in list(
                itertools.product(sequences, sample_set_idxs)
            )
        }
        block_sample_set_idxs = []
        for sample_set_idx, key in key_by_sample_set_idx.items():
            block_sample_set_idxs.append(
                np.full(self.data[key].shape[0], int(sample_set_idx))
            )
        block_sample_set_idxs = np.concatenate(block_sample_set_idxs, axis=0)
        return block_sample_set_idxs

    def _make_windows(self, config):
        ### meta_seqs['seq_names'] => order
        meta_blocks = self._get_meta("blocks")
        if not meta_blocks:
            sys.exit("[X] No blocks found.")
        meta_seqs = self._get_meta("seqs")
        sample_set_idxs = np.array(
            self._get_sample_set_idxs(query=config["sample_sets"]), dtype=np.int64
        )
        blockable_seqs, unblockable_seqs = [], []
        for seq_name in meta_seqs["seq_names"]:
            # problem here is that meta_blocks["count_by_sequence"] contains ALL blocks ... 
            # should contain only blocks from inter-pop sample-sets
            if meta_blocks["count_by_sequence"].get(seq_name, 0) >= config["window_size"]:
                blockable_seqs.append(seq_name)
            else:
                unblockable_seqs.append(seq_name)
        if not blockable_seqs:
            sys.exit(
                "[X] Not enough blocks to make windows of this size (%s)."
                % (config["window_size"])
            )
        print(
            "[+] Making windows along %s sequences (%s sequences excluded)"
            % (len(blockable_seqs), len(unblockable_seqs))
        )
        for seq_name in tqdm(
            blockable_seqs, total=len(blockable_seqs), desc="[%] Progress", ncols=100,
        ):
            block_sample_set_idxs = self._get_block_sample_set_idxs(
                sample_sets=config["sample_sets"], sequences=[seq_name]
            )
            if block_sample_set_idxs.size >= config["window_size"]:

                block_variation = self._get_variation(
                    data_type="blocks",
                    sample_sets=config["sample_sets"],
                    sequences=[seq_name],
                )
                block_starts, block_ends = self._get_block_coordinates(
                    sample_sets=config["sample_sets"], sequences=[seq_name]
                )
                windows = blocks_to_windows(
                    sample_set_idxs,
                    block_variation,
                    block_starts,
                    block_ends,
                    block_sample_set_idxs,
                    config["window_size"],
                    config["window_step"],
                )
                (
                    window_variation,
                    window_starts,
                    window_ends,
                    window_pos_mean,
                    window_pos_median,
                    balance,
                    mse_sample_set_cov,
                ) = windows
                config["window_count"] += self._set_windows(
                    seq_name,
                    window_variation,
                    window_starts,
                    window_ends,
                    window_pos_mean,
                    window_pos_median,
                    balance,
                    mse_sample_set_cov,
                )
        print(
            "[+] Made %s window(s)." % format_count(config["window_count"])
        )
        return config

    def _set_windows(
        self,
        seq_name,
        window_variation,
        window_starts,
        window_ends,
        window_pos_mean,
        window_pos_median,
        balance,
        mse_sample_set_cov,
    ):
        self.data.create_dataset(
            "windows/%s/variation" % seq_name, data=window_variation, overwrite=True
        )
        self.data.create_dataset(
            "windows/%s/starts" % seq_name, data=window_starts, overwrite=True
        )
        self.data.create_dataset(
            "windows/%s/ends" % seq_name, data=window_ends, overwrite=True
        )
        self.data.create_dataset(
            "windows/%s/pos_mean" % seq_name, data=window_pos_mean, overwrite=True
        )
        self.data.create_dataset(
            "windows/%s/pos_median" % seq_name, data=window_pos_median, overwrite=True
        )
        self.data.create_dataset(
            "windows/%s/balance" % seq_name, data=balance, overwrite=True
        )
        self.data.create_dataset(
            "windows/%s/mse_sample_set_cov" % seq_name,
            data=mse_sample_set_cov,
            overwrite=True,
        )
        return window_variation.shape[0]

    ####################### REPORTS ######################

    def _get_parse_report(self, width):
        meta_seqs = self._get_meta("seqs")
        reportObj = ReportObj(width=width)
        reportObj.add_line(
            prefix="[+]", left="[", center="Parsed data", right="]", fill="="
        )
        reportObj.add_line(prefix="[+]", left="seqs")
        right = "%s in %s sequence(s) (n50 = %s)" % (
            format_bases(sum(meta_seqs["seq_lengths"])),
            format_count(len(meta_seqs["seq_lengths"])),
            format_bases(meta_seqs["seq_n50"]),
        )
        reportObj.add_line(
            prefix="[+]", branch="T", fill=".", left="genome", right=right
        )
        right = "%s samples in %s populations" % (
            format_count(len(meta_seqs["samples"])),
            format_count(len(meta_seqs["population_ids"])),
        )
        reportObj.add_line(
            prefix="[+]", branch="T", left="populations", fill=".", right=right
        )
        sample_counts_by_population = collections.Counter(meta_seqs["populations"])
        for idx, (letter, population_id) in enumerate(
            meta_seqs["population_by_letter"].items()
        ):
            left = "%s = %s" % (letter, population_id)
            right = "%s" % format_count(sample_counts_by_population[population_id])
            branch = "P%s" % (
                "F" if idx == len(meta_seqs["population_by_letter"]) - 1 else "T"
            )
            reportObj.add_line(prefix="[+]", branch=branch, left=left, right=right)
        reportObj.add_line(
            prefix="[+]",
            branch="T",
            left="sample sets",
            fill=".",
            right=format_count(len(meta_seqs["sample_sets"])),
        )
        reportObj.add_line(
            prefix="[+]",
            branch="PT",
            left="INTER-population sample-sets (X)",
            right=format_count(len(self._get_sample_set_idxs(query="X"))),
        )
        reportObj.add_line(
            prefix="[+]",
            branch="PT",
            left="INTRA-population sample-sets (A)",
            right=format_count(len(self._get_sample_set_idxs(query="A"))),
        )
        reportObj.add_line(
            prefix="[+]",
            branch="PF",
            left="INTRA-population sample-sets (B)",
            right=format_count(len(self._get_sample_set_idxs(query="B"))),
        )
        reportObj.add_line(
            prefix="[+]",
            branch="T",
            left="variants",
            fill=".",
            right=(
                "%s (%s per 1 kb)"
                % (
                    format_count(meta_seqs["variants_counts"]),
                    format_proportion(
                        1000
                        * meta_seqs["variants_counts"]
                        / sum(meta_seqs["seq_lengths"])
                    ),
                )
            ),
        )
        reportObj.add_line(
            prefix="[+]",
            branch="PP",
            right="".join([c.rjust(8) for c in ["HOMREF", "HOMALT", "HET", "MISS"]]),
        )
        for idx, sample in enumerate(meta_seqs["samples"]):
            variant_idx = meta_seqs["variants_idx_by_sample"][sample]
            branch = (
                "PF" if idx == len(meta_seqs["variants_idx_by_sample"]) - 1 else "PT"
            )
            left = sample
            right = "%s %s %s %s" % (
                (
                    format_percentage(
                        meta_seqs["variants_counts_hom_ref"][variant_idx]
                        / meta_seqs["variants_counts"]
                    )
                    if meta_seqs["variants_counts"]
                    else format_percentage(0)
                ).rjust(7),
                (
                    format_percentage(
                        meta_seqs["variants_counts_hom_alt"][variant_idx]
                        / meta_seqs["variants_counts"]
                    )
                    if meta_seqs["variants_counts"]
                    else format_percentage(0)
                ).rjust(7),
                (
                    format_percentage(
                        meta_seqs["variants_counts_het"][variant_idx]
                        / meta_seqs["variants_counts"]
                    )
                    if meta_seqs["variants_counts"]
                    else format_percentage(0)
                ).rjust(7),
                (
                    format_percentage(
                        meta_seqs["variants_counts_missing"][variant_idx]
                        / meta_seqs["variants_counts"]
                    )
                    if meta_seqs["variants_counts"]
                    else format_percentage(0)
                ).rjust(7),
            )
            reportObj.add_line(prefix="[+]", branch=branch, left=left, right=right)
        reportObj.add_line(
            prefix="[+]",
            branch="F",
            left="intervals",
            fill=".",
            right="%s intervals across %s (%s of genome)"
            % (
                format_count(meta_seqs["intervals_count"]),
                format_bases(meta_seqs["intervals_span"]),
                format_percentage(
                    meta_seqs["intervals_span"] / sum(meta_seqs["seq_lengths"])
                ),
            ),
        )
        return reportObj

    def _get_storage_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(
            prefix="[+]", left="[", center="Storage", right="]", fill="="
        )
        total_size = recursive_get_size(self.path)
        total_right = "%s | %s" % (
            format_bytes(total_size),
            format_percentage(total_size / total_size),
        )
        reportObj.add_line(
            prefix="[+]", left=pathlib.Path(self.path).name, right=total_right, fill="."
        )
        for idx, group in enumerate(self.data):
            size = recursive_get_size(pathlib.Path(self.path) / pathlib.Path(group))
            percentage = format_percentage(size / total_size)
            branch = branch = "T" if idx < len(self.data) - 1 else "F"
            reportObj.add_line(
                prefix="[+]",
                branch=branch,
                left=group,
                right="%s | %s" % (format_bytes(size), percentage.rjust(7)),
                fill=".",
            )
        return reportObj

    def _get_blocks_report_metrics(self):
        meta_blocks = self._get_meta("blocks")
        block_length = meta_blocks["length"]
        meta_seqs = self._get_meta("seqs")
        intervals_span = meta_seqs["intervals_span"]
        BRMs = {}
        for sample_sets in ["X", "A", "B"]:
            sample_sets_count = len(self._get_sample_set_idxs(sample_sets))
            tally = tally_variation(
                self._get_variation(data_type="blocks", sample_sets=sample_sets),
                form="tally",
            )
            BRM = calculate_blocks_report_metrics(
                tally, sample_sets_count, block_length, intervals_span
            )
            for k, v in BRM.items():
                if k in ["interval_coverage", "blocks_invariant", "blocks_fgv"]:
                    BRM[k] = format_percentage(v, precision=2)
                elif k in ["blocks_total"]:
                    BRM[k] = format_count(v)
                else:
                    BRM[k] = format_proportion(v, precision=5)
            if sample_sets == "X":
                del BRM["pi"]
                del BRM["watterson_theta"]
            if sample_sets == "A" or sample_sets == "B":
                if sample_sets == "A":
                    del BRM["heterozygosity_B"]
                    BRM["heterozygosity_A"] = BRM["heterozygosity_intra"]
                if sample_sets == "B":
                    del BRM["heterozygosity_A"]
                    BRM["heterozygosity_B"] = BRM["heterozygosity_intra"]
                del BRM["dxy"]
                del BRM["fst"]
            BRMs[sample_sets] = BRM
        return BRMs

    def _get_blocks_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left="[", center="Blocks", right="]", fill="=")
        if self.has_stage("blocks"):
            column_just = 14
            meta_blocks = self._get_meta("blocks")
            BRMs = self._get_blocks_report_metrics()
            reportObj.add_line(prefix="[+]", left="blocks/")
            reportObj.add_line(
                prefix="[+]",
                branch="T",
                fill=".",
                left="'-l %s -m %s -u %s -i %s'"
                % (
                    meta_blocks["length"],
                    meta_blocks["span"],
                    meta_blocks["max_missing"],
                    meta_blocks["max_multiallelic"],
                ),
                right=" %s blocks (%s discarded)"
                % (
                    format_count(meta_blocks["count_total"]),
                    format_percentage(
                        1
                        - (meta_blocks["count_total"] / meta_blocks["count_total_raw"])
                    ),
                ),
            )
            reportObj.add_line(
                prefix="[+]",
                branch="P",
                right="".join([c.rjust(column_just) for c in ["X", "A", "B"]]),
            )
            for label, key, branch in [
                (
                    "BED interval sites in blocks (estimated %)",
                    "interval_coverage",
                    "T",
                ),
                ("Total blocks", "blocks_total", "T"),
                ("Invariant blocks", "blocks_invariant", "T"),
                ("Four-gamete-violation blocks", "blocks_fgv", "T"),
                ("Heterozygosity (population A)", "heterozygosity_A", "T"),
                ("Heterozygosity (population B)", "heterozygosity_B", "T"),
                ("D_xy", "dxy", "T"),
                ("F_st", "fst", "T"),
                ("Pi", "pi", "T"),
                ("Watterson theta", "watterson_theta", "F"),
            ]:
                reportObj.add_line(
                    prefix="[+]",
                    branch=branch,
                    left=label,
                    right="".join(
                        [
                            c.rjust(column_just)
                            for c in [BRMs["X"][key], BRMs["A"][key], BRMs["B"][key]]
                        ]
                    ),
                )
        return reportObj

    def _get_windows_report(self, width):
        meta_windows = self._get_meta("windows")
        reportObj = ReportObj(width=width)
        reportObj.add_line(
            prefix="[+]", left="[", center="Windows", right="]", fill="="
        )
        if self.has_stage("windows"):
            window_size = meta_windows.get(
                "window_size", meta_windows.get("size", "N/A")
            )  # fallback for previous metas
            window_step = meta_windows.get(
                "window_step", meta_windows.get("step", "N/A")
            )  # fallback for previous metas
            window_count = meta_windows.get(
                "window_count", meta_windows.get("count", "N/A")
            )  # fallback for previous metas
            reportObj.add_line(prefix="[+]", left="windows/")
            reportObj.add_line(
                prefix="[+]",
                branch="F",
                fill=".",
                left="'-w %s -s %s'" % (window_size, window_step),
                right="%s windows of inter-population (X) blocks"
                % (format_count(window_count)),
            )
        return reportObj

    def _get_tally_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(
            prefix="[+]", left="[", center="Tallies", right="]", fill="="
        )
        reportObj.add_line(prefix="[+]", left="tally/")
        if "tally/" in self.data:
            shape_by_tally = {tally_label: tally.shape for tally_label, tally in self.data["tally/"].arrays() if not tally_label in set(['blocks_raw', 'windows_raw'])}
            for idx, (tally_label, tally_shape) in enumerate(shape_by_tally.items()):
                branch = "F" if idx == (len(shape_by_tally)-1) else "T"
                reportObj.add_line(
                        prefix="[+]",
                        branch=branch,
                        fill=".",
                        left="'tally/%s'" % (tally_label),
                        right=" shape %s " % (str(tally_shape)),
                    )
        return reportObj

    def info(self, version=None, tree=False):
        width = 100
        if tree:
            return self.data.tree()
        print("[+] \t Getting storage info ...")
        report = self._get_storage_report(width)
        print("[+] \t Getting parse info...")
        report += self._get_parse_report(width)
        print("[+] \t Getting blocks info...")
        report += self._get_blocks_report(width)
        print("[+] \t Getting windows info...")
        report += self._get_windows_report(width)
        print("[+] \t Getting tally info...")
        report += self._get_tally_report(width)
        print("[+] \t Getting simulate info...")
        report += self._get_sims_report(width)
        # report += self._get_optimize_report(width)
        # report += self._get_grids_report(width)
        # report += self._get_lncls_report(width)
        # report += self._get_bsfs_report(width)
        write_info_report(version, report, self.prefix)
        return report

    def _get_grids_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(prefix="[+]", left="[", center="Grids", right="]", fill="=")
        for grid_id in self.data["grids/"]:
            grid_path = "grids/%s" % grid_id
            grid_dict = self.data[grid_path].attrs.asdict()
            reportObj.add_line(
                prefix="[+]",
                branch="W",
                fill=".",
                left=grid_id,
                right=" %s grid points" % format_count(len(grid_dict)),
            )
            # also needs kmax, populations? sync?
            table = []
            for grid_point, grid_params in grid_dict.items():
                table.append([grid_point] + list(grid_params.values()))
            rows = tabulate.tabulate(
                table, numalign="right", headers=["i"] + list(grid_params.keys())
            ).split("\n")
            for i, row in enumerate(rows):
                branch = "F" if (i + 1) == len(rows) else "P"
                reportObj.add_line(
                    prefix="[+]", branch=branch, fill="", left=row, right=""
                )
        return reportObj

    def _get_sims_report(self, width):
        reportObj = ReportObj(width=width)
        reportObj.add_line(
            prefix="[+]", left="[", center="Simulations", right="]", fill="="
        )
        reportObj.add_line(prefix="[+]", left="simulate/")
        if "simulate/" in self.data:
            simulate_labels = list(self.data["simulate"].group_keys())
            for idx, simulate_label in enumerate(simulate_labels):
                branch = "F" if idx == (len(simulate_labels)-1) else "T"
                simulate_key = "simulate/%s" % (simulate_label)
                simulate_meta = self._get_meta(simulate_key)
                reportObj.add_line(
                    prefix="[+]",
                    branch=branch,
                    fill=".",
                    left="%r" % (simulate_key),
                    right="-a %s -b %s -r %s -w %s -b %s -l %s" % (
                        simulate_meta.get('samples_A', "?"),
                        simulate_meta.get('samples_B', "?"),
                        simulate_meta.get('replicates', "?"),
                        simulate_meta.get('windows', "?"),
                        simulate_meta.get('blocks', "?"),
                        simulate_meta.get('block_length', "?"),
                        ),
                )    
        return reportObj
