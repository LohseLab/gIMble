import itertools
import string
import sys
import pathlib
import numpy as np

PARAMETERS_OF_MODEL = {
    "DIV": ['Ne_A_B', 'Ne_A', 'Ne_B', 'T'],
    "MIG_AB": ['Ne_A', 'Ne_B', 'me'],
    "MIG_BA": ['Ne_A', 'Ne_B', 'me'],
    "IM_AB": ['Ne_A_B', 'Ne_A', 'Ne_B', 'me', 'T'],
    "IM_BA": ['Ne_A_B', 'Ne_A', 'Ne_B', 'me', 'T'],
}

MODELS = PARAMETERS_OF_MODEL.keys()

def format_time(seconds):
    hours, remainder = divmod(seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return "{:0>2}h:{:0>2}m:{:06.3f}s".format(int(hours), int(minutes), seconds)
    
class RunArgs(object):
    """Superclass RunArgs"""

    def __init__(self, params):
        self._PATH = params["path"]
        self._VERSION = params["version"]
        self._MODULE = params["module"]
        self._CWD = params["cwd"]

    def __repr__(self):
        return "[+] VER := %s\n[+] CWD := %s\n[+] CMD := %s\n" % (
            self._VERSION,
            self._CWD,
            self._get_cmd(),
        )

    def _get_int(self, string, ret_none=False):
        try:
            return int(string)
        except TypeError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to integer." % string)
            return None
        except ValueError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to integer." % string)
            return None

    def _get_float(self, string, ret_none=False):
        try:
            return float(string)
        except TypeError:
            if not ret_none:
                sys.exit("[X] %r can't be converted to float." % string)
            return None

    def _get_cmd(self):
        return "%s/gimble %s %s" % (
            self._PATH,
            self._MODULE,
            "".join(
                [
                    "--%s " % " ".join((k, str(v)))
                    for k, v in self.__dict__.items()
                    if not k.startswith("_")
                ]
            ),
        )

    def _check_kmax(self, arg):
        l = arg.split(",")
        if len(l) == 4:
            try:
                return np.array([int(_) for _ in l])
            except ValueError:
                pass
        sys.exit("[X] --kmax must be list of four digits, not: %s" % kmax)

    def _check_model(self, model, ret_none=False):
        parameters = {'Ne_A': self.Ne_A, 'Ne_B': self.Ne_B, 'Ne_A_B': self.Ne_A_B, 'me': self.me, 'T': self.T}
        if not model in MODELS:
            if model is None and ret_none:
                return model
            sys.exit("[X] Model %r not supported. Supported models are: %s" % (model, ", ".join(MODELS)))
        PARAMETERS = PARAMETERS_OF_MODEL[model]
        missing_parameters = [parameter for parameter in PARAMETERS if parameters[parameter] is None]
        extra_parameters = [k for k, v in parameters.items() if v is not None and k not in set(PARAMETERS)]
        if missing_parameters:
            print('[X] Model %r requires values for the following parameter(s): %s' % (model, ", ".join(missing_parameters)))
        if extra_parameters:
            print('[X] Model %r does not need the following parameter(s): %s' % (model, ", ".join(extra_parameters)))
        is_me_zero = False
        if model.startswith("MIG"):
            #print('parameters', parameters)
            if isinstance(parameters['me'], list):
                is_me_zero = (parameters['me'][0] == 0)
            else: 
                is_me_zero = (parameters['me'] == 0)
            if is_me_zero:
                print('[X] Model %r only supports non-zero migration rates.' % model)
        if missing_parameters or extra_parameters or is_me_zero:
            sys.exit(1)
        return model

    def _get_path(self, infile, path=False):
        if infile is None:
            return None
        _path = pathlib.Path(infile).resolve()
        if not _path.exists():
            sys.exit("[X] File not found: %r" % str(infile))
        if path:
            return _path
        return str(_path)

    def _get_prefix(self, prefix):
        if prefix is None:
            return None
        path = pathlib.Path(prefix).resolve()
        new_path = pathlib.Path(prefix + ".z").resolve()
        if new_path.exists():
            sys.exit(
                "[X] zstore already exists. Specify using -z or provide new prefix."
            )
        parent_path = pathlib.Path(prefix).resolve().parent
        if not parent_path.exists():
            sys.exit("[X] Path does not exist: %r" % str(parent_path))
        return str(path)