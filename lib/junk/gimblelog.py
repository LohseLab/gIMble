import logging
import sys

# CONSTANTS
SIMPLE_LOG_FMT = '%(message)s'
DEBUG_LOG_DATEFMT = '%m/%d/%Y %H:%M:%S'
DEBUG_LOG_FMT = '[%(asctime)s] [%(levelname)-5s] [%(name)s] %(message)s'

def get_logger(run_params, debug):
    if debug:
        logging.basicConfig(
            level=logging.DEBUG, 
            format=DEBUG_LOG_FMT, 
            datefmt=DEBUG_LOG_DATEFMT, 
            filename="gIMble_v%s.%s.debug.log" % (run_params['version'], run_params['module']), 
            filemode='w')
    else:
        logging.basicConfig(
            level=logging.INFO, 
            format=SIMPLE_LOG_FMT, 
            filename="gIMble_v%s.%s.log" % (run_params['version'], run_params['module']),
            filemode='w')
    log = logging.getLogger(run_params['module'])
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter(SIMPLE_LOG_FMT))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    log.info("[V] gIMble v%s" % (run_params['version']))
    return log