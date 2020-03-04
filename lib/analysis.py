class ParameterObj(object):
    def __init__(self, args):
        print(args)
        self.module = 'analysis'
        self.zstore = args['--zarr']