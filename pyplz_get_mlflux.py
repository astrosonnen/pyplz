import numpy as np
import sys
import h5py
import os
import debug_mlflux_objects as pyplz_mlflux_objects
import emcee
from scipy.stats import truncnorm


nargv = len(sys.argv)
configfile = None
old_chainname = None

nkeys = (len(sys.argv) - 1)/2

allowed_keys = ['-c']

for i in range(nkeys):
    key = sys.argv[2*i+1]
    arg = sys.argv[2*i+2]
    if not key in allowed_keys:
        print 'unrecognized option: %s'%key
        df
    else:
        if key == '-c':
            configfile = arg

if configfile is not None:
    print configfile

    config = pyplz_mlflux_objects.read_config(configfile)
    model = pyplz_mlflux_objects.PyPLZModel(config)

    model.update()
    model.optimize_amp()

    modelname = configfile+'_ML'
    model.save(modelname)
    model.write_config_file(config, modelname)

