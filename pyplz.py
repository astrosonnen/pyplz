import numpy as np
import sys
import h5py
import os
from pyplz_configtools import read_config
import pyplz_basic, pyplz_sedfit, pyplz_fitters


inputfound = False
nargv = len(sys.argv)
i = 1
key = None
saveind = None
catlist = False
action = None
while i <= nargv and not inputfound:
    argv = sys.argv[i]
    if argv[0] == '-':
        key = argv[1:]
        if key == 'c':
            saveind = int(sys.argv[i+1])
            inputfile = sys.argv[i+2]
            inputfound = True
        elif key == 's':
            inputfile = sys.argv[i+1]
            inputfound = True
        elif key == 'l':
            catlist = True
            inputfile = sys.argv[i+1]
            inputfound = True
        elif key == 'M':
            inputfile = sys.argv[i+1]
            inputfound = True
        else:
            sys.exit('unknown option -%s'%key)

    else:
        sys.exit('must specify mode: -s, -l, -M, -c')

    i += 1

confnames = []
if catlist:
    f = open(inputfile, 'r')
    lines = f.readlines()
    f.close()

    for line in lines:
        confnames.append(line.rstrip())
else:
    confnames = [inputfile]

for configfile in confnames:

    print(configfile)

    config = read_config(configfile)
    if config['modeltype'] == 'standard':
        model = pyplz_basic.PyPLZModel(config)
    elif config['modeltype'] == 'SED':
        model = pyplz_sedfit.PyPLZModel(config)
    else:
        sys.exit("'modeltype' can either be 'standard' or 'SED'")

    if key == 'M' or key == 'l':
   
        chainname = config['output_dir']+configfile+'_chain.hdf5'
        if len(sys.argv) > 3:
            npars = len(model.pars)
            old_chainname = sys.argv[3]
            print('starting model from last iteration in %s'%old_chainname)
            old_chain = h5py.File(old_chainname, 'r')
            walkers = np.zeros((config['Nwalkers'], npars))
            for i in range(npars):
                walkers[:, i] = old_chain[model.index2par[i]][:, -1]
        else:
            walkers = config['Nwalkers']

        pyplz_fitters.run_mcmc(model, chainname, walkers, config['Nsteps'])

        chain = h5py.File(chainname, 'r')

        ML = chain['logp'][()].argmax()
        modelname = config['output_dir']+configfile+'_ML'
        for i in range(len(model.pars)):
            model.pars[i].value = chain['%s'%model.index2par[i]][()].flatten()[ML]

        model.update()
        model.optimize_amp()

        model.save(modelname, config)
        model.write_config_file(config, modelname)

    elif key == 'c':
    
        chain = h5py.File(config['output_dir']+configfile+'_chain.hdf5', 'r')

        modelname = config['output_dir']+configfile+'_%06d'%saveind
        for i in range(len(model.pars)):
            model.pars[i].value = chain['%s'%model.index2par[i]][()].flatten()[saveind]

        model.update()
        model.optimize_amp()

        model.write_config_file(config, modelname)
  
    elif key == 's':

        modelname = config['output_dir']+configfile

        model.update()
        model.optimize_amp()

        model.save(modelname, config)

