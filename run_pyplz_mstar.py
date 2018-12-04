import numpy as np
import sys
import h5py
import os
import pyplz_mstar_objects, pyplz_mstar_fitters


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
            print 'unknown option -%s'%key
            df
    else:
        print 'must specify mode: -s, -l, -M, -c'
        df
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

    print configfile

    config = pyplz_mstar_objects.read_config(configfile)
    model = pyplz_mstar_objects.PyPLZModel(config)

    if key == 'M' or key == 'l':
   
        chainname = config['output_dir']+configfile+'_chain.hdf5'
        if len(sys.argv) > 3:
            npars = len(model.pars)
            old_chainname = sys.argv[3]
            print 'starting model from last iteration in %s'%old_chainname
            old_chain = h5py.File(old_chainname, 'r')
            walkers = np.zeros((config['Nwalkers'], npars))
            for i in range(npars):
                walkers[:, i] = old_chain[model.index2par[i]][:, -1]
        else:
            walkers = config['Nwalkers']

        pyplz_mstar_fitters.run_mcmc(model, chainname, walkers, config['Nsteps'])

        chain = h5py.File(chainname, 'r')

        ML = chain['logp'].value.argmax()
        modelname = config['output_dir']+configfile+'_ML'
        for i in range(len(model.pars)):
            model.pars[i].value = chain['%s'%model.index2par[i]].value.flatten()[ML]

        model.update()

        model.save(modelname)
        model.write_config_file(config, modelname)

    elif key == 'c':
    
        chain = h5py.File(config['output_dir']+configfile+'_chain.hdf5', 'r')

        modelname = config['output_dir']+configfile+'_%06d'%saveind
        for i in range(len(model.pars)):
            model.pars[i].value = chain['%s'%model.index2par[i]].value.flatten()[saveind]

        model.update()

        model.write_config_file(config, modelname)
  
    elif key == 's':

        modelname = config['output_dir']+configfile

        model.update()

        model.save(modelname)


