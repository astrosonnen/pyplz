import numpy as np
import sys
import h5py
import os
import pyplz_objects, pyplz_fitters


inputfound = False
nargv = len(sys.argv)
i = 1
key = None
saveind = None
catlist = False
dofit = True
while i <= nargv and not inputfound:
    argv = sys.argv[i]
    if argv[0] == '-':
        key = argv[1:]
        if key == 's':
            saveind = int(sys.argv[i+1])
            inputfile = sys.argv[i+2]
            inputfound = True
            dofit = False
        elif key == 'l':
            catlist = True
            inputfile = sys.argv[i+1]
            inputfound = True
    else:
        inputfile = argv
        inputfound = True
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

    config = pyplz_objects.read_config(configfile)
    model = pyplz_objects.PyPLZModel(config)

    if dofit:
   
        chainname = config['output_dir']+configfile+'_chain.hdf5'
        pyplz_fitters.run_mcmc(model, chainname, config['Nwalkers'], config['Nsteps'])

        chain = h5py.File(chainname, 'r')

        ind = chain['logp'].value.argmax()
        modelname = config['output_dir']+configfile+'_ML'

    elif key == 's':
    
        chain = h5py.File(config['output_dir']+configfile+'_chain.hdf5', 'r')

        modelname = config['output_dir']+configfile+'_%06d'%ind

    for i in range(len(model.pars)):
        model.pars[i].value = chain['%s'%model.index2par[i]].value.flatten()[ind]

    model.update()
    model.optimize_amp()
    
    model.save(modelname)
    model.write_config_file(config, modelname)

