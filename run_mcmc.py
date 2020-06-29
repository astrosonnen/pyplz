import numpy as np
import sys
import h5py
import os
import pyplz_basic
import emcee
from scipy.stats import truncnorm


nargv = len(sys.argv)
configfile = None
old_chainname = None

nkeys = (len(sys.argv) - 1)//2

allowed_keys = ['-M', '-o']

for i in range(nkeys):
    key = sys.argv[2*i+1]
    arg = sys.argv[2*i+2]
    if not key in allowed_keys:
        print('unrecognized option: %s'%key)
        df
    else:
        if key == '-M':
            configfile = arg
        elif key == '-o':
            old_chainname = arg

if configfile is not None:
    print(configfile)

    config = pyplz_basic.read_config(configfile)
    model = pyplz_basic.PyPLZModel(config)

    chainname = config['output_dir']+configfile+'_chain.hdf5'

    npars = len(model.pars)
    nbands = len(model.bands)

    if old_chainname is not None:
        print('starting model from last iteration in %s'%old_chainname)
        old_chain = h5py.File(old_chainname, 'r')
        start = np.zeros((config['Nwalkers'], npars))
        for i in range(npars):
            start[:, i] = old_chain[model.index2par[i]][:, -1]
    else:
        nwalkers = config['Nwalkers']
        start = []
        for j in range(npars):
            a, b = (model.pars[j].lower - model.pars[j].value)/model.pars[j].step, (model.pars[j].upper - model.pars[j].value)/model.pars[j].step
            tmp = truncnorm.rvs(a, b, size=nwalkers)*model.pars[j].step + model.pars[j].value
            start.append(tmp)
        start = np.array(start).T

    def logprior(allpars):
        for i in range(npars):
            if allpars[i] < model.pars[i].lower or allpars[i] > model.pars[i].upper:
                return -np.inf
        return 0.

    fakemags = np.zeros((model.nlight+model.nsource, nbands))

    for i in range(model.nlight):
        #magdic = {}
        #for band in model.bands:
        #    magdic[band] = 99.
        #fakemags.append(magdic)
        for n in range(nbands):
            fakemags[i, n] = 99.

    for i in range(model.nsource):
        #magdic = {}
        #for band in model.bands:
        #    magdic[band] = 99.
        #fakemags.append(magdic)
        for n in range(nbands):
             fakemags[model.nlight+i, n] = 99.

    def logpfunc(allpars):
        lp = logprior(allpars)
        if not np.isfinite(lp):
            return -np.inf, fakemags

        for j in range(npars):
            model.pars[j].value = allpars[j]

        model.update()

        chi2 = model.optimize_amp()
        logp = -0.5*chi2

        if logp != logp:
            return -np.inf, fakemags

        allmags = np.zeros((model.nlight+model.nsource, nbands))

        for i in range(model.nlight):
            #magdic = {}
            #for band in model.light_mags[i]:
            #    magdic[band] = model.light_mags[i][band]
            #allmags.append(magdic)
            for n in range(nbands):
                allmags[i, n] = model.light_mags[i][model.bands[n]]
    
        for i in range(model.nsource):
            #magdic = {}
            #for band in model.source_mags[i]:
            #    magdic[band] = model.source_mags[i][band]
            #allmags.append(magdic)
            for n in range(nbands):
                allmags[model.nlight+i, n] = model.source_mags[i][model.bands[n]]
     
        return logp, allmags

    sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc, threads=config['Nthreads'])

    print("fitting model...")

    sampler.run_mcmc(start, config['Nsteps'])

    chain = sampler.chain
    magschain = sampler.blobs

    outchain = {}
    outchain['logp'] = sampler.lnprobability

    for i in range(npars):
        outchain[model.index2par[i]] = chain[:, :, i]

    for i in range(model.nlight):
        for band in model.bands:
            outchain['light%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, config['Nsteps']))

    for i in range(model.nsource):
        for band in model.bands:
            outchain['source%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, config['Nsteps']))

    for i in range(config['Nsteps']):
        for j in range(nwalkers):
            #for band in model.bands:
            for n in range(nbands):
                for l in range(model.nlight):
                    #outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][l][band]
                    outchain['light%d.mag_%s'%(l+1, model.bands[n])][j, i] = magschain[i][j][l*nbands+n]
                if model.nsource > 0:
                    for s in range(model.nsource):
                        #outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][model.nlight+s][band]
                        outchain['source%d.mag_%s'%(s+1, model.bands[n])][j, i] = magschain[i][j][(model.nlight+s)*nbands+n]

    h5py_file = h5py.File(chainname, 'w')
    for par in outchain:
        h5py_file.create_dataset(par, data=outchain[par])

    ML = outchain['logp'].argmax()

    modelname = config['output_dir']+configfile+'_ML'
    for i in range(len(model.pars)):
        model.pars[i].value = outchain['%s'%model.index2par[i]].flatten()[ML]

    model.update()
    model.optimize_amp()

    model.save(modelname)
    model.write_config_file(config, modelname)

