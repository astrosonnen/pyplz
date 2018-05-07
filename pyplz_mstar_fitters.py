import emcee
import pymc
import numpy as np
import h5py
from scipy.stats import truncnorm
from scipy.optimize import basinhopping
from scipy.interpolate import splrep, splev, splint
import luminosity_functions
import pyplz_cosmology
from pyplz_cosmology import omegaL, omegaM


def run_mcmc(model, chainname, nwalkers=100, nsteps=1000):

    start = []
    npars = len(model.pars)

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

    fakemags = []

    for i in range(model.nlight):
        magdic = {}
        for band in model.bands:
            magdic[band] = 99.
        fakemags.append(magdic)

    for i in range(model.nsource):
        magdic = {}
        for band in model.bands:
            magdic[band] = 99.
        fakemags.append(magdic)

    def logpfunc(allpars):
        lp = logprior(allpars)
        if not np.isfinite(lp):
            return -np.inf, fakemags

        for j in range(npars):
            model.pars[j].value = allpars[j]

        model.update()

        #logp = -0.5*chi2

        if model.logp != model.logp:
            return -np.inf, fakemags

        allmags = []

        for i in range(model.nlight):
            magdic = {}
            for band in model.bands:
                magdic[band] = model.light_sed_models[i].mags[band]
            allmags.append(magdic)
    
        for i in range(model.nsource):
            magdic = {}
            for band in model.bands:
                magdic[band] = model.source_sed_models[i].mags[band]
            allmags.append(magdic)
     
        return model.logp, allmags

    sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc)

    print "fitting model..."

    sampler.run_mcmc(start, nsteps)

    chain = sampler.chain
    magschain = sampler.blobs

    outchain = {}
    outchain['logp'] = sampler.lnprobability

    for i in range(npars):
        outchain[model.index2par[i]] = chain[:, :, i]

    for i in range(model.nlight):
        for band in model.bands:
            outchain['light%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, nsteps))

    for i in range(model.nsource):
        for band in model.bands:
            outchain['source%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, nsteps))

    for i in range(nsteps):
        for j in range(nwalkers):
            for band in model.bands:
                for l in range(model.nlight):
                    outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][l][band]
                for s in range(model.nsource):
                    outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][model.nlight+s][band]

    h5py_file = h5py.File(chainname, 'w')
    for par in outchain:
        h5py_file.create_dataset(par, data=outchain[par])



