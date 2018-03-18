import emcee
import numpy as np
import h5py
from scipy.stats import truncnorm
from scipy.optimize import basinhopping


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

    nlight = len(model.light_sb_models)
    nsource = len(model.source_sb_models)
    ncomp = nlight + nsource

    fakemags = []

    for i in range(ncomp):
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

        chi2 = model.optimize_amp()

        logp = -0.5*chi2

        if logp != logp:
            return -np.inf, fakemags

        return logp, model.light_mags + model.source_mags

    sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc)

    print "fitting model..."

    sampler.run_mcmc(start, nsteps)

    chain = sampler.chain
    magschain = sampler.blobs

    outchain = {}
    outchain['logp'] = sampler.lnprobability

    for i in range(npars):
        outchain[model.index2par[i]] = chain[:, :, i]

    for i in range(nlight):
        for band in model.bands:
            outchain['light%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, nsteps))

    for i in range(nsource):
        for band in model.bands:
            outchain['source%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, nsteps))

    for i in range(nsteps):
        for j in range(nwalkers):
            for band in model.bands:
                for l in range(nlight):
                    outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][l][band]
                for s in range(nsource):
                    outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][nlight+s][band]
    
    h5py_file = h5py.File(chainname, 'w')
    for par in outchain:
        h5py_file.create_dataset(par, data=outchain[par])


def optimize(model, niter=1000):

    start = []
    bounds = []
    npars = len(model.pars)

    for j in range(npars):
        start.append(model.pars[j].value)
        bounds.append((model.pars[j].lower, model.pars[j].upper))

    start = np.array(start)
    bounds = np.array(bounds)

    scale_free_bounds = 0. * bounds
    scale_free_bounds[:, 1] = 1.

    scale_free_guess = (start - bounds[:, 0])/(bounds[:, 1] - bounds[:, 0])

    minimizer_kwargs = dict(method="L-BFGS-B", bounds=scale_free_bounds, tol=1.)

    nlight = len(model.light_sb_models)
    nsource = len(model.source_sb_models)
    ncomp = nlight + nsource

    def nlogpfunc(scaledp):

        p = scaledp * (bounds[:, 1] - bounds[:, 0]) + bounds[:, 0]
        for j in range(npars):
            model.pars[j].value = p[j]

        model.update()

        chi2 = model.optimize_amp()

        logp = -0.5*chi2

        if logp != logp:
            return np.inf

        return -logp

    res = basinhopping(nlogpfunc, scale_free_guess, stepsize=0.1, niter=niter, minimizer_kwargs=minimizer_kwargs, interval=30, T=3.)

    ML_pars = res.x*(bounds[:, 1] - bounds[:, 0]) + bounds[:, 0]

    for i in range(npars):
        model.pars[i].value = ML_pars[i]

