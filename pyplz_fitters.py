import emcee
import pymc
import numpy as np
import h5py
from scipy.stats import truncnorm
from scipy.optimize import basinhopping
from scipy.interpolate import splrep, splev, splint
import luminosity_functions
from luminosity_functions import omegaL, omegaM


nz = 101
z_grid = np.linspace(0., 5., nz)
comovd_grid = 0. * z_grid
for i in range(nz):
    comovd_grid[i] = luminosity_functions.comovd(z_grid[i])

comovd_spline = splrep(z_grid, comovd_grid)
 
def run_mcmc(model, chainname, nwalkers=100, nsteps=1000):

    start = []
    npars = len(model.pars)

    light_photoz_zgrids = {}
    source_photoz_zgrids = {}
    light_photoz_dgrids = {}
    source_photoz_dgrids = {}

    for j in range(npars):

        a, b = (model.pars[j].lower - model.pars[j].value)/model.pars[j].step, (model.pars[j].upper - model.pars[j].value)/model.pars[j].step
        tmp = truncnorm.rvs(a, b, size=nwalkers)*model.pars[j].step + model.pars[j].value
        start.append(tmp)

    start = np.array(start).T

    # checks if there are photoz-s among the free parameters, and prepares grids
    for n in range(model.nlight):
        parname = 'light%d.zs'%(n+1)
        if parname in model.par2index:
            zgrid_here = np.linspace(max(0.001, model.pars[model.par2index[parname]].lower), model.pars[model.par2index[parname]].upper, nz)
            dgrid_here = 0. * zgrid_here
            light_photoz_zgrids[n] = zgrid_here
            for i in range(nz):
                dgrid_here[i] = luminosity_functions.comovd(zgrid_here[i])
            light_photoz_dgrids[n] = dgrid_here

    for n in range(model.nsource):
        parname = 'source%d.zs'%(n+1)
        if parname in model.par2index:
            zgrid_here = np.linspace(max(0.001, model.pars[model.par2index[parname]].lower), model.pars[model.par2index[parname]].upper, nz)
            dgrid_here = 0. * zgrid_here
            source_photoz_zgrids[n] = zgrid_here
            for i in range(nz):
                dgrid_here[i] = luminosity_functions.comovd(zgrid_here[i])
            source_photoz_dgrids[n] = dgrid_here

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
        if i in light_photoz_zgrids:
            magdic['UV'] = 99.
        fakemags.append(magdic)
    for l in light_photoz_zgrids:
        fakemags.append(-np.inf)

    for i in range(model.nsource):
        magdic = {}
        for band in model.bands:
            magdic[band] = 99.
        if i in source_photoz_zgrids:
            magdic['UV'] = 99.
        fakemags.append(magdic)
    for s in source_photoz_zgrids:
        fakemags.append(-np.inf)

    def logpfunc(allpars):
        lp = logprior(allpars)
        if not np.isfinite(lp):
            return -np.inf, fakemags

        for j in range(npars):
            model.pars[j].value = allpars[j]

        model.update()

        chi2 = model.optimize_amp()
        logp = -0.5*chi2

        light_uvmags, source_uvmags = model.get_restUV_mags()

        alluvpriors = []

        for lcomp in light_photoz_zgrids:
            Mgrid_here = light_uvmags[lcomp] - 5.*np.log10(splev(light_photoz_zgrids[lcomp], comovd_spline)*(1.+light_photoz_zgrids[lcomp])*1e5)
            integrand_grid = 1./(omegaL + omegaM*(1.+light_photoz_zgrids[lcomp])**3)**0.5 * light_photoz_dgrids[lcomp]**2 * luminosity_functions.phi(Mgrid_here, light_photoz_zgrids[lcomp])
            integrand_spline = splrep(light_photoz_zgrids[lcomp], integrand_grid)
            norm = splint(light_photoz_zgrids[lcomp][0], light_photoz_zgrids[lcomp][-1], integrand_spline)
            Muv = light_uvmags[lcomp] - 5.*np.log10(splev(model.light_sed_models[lcomp].zs, comovd_spline)*(1.+model.light_sed_models[lcomp].zs)*1e5)
            model.light_mags[lcomp]['UV'] = Muv

            prior = luminosity_functions.phi(Muv, model.light_sed_models[lcomp].zs) * 1./(omegaL + omegaM*(1.+model.light_sed_models[lcomp].zs)**3)**0.5 * splev(model.light_sed_models[lcomp].zs, comovd_spline)**2 / norm
            alluvpriors.append(np.log(prior))
            logp += np.log(prior)

        for scomp in source_photoz_zgrids:
            Mgrid_here = source_uvmags[scomp] - 5.*np.log10(splev(source_photoz_zgrids[scomp], comovd_spline)*(1.+source_photoz_zgrids[scomp])*1e5)
            integrand_grid = 1./(omegaL + omegaM*(1.+source_photoz_zgrids[scomp])**3)**0.5 * source_photoz_dgrids[scomp]**2 * luminosity_functions.phi(Mgrid_here, source_photoz_zgrids[scomp])
            integrand_spline = splrep(source_photoz_zgrids[scomp], integrand_grid)
            norm = splint(source_photoz_zgrids[scomp][0], source_photoz_zgrids[scomp][-1], integrand_spline)
            Muv = source_uvmags[scomp] - 5.*np.log10(splev(model.source_sed_models[scomp].zs, comovd_spline)*(1.+model.source_sed_models[scomp].zs)*1e5)
            model.source_mags[scomp]['UV'] = Muv

            prior = luminosity_functions.phi(Muv, model.source_sed_models[scomp].zs) * 1./(omegaL + omegaM*(1.+model.source_sed_models[scomp].zs)**3)**0.5 * splev(model.source_sed_models[scomp].zs, comovd_spline)**2 / norm
            alluvpriors.append(np.log(prior))
            logp += np.log(prior)

        if logp != logp:
            return -np.inf, fakemags

        allmags = []

        for i in range(model.nlight):
            magdic = {}
            for band in model.light_mags[i]:
                magdic[band] = model.light_mags[i][band]
            allmags.append(magdic)
    
        for i in range(model.nsource):
            magdic = {}
            for band in model.source_mags[i]:
                magdic[band] = model.source_mags[i][band]
            allmags.append(magdic)
     
        return logp, allmags + alluvpriors

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
        if i in light_photoz_zgrids:
            outchain['light%d.M_UV'%(i+1)] = np.zeros((nwalkers, nsteps))
            outchain['light%d.loguvprior'%(i+1)] = np.zeros((nwalkers, nsteps))

    for i in range(model.nsource):
        for band in model.bands:
            outchain['source%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, nsteps))
        if i in source_photoz_zgrids:
            outchain['source%d.M_UV'%(i+1)] = np.zeros((nwalkers, nsteps))
            outchain['source%d.loguvprior'%(i+1)] = np.zeros((nwalkers, nsteps))

    for i in range(nsteps):
        for j in range(nwalkers):
            for band in model.bands:
                for l in range(model.nlight):
                    outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][l][band]
                for s in range(model.nsource):
                    outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][model.nlight+s][band]
            for n in light_photoz_zgrids:
                outchain['light%d.M_UV'%(n+1)][j, i] = magschain[i][j][n]['UV']
                outchain['light%d.loguvprior'%(n+1)][j, i] = magschain[i][j][model.nlight+model.nsource+n]
            for n in source_photoz_zgrids:
                outchain['source%d.M_UV'%(n+1)][j, i] = magschain[i][j][model.nlight+n]['UV']
                outchain['source%d.loguvprior'%(n+1)][j, i] = magschain[i][j][model.nlight+model.nsource+len(light_photoz_zgrids)+n]

    h5py_file = h5py.File(chainname, 'w')
    for par in outchain:
        h5py_file.create_dataset(par, data=outchain[par])


def run_pymc(model, chainname, nsteps=11000, burnin=1000):

    npars = len(model.pars)

    light_photoz_zgrids = {}
    source_photoz_zgrids = {}
    light_photoz_dgrids = {}
    source_photoz_dgrids = {}

    pars = []
    for j in range(npars):
        pars.append(pymc.Uniform(model.index2par[j], lower=model.pars[j].lower, upper=model.pars[j].upper, value=model.pars[j].value))

    # checks if there are photoz-s among the free parameters, and prepares grids
    for n in range(model.nlight):
        parname = 'light%d.zs'%(n+1)
        if parname in model.par2index:
            zgrid_here = np.linspace(max(0.001, model.pars[model.par2index[parname]].lower), model.pars[model.par2index[parname]].upper, nz)
            dgrid_here = 0. * zgrid_here
            light_photoz_zgrids[n] = zgrid_here
            for i in range(nz):
                dgrid_here[i] = luminosity_functions.comovd(zgrid_here[i])
            light_photoz_dgrids[n] = dgrid_here

    for n in range(model.nsource):
        parname = 'source%d.zs'%(n+1)
        if parname in model.par2index:
            zgrid_here = np.linspace(max(0.001, model.pars[model.par2index[parname]].lower), model.pars[model.par2index[parname]].upper, nz)
            dgrid_here = 0. * zgrid_here
            source_photoz_zgrids[n] = zgrid_here
            for i in range(nz):
                dgrid_here[i] = luminosity_functions.comovd(zgrid_here[i])
            source_photoz_dgrids[n] = dgrid_here

    @pymc.deterministic()
    def logpfunc(p=pars):

        for j in range(npars):
            model.pars[j].value = p[j]

        model.update()

        chi2 = model.optimize_amp()
        sumlogp = -0.5*chi2

        light_uvmags, source_uvmags = model.get_restUV_mags()

        for lcomp in light_photoz_zgrids:
            Mgrid_here = light_uvmags[lcomp] - 5.*np.log10(splev(light_photoz_zgrids[lcomp], comovd_spline)*(1.+light_photoz_zgrids[lcomp])*1e5)
            integrand_grid = 1./(omegaL + omegaM*(1.+light_photoz_zgrids[lcomp])**3)**0.5 * light_photoz_dgrids[lcomp]**2 * luminosity_functions.phi(Mgrid_here, light_photoz_zgrids[lcomp])
            integrand_spline = splrep(light_photoz_zgrids[lcomp], integrand_grid)
            norm = splint(light_photoz_zgrids[lcomp][0], light_photoz_zgrids[lcomp][-1], integrand_spline)
            Muv = light_uvmags[lcomp] - 5.*np.log10(splev(model.light_sed_models[lcomp].zs, comovd_spline)*(1.+model.light_sed_models[lcomp].zs)*1e5)
            model.light_mags[lcomp]['UV'] = Muv

            prior = luminosity_functions.phi(Muv, model.light_sed_models[lcomp].zs) * 1./(omegaL + omegaM*(1.+model.light_sed_models[lcomp].zs)**3)**0.5 * splev(model.light_sed_models[lcomp].zs, comovd_spline)**2 / norm
            sumlogp += np.log(prior)

        if sumlogp != sumlogp:
            return -1e200

        return sumlogp

    @pymc.stochastic()
    def logp(p=pars, value=0., observed=True):
        return logpfunc

    M = pymc.MCMC(pars + [logpfunc])
    M.use_step_method(pymc.AdaptiveMetropolis, pars)
    print "fitting model..."

    M.sample(nsteps, burnin)

    outchain = {}
    outchain['logp'] = M.trace('logpfunc')[:]
    for par in pars:
        outchain[str(par)] = M.trace(par)[:]
    
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

    minimizer_kwargs = dict(method="L-BFGS-B", bounds=scale_free_bounds, tol=100.)

    ncomp = model.nlight + model.nsource

    light_photoz_zgrids = {}
    source_photoz_zgrids = {}
    light_photoz_dgrids = {}
    source_photoz_dgrids = {}

    # checks if there are photoz-s among the free parameters, and prepares grids
    for n in range(model.nlight):
        parname = 'light%d.zs'%(n+1)
        if parname in model.par2index:
            zgrid_here = np.linspace(max(0.001, model.pars[model.par2index[parname]].lower), model.pars[model.par2index[parname]].upper, nz)
            dgrid_here = 0. * zgrid_here
            light_photoz_zgrids[n] = zgrid_here
            for i in range(nz):
                dgrid_here[i] = luminosity_functions.comovd(zgrid_here[i])
            light_photoz_dgrids[n] = dgrid_here

    for n in range(model.nsource):
        parname = 'source%d.zs'%(n+1)
        if parname in model.par2index:
            zgrid_here = np.linspace(max(0.001, model.pars[model.par2index[parname]].lower), model.pars[model.par2index[parname]].upper, nz)
            dgrid_here = 0. * zgrid_here
            source_photoz_zgrids[n] = zgrid_here
            for i in range(nz):
                dgrid_here[i] = luminosity_functions.comovd(zgrid_here[i])
            source_photoz_dgrids[n] = dgrid_here

    def nlogpfunc(scaledp):

        p = scaledp * (bounds[:, 1] - bounds[:, 0]) + bounds[:, 0]
        for j in range(npars):
            model.pars[j].value = p[j]

        model.update()

        chi2 = model.optimize_amp()

        logp = -0.5*chi2

        light_uvmags, source_uvmags = model.get_restUV_mags()

        for lcomp in light_photoz_zgrids:
            Mgrid_here = light_uvmags[lcomp] - 5.*np.log10(splev(light_photoz_zgrids[lcomp], comovd_spline)*(1.+light_photoz_zgrids[lcomp])*1e5)
            integrand_grid = 1./(omegaL + omegaM*(1.+light_photoz_zgrids[lcomp])**3)**0.5 * light_photoz_dgrids[lcomp]**2 * luminosity_functions.phi(Mgrid_here, light_photoz_zgrids[lcomp])
            integrand_spline = splrep(light_photoz_zgrids[lcomp], integrand_grid)
            norm = splint(light_photoz_zgrids[lcomp][0], light_photoz_zgrids[lcomp][-1], integrand_spline)
            Muv = light_uvmags[lcomp] - 5.*np.log10(splev(model.light_sed_models[lcomp].zs, comovd_spline)*(1.+model.light_sed_models[lcomp].zs)*1e5)

            prior = luminosity_functions.phi(Muv, model.light_sed_models[lcomp].zs) * 1./(omegaL + omegaM*(1.+model.light_sed_models[lcomp].zs)**3)**0.5 * splev(model.light_sed_models[lcomp].zs, comovd_spline)**2 / norm
            logp += np.log(prior)

        if logp != logp:
            return 1e300
        if not np.isfinite(logp):
            return 1e300
        
        return -logp

    res = basinhopping(nlogpfunc, scale_free_guess, stepsize=0.1, niter=niter, minimizer_kwargs=minimizer_kwargs, interval=50, T=1.)

    ML_pars = res.x*(bounds[:, 1] - bounds[:, 0]) + bounds[:, 0]

    for i in range(npars):
        model.pars[i].value = ML_pars[i]

