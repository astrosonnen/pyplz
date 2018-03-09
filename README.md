# MCMC-based photometry fitting, lens modeling and photo-z code

Based on Matt Auger's pylens and imageSim, and on BPZ. Fits Sersic profiles, lens models and SED templates to multi-band data.


USAGE: python pyplz.py configfile
- need to define an environmental variable "PYPLZDIR", pointing at this directory.

```
python $PYLENSPATH/pyplz.py configfile
```

## Example

In the directory example/ there is one working example.
It fits a sersic profile to CFHT u, g, r, i, z data of a strong lens.
To run the example, from the example/ directory run

```
python ../pyplz.py sersic_light
```

It will produce a bunch of output files.

sersic_light.out is a new configuration file with the same setting as the original, except for the starting point of the model, which is set to the maximum likelihood of the MCMC chain.

The MCMC is done using emcee. emcee uses a number of walkers to explore the posterior PDF (specified by keyword Nwalkers in configuration file).
In the end the chain will have shape (Nwalkers, Nsteps).


