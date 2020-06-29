# MCMC-based photometry fitting, lens modeling and photo-z code

Based on Matt Auger's pylens and imageSim, and on BPZ. Fits Sersic profiles, lens models and SED templates to multi-band data.

Requirements:
- emcee
- Astropy (only for pyfits)
- PIL (or Pillow)

- need to define an environmental variable "PYPLZDIR", pointing at this directory.

USAGE: 

```
python $PYPLZDIR/pyplz.py OPTION configfile
```

A configuration file consists of a preamble, where keywords indicating the names of the files, filters, length of the MCMC chain etc. are defined, and a models section.

## Options

- -M: runs MCMC using emcee
- -l: reads a list of names of configuration files and runs an MCMC for each one
- -s: saves .fits file with model specified by the configuration file
- -c: reads point from the MCMC chain of the corresponding configuration file, and writes a new configuration file with it as starting point

## Model type: standard or SED

By default, PYPLZ optimizes for the amplitude of each component in each band independently: for a given set of values of the structural parameters of the components (i.e. centroid, half-light radius, axis ratio, position angle, etc.), it finds for each filter the amplitudes of the light components that minimize the chi-squared between the model and the data (in the case of a single light component, this is a linear problem with analytical solution).

Alternatively, one can choose to fit model a spectral energy distribution (SED) to the data. In that case, while the overall amplitude of each component is still optimized, the colors are set by the model SED (they can be varied by means of additional parameters). If one desires to use SED models for the inference, they need to specify it in the preamble of the configuration file, by adding the line 'modeltype: SED'. Since this introduces new parameters on top of the structural ones, convergence will be slower with respect to the default approach.

## Models
The models section can consist of photometric components, mass models (lenses), and lensed sources, in any combination. Photometric components are declared by the keyword 'light_model', followed by the surface brightness profile family (only Sersic is available at the moment), and the keyword specifying the family of SED model, if 'modeltype' is set to 'SED'. The SED model can either be 'colorpars' or 'template'. 'colorpars', which lets one set a uniform prior on the colors of a component with respect to a reference band.

## Parameters
Lines following 'light_model' should list all model parameters describing the component. For Sersic, these are x, y, q (axis ratio), pa, re, n (Sersic index).

Following each parameter name, there should be 5 numbers. The first number indicates the starting value of the parameter in the MCMC (actually, it's the mean of the truncated Gaussian from which the starting values of the parameter are drawn). The second and third number indicate the lower and upper bounds on the parameter. The fourth number indicates the dispersion of the Gaussian from which starting values of the parameter are drawn. The fifth value is an integer: if 0 it means that the parameter is fixed, if 1 the parameter is included in the MCMC.
Optionally, a label indicating a link between parameters can be added after the five numbers.

If 'modeltype' is set to 'SED', each light component will require parameters describing its SED. For 'colorpars', these are the colors in all bands relative to the band specified by the 'reference_band' keyword in the preamble. 
For 'template', the parameters are the redshift zs and the template number tn, which can take any real value between 0 and 7 (non-integer values correspond to linear combinations of two templates). These are the same templates used by BPZ.

## Links
Links to a variable belonging to a light_model are written as 'lightN.par', where N indicates the N-th light_model component in the configuration file, and par is a parameter.
For instance, one may wish to fit two Sersic components to a given galaxy, and impose that the two components have the same centroid.
This can be achieved as follows
light_model Sersic freecolors
x -1 -1 -1 -1 -1 light1.x
y -1 -1 -1 -1 -1 light1.y
In this case the values of the five numbers following the parameter name are ignored, since the values of the corresponding parameters in component light1 are used instead.

## Examples

In the directory examples/ there are three working examples.
The examples are based on a strong lens system, SL2SJ220506+014703 (Sonnenfeld et al. 2013).
HSC g, r, i, z, y-band data of this lens are provided in the examples/data directory.
By running

```
python $PYPLZDIR/pyplz.py -M foreground_standard
```

pyplz will fit a Sersic profile to the lens galaxy, while the background source is masked out. The flux in each filter is optimized independently (no restrictions on colors are imposed).
After the MCMC is completed (it should take around 10 minutes), pyplz outputs a bunch of files.
'foreground_standard_ML' is a new configuration file, where the model is initialized at the maximum likelihood value of the MCMC. 'foreground_standard_ML_rgb.png' and 'foreground_standard_ML.fits' are a PNG and a .fits corresponding to the maximum likelihood model.
The file 'foreground_standard_chain.hdf5' contains the whole chain.
The MCMC is done using emcee. emcee uses a number of walkers to explore the posterior PDF (specified by keyword Nwalkers in configuration file).
In the end the chain will have shape (Nwalkers, Nsteps).

```
python $PYPLZDIR/pyplz.py -M foreground_colorpars
```

same as the previous example, but now the colors of the lens galaxy are explicitly parameterized.

```
python $PYPLZDIR/pyplz.py -M foreground_photoz
```
will fit a Sersic profile to the same galaxy, but this time the model SED is described by a spectral template, and the redshift is a free parameter (zs in the configuration file).
The inferred redshift distribution (stored as 'light1.zs' in the corresponding .hdf5 file), is bimodal. One peak corresponds to the correct solution (the spec-z of this lens is 0.476), while the second peak corresponds to an alternative solution at redshift 3.8.
This second solution should be downweighted by applying a prior on the magnitude-redshift distribution, as done by BPZ, but it's not implemented yet.

