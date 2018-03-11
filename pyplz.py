import numpy as np
import sys
from pylens import pylens, MassModels
from imageSim import SBModels, SEDModels, convolve, plotting_tools
import h5py
from scipy.stats import truncnorm
from scipy.interpolate import splrep, splev
import emcee
from scipy.optimize import nnls
import os
from astropy.io import fits as pyfits


inputfound = False
nargv = len(sys.argv)
i = 1
key = None
saveind = None
catlist = False
while i <= nargv and not inputfound:
    argv = sys.argv[i]
    if argv[0] == '-':
        key = argv[1:]
        if key == 's':
            saveind = int(sys.argv[i+1])
            inputfile = sys.argv[i+2]
            inputfound = True
        elif key == 'l':
            catlist = True
            inputfile = sys.argv[i+1]
            inputfound = True
    else:
        inputfile = argv
        inputfound = True
    i += 1

rootdir = os.environ.get('PYLENSDIR')

class Par:

    def __init__(self, name, lower=0., upper=1., value=0.):

        self.name = name
        self.lower = lower
        self.upper = upper
        self.value = value

def read_config(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    config = {'data_dir':'./', 'mask_dir': None, 'output_dir':'./', 'filters': None, 'main_band': None, 'fitbands': None, \
              'zeropoints': None, \
              'filename': None, 'filter_prefix': '', 'filter_suffix': '', 'science_tag':'_sci.fits', 'err_tag':'_var.fits', 'err_type': 'VAR', 'psf_tag':'_psf.fits', \
              'rmax': None, 'Nsteps':300, 'Nwalkers':30, 'burnin':None, 'maskname':None, \
              'rgbcuts': None, 'outname': None}

    preamble = True

    i = 0
    while preamble and i < len(lines):
        if '#' in lines[i] and 'MODELS' in lines[i]:
            preamble = False
        else:
            line = lines[i].split('#')[0].split()
            if len(line) > 0:
                parname = line[0].split(':')[0]
                if parname in config:
                    config[parname] = lines[i].split('#')[0].split(':')[1].split('\n')[0].lstrip().rstrip()
        i += 1

    filtlist = []
    filternames = config['filters'].split(',')
    for name in filternames:
        filtlist.append(name.lstrip())
    config['filters'] = filtlist
    if config['fitbands'] is not None:
        filtlist = []
        filternames = config['fitbands'].split(',')
        for name in filternames:
            filtlist.append(name.lstrip())
        config['fitbands'] = filtlist
    else:
        config['fitbands'] = config['filters']

    config['colors'] = []
    for band in config['fitbands']:
        if band != config['main_band']:
            config['colors'].append('%s-%s'%(band, config['main_band']))

    if config['rgbcuts'] is not None:
        cutlist = []
        cuts = config['rgbcuts'].split(',')
        for cut in cuts:
            cutlist.append(float(cut))

        config['rgbcuts'] = cutlist
    else:
        config['rgbcuts'] = (99., 99., 99.)

    if config['outname'] is None:
        config['outname'] = configfile+'.output'

    config['zeropoints'] = np.array(config['zeropoints'].split(','), dtype='float')
    config['Nsteps'] = int(config['Nsteps'])
    config['Nwalkers'] = int(config['Nwalkers'])
    config['rmax'] = float(config['rmax'])

    light_components = []
    lens_components = []
    source_components = []

    while i < len(lines):

        line = lines[i].split()
        if len(line) > 0:
            if line[0] == 'light_model':
                model_class = line[1]
                sed_class = line[2]

                if model_class == 'Sersic':
                    parnames = ['x', 'y', 'q', 'pa', 're', 'n']
                else:
                    df

                if sed_class == 'freecolors':
                    parnames += config['colors']
                elif sed_class == 'template':
                    parnames += ['zs', 'tn']
                else:
                    df

                npars = len(parnames)

                comp = {'class':model_class, 'pars':{}, 'sed': sed_class}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print 'not all parameters found!'
                else:
                    light_components.append(comp)

            elif line[0] == 'source_model':
                model_class = line[1].lstrip()
                sed_class = line[2]

                if model_class == 'Sersic':
                    parnames = ['x', 'y', 'q', 'pa', 're', 'n']
                else:
                    df

                if sed_class == 'freecolors':
                    parnames += config['colors']
                else:
                    parnames += ['zs', 'tn']

                npars = len(parnames)

                comp = {'class':model_class, 'sed': sed_class, 'pars':{}}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print 'not all parameters found!'
                else:
                    source_components.append(comp)

            elif 'lens_model' in line[0]:
                model_class = line[1].lstrip()

                if model_class == 'Powerlaw':
                    npars = 6
                    parnames = ['x', 'y', 'q', 'pa', 'b', 'eta']
                else:
                    df

                comp = {'class': model_class, 'sed': sed_class, 'pars':{}}

                foundpars = 0
                j = 1

                while foundpars < npars and j+i < len(lines):
                    line = lines[j+i].split()
                    if lines[j+i][0] != '#' and len(line) > 0:
                        if line[0] in parnames:
                            foundpars += 1
                            par = line[0]
                            link = None
                            if len(line) > 6:
                                link = line[6]
                            tmp_par = {'value': float(line[1]), 'low': float(line[2]), 'up': float(line[3]), \
                           'step': float(line[4]), 'var': int(line[5]), 'link':link}
                            comp['pars'][par] = tmp_par
                    j += 1

                i += j

                if foundpars < npars:
                    print 'not all parameters found!'
                else:
                    lens_components.append(comp)

            else:
                i += 1
        else:
            i += 1

    config['light_components'] = light_components
    config['source_components'] = source_components
    config['lens_components'] = lens_components

    return config

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

    config = read_config(configfile)
    
    images = {}
    sigmas = {}
    convol_matrix = {}
    light_models = []
    source_models = []
    lens_models = []
    light_seds = []
    source_seds = []
    zp = {}
    pars = []
    bounds = []
    steps = []
    filters = [filt for filt in config['filters']]
    fitbands = [band for band in config['fitbands']]
    
    #defines model parameters
    par2index = {}
    index2par = {}
    
    nlight = len(config['light_components'])
    nsource = len(config['source_components'])
    nlens = len(config['lens_components'])
    
    ncomp = 0
    npar = 0
    for comp in config['light_components']:
        ncomp += 1
        for par in comp['pars']:
            parpar = comp['pars'][par]
            if parpar['link'] is None and parpar['var'] == 1:
                pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
                bounds.append((parpar['low'], parpar['up']))
                steps.append(parpar['step'])
                par2index['light'+str(ncomp)+'.'+par] = npar
                index2par[npar] = 'light'+str(ncomp)+'.'+par
                npar += 1
    
    ncomp = 0
    for comp in config['source_components']:
        ncomp += 1
        for par in comp['pars']:
            parpar = comp['pars'][par]
            if parpar['link'] is None and parpar['var'] == 1:
                pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
                bounds.append((parpar['low'], parpar['up']))
                steps.append(parpar['step'])
                par2index['source'+str(ncomp)+'.'+par] = npar
                index2par[npar] = 'source'+str(ncomp)+'.'+par
                npar += 1
    
    ncomp = 0
    for comp in config['lens_components']:
        ncomp += 1
        for par in comp['pars']:
            parpar = comp['pars'][par]
            if parpar['link'] is None and parpar['var'] == 1:
                pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value']))
                bounds.append((parpar['low'], parpar['up']))
                steps.append(parpar['step'])
                par2index['lens'+str(ncomp)+'.'+par] = npar
                index2par[npar] = 'lens'+str(ncomp)+'.'+par
                npar += 1
    
    npars = len(pars)
    
    i = 0
    
    filtdic = {}
    for band in config['filters']:
    
        zp[band] = config['zeropoints'][i]
    
        filtname = rootdir+'pylens/filters/%s%s%s'%(config['filter_prefix'], band, config['filter_suffix'])
    
        f = open(filtname, 'r')
        ftable = np.loadtxt(f)
        f.close()
    
        filtdic[band] = (ftable[:, 0], ftable[:, 1])
    
        hdu = pyfits.open(config['data_dir']+'/'+config['filename']+'_%s'%band+config['science_tag'])[0]
    
        img = hdu.data.copy()
        subimg = img.copy()
        images[band] = subimg
        suberr = pyfits.open(config['data_dir']+'/'+config['filename']+'_%s'%band+config['err_tag'])[0].data.copy()
    
        if config['err_type'] == 'VAR':
            sigmas[band] = suberr**0.5
        elif config['err_type'] == 'SIGMA':
            sigmas[band] = suberr.copy()
        else:
            df
    
        psf_file = pyfits.open(config['data_dir']+'/'+config['filename']+'_%s'%band+config['psf_tag'])
        nhdu = len(psf_file)
        found = False
        n = 0
        while not found and n < nhdu:
            if psf_file[n].data is not None:
                psf = psf_file[n].data.copy()
                found = True
            else:
                n += 1
        if not found:
            df
    
        m = (psf[:2].mean()+psf[-2:].mean()+psf[:,:2].mean()+psf[:,-2:].mean())/4.
        psf -= m
        psf /= psf.sum()
    
        convol_matrix[band] = convolve.convolve(images[band], psf)[1]
    
    ncomp = 1
    for comp in config['lens_components']:
    
        name = 'lens%d'%ncomp
    
        pars_here = {}
        for par in comp['pars']:
            if comp['pars'][par]['link'] is None:
                if comp['pars'][par]['var'] == 1:
                    pars_here[par] = pars[par2index[name+'.'+par]]
                elif comp['pars'][par]['var'] == 0:
                    pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                else:
                    df
            else:
                pars_here[par] = pars[par2index[comp['pars'][par]['link']]]
    
        ncomp += 1
    
        lens = MassModels.PowerLaw(name, pars_here)
        lens_models.append(lens)
    
    ncomp = 1
    light_pardicts = []
    light_seddicts = []
    sersicpars = ['x', 'y', 'q', 're', 'pa', 'n']
    for comp in config['light_components']:
    
        name = 'light%d'%ncomp
    
        pars_here = {}
        sed_here = {}
        for par in comp['pars']:
            if par in sersicpars:
                if comp['pars'][par]['link'] is None:
                    if comp['pars'][par]['var'] == 1:
                        pars_here[par] = pars[par2index[name+'.'+par]]
                    else:
                        pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                else:
                    pars_here[par] = pars[par2index[comp['pars'][par]['link']]]
            else:
                if comp['pars'][par]['link'] is None:
                    if comp['pars'][par]['var'] == 1:
                        sed_here[par] = pars[par2index[name+'.'+par]]
                    else:
                        sed_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                else:
                    sed_here[par] = pars[par2index[comp['pars'][par]['link']]]
    
        ncomp += 1
    
        light_pardicts.append(pars_here)
        light_seddicts.append(sed_here)
    
    ncomp = 1
    source_pardicts = []
    source_seddicts = []
    for comp in config['source_components']:
    
        name = 'source%d'%ncomp
    
        pars_here = {}
        for par in comp['pars']:
            if par in sersicpars:
                if comp['pars'][par]['link'] is None:
                    if comp['pars'][par]['var'] == 1:
                        pars_here[par] = pars[par2index[name+'.'+par]]
                    else:
                        pars_here[par] = comp['pars'][par]['value']
                else:
                    pars_here[par] = pars[par2index[comp['pars'][par]['link']]]
            else:
                if comp['pars'][par]['link'] is None:
                    if comp['pars'][par]['var'] == 1:
                        sed_here[par] = pars[par2index[name+'.'+par]]
                    else:
                        sed_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                else:
                    sed_here[par] = pars[par2index[comp['pars'][par]['link']]]
    
        ncomp += 1
    
        source_pardicts.append(pars_here)
        source_seddicts.append(sed_here)
    
    ncomp = 1
    for pardict, seddict, comp in zip(light_pardicts, light_seddicts, config['light_components']):
    
        light = SBModels.Sersic('light%d'%ncomp, pardict)
        light_models.append(light)
        if comp['sed'] == 'freecolors':
            sed = SEDModels.Colors('light_sed%d'%ncomp, seddict, zp)
        elif comp['sed'] == 'template':
            sed = SEDModels.Template('light_sed%d'%ncomp, seddict, filtdic, zp)
        light_seds.append(sed)
        ncomp += 1
    
    ncomp = 1
    for pardict, seddict, comp in zip(source_pardicts, source_seddicts, config['source_components']):
    
        source = SBModels.Sersic('source%d'%ncomp, pardict)
        source_models.append(source)
        if comp['sed'] == 'freecolors':
            sed = SEDModels.Colors('light_sed%d'%ncomp, seddict, filtdic, zp=zp)
        elif comp['sed'] == 'template':
            sed = SEDModels.Template('light_sed%d'%ncomp, seddict, filtdic, zp=zp)
        source_seds.append(sed)
        ncomp += 1
    
    ny, nx = images[filters[0]].shape
    X, Y = np.meshgrid(np.arange(1.*nx), np.arange(1.*ny))
    R = ((X - nx/2)**2 + (Y - ny/2)**2)**0.5
    
    if config['mask_dir'] is None:
        mask_dir = config['data_dir']
    else:
        mask_dir = './'
    
    if config['maskname'] is not None:
        MASK = pyfits.open(mask_dir+config['maskname'])[0].data.copy()
    else:
        MASK = np.ones(X.shape, dtype=int)
    
    if config['rmax'] is not None:
        MASK[R > config['rmax']] = 0
    
    mask = MASK > 0
    mask_r = mask.ravel()
    
    output = {}
    
    nfitbands = len(fitbands)
    
    modelstack = np.zeros((nfitbands * ny, nx))
    datastack = 0. * modelstack
    sigmastack = 0. * modelstack
    maskstack = np.zeros((nfitbands * ny, nx), dtype=bool)
    for i in range(nfitbands):
        datastack[i*ny: (i+1)*ny, :] = images[fitbands[i]]
        sigmastack[i*ny: (i+1)*ny, :] = sigmas[fitbands[i]]
        maskstack[i*ny: (i+1)*ny, :] = mask
    maskstack_r = maskstack.ravel()
    
    def save_model(config, chain, ind=None):
    
        if ind is None:
            ind = chain['logp'].value.argmax()
            fitsname = config['output_dir']+configfile+'_ML.fits'
            rgbname = config['output_dir']+configfile+'_ML_rgb.png'
        else:
            fitsname = config['output_dir']+configfile+'_%06d.fits'%ind
            rgbname = config['output_dir']+configfile+'_%06d_rgb.png'%ind
    
        for par in par2index:
            pars[par2index[par]].value = chain[par].value.flatten()[ind]
    
        hdr = pyfits.Header()
    
        light_ind_model = []
        source_ind_model = []
        
        # saves best fit model images
        for lens in lens_models:
            lens.setPars()
        
        xl, yl = pylens.getDeflections(lens_models, (X, Y))
        
        modlist = []
        
        ntotbands = len(filters)
        
        modelstack = np.zeros((ntotbands * ny, nx))
        datastack = 0. * modelstack
        sigmastack = 0. * modelstack
        maskstack = np.zeros((ntotbands * ny, nx), dtype=bool)
        
        for i in range(ntotbands):
            datastack[i*ny: (i+1)*ny, :] = images[filters[i]]
            sigmastack[i*ny: (i+1)*ny, :] = sigmas[filters[i]]
            maskstack[i*ny: (i+1)*ny, :] = mask
        
        maskstack_r = maskstack.ravel()
        
        light_maglist = []
    
        for light, sed in zip(light_models, light_seds):
            light.setPars()
            sed.setPars()
            mags = {}
            light.amp = 1.
            lpix = light.pixeval(X, Y)
            lmodel = 0. * datastack
            for i in range(ntotbands):
                scale = sed.scale(filters[i], config['main_band'])
                lmodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(lpix, convol_matrix[filters[i]], False)[0]
                mags[filters[i]] = -2.5*np.log10(scale) + zp[filters[i]] - zp[config['main_band']]
    
            modlist.append((lmodel/sigmastack).ravel()[maskstack_r])
            light_maglist.append(mags)
    
        source_maglist = []
        for source, sed in zip(source_models, source_seds):
            source.setPars()
            sed.setPars()
            mags = {}
            source.amp = 1.
            spix = source.pixeval(xl, yl)
            smodel = 0. * datastack
            for i in range(ntotbands):
                scale = sed.scale(filters[i], config['main_band'])
                smodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(spix, convol_matrix[filters[i]], False)[0]
                mags[filters[i]] = -2.5*np.log10(scale) + zp[filters[i]] - zp[config['main_band']]
            modlist.append((smodel/sigmastack).ravel()[maskstack_r])
            source_maglist.append(mags)
        
        modarr = np.array(modlist).T
        
        amps, chi = nnls(modarr, (datastack/sigmastack).ravel()[maskstack_r])
        
        n = 0
        for light, sed, mags in zip(light_models, light_seds, light_maglist):
        
            light_ind_dic = {}
    
            light.amp = amps[n]
            lpix = light.pixeval(X, Y)
            mainmag = light.Mag(zp[config['main_band']])
            for band in filters:
                scale = sed.scale(band, config['main_band'])
                light_ind_here = scale * convolve.convolve(lpix, convol_matrix[band], False)[0]
                mags[band] += mainmag
                hdr['%s.mag_%s'%(light.name, band)] = mags[band]
                
                light_ind_dic[band] = light_ind_here
    
            light_ind_model.append(light_ind_dic)
    
            n += 1
    
        for source, sed, mags in zip(source_models, source_seds, source_maglist):
        
            source_ind_dic = {}
    
            source.amp = amps[n]
            spix = source.pixeval(xl, yl)
            mainmag = source.Mag(zp[config['main_band']])
        
            for band in filters:
                scale = sed.scale(band, config['main_band'])
                mags[band] += mainmag
                source_ind_here = scale * convolve.convolve(spix, convol_matrix[band], False)[0]
                hdr['%s.mag_%s'%(source.name, band)] = mags[band]
    
                source_ind_dic[band] = source_ind_here
    
            source_ind_model.append(source_ind_dic)
        
            n += 1
    
        hdr['CHI2'] = chi**2
    
        phdu = pyfits.PrimaryHDU(header=hdr)
    
        hdulist = pyfits.HDUList([phdu])
    
        for light, light_ind_dic in zip(light_models, light_ind_model):
            
            for band in filters:
                hdu_here = pyfits.ImageHDU(data=light_ind_dic[band])
                hdu_here.header['EXTNAME'] = '%s_%s'%(light.name, band)
    
                hdulist.append(hdu_here)
    
        for source, source_ind_dic in zip(source_models, source_ind_model):
            for band in filters:
                hdu_here = pyfits.ImageHDU(data=source_ind_dic[band])
                hdu_here.header['EXTNAME'] = '%s_%s'%(source.name, band)
    
                hdulist.append(hdu_here)
    
        hdulist.writeto(fitsname, overwrite=True)
    
        modlist = []
        for light_ind_dic in light_ind_model:
            lmodel = 0. * datastack
            for i in range(ntotbands):
                lmodel[i*ny: (i+1)*ny, :] = light_ind_dic[filters[i]]
            modlist.append((lmodel/sigmastack).ravel()[maskstack_r])
    
        for source_ind_dic in source_ind_model:
            smodel = 0. * datastack
            for i in range(ntotbands):
                smodel[i*ny: (i+1)*ny, :] = source_ind_dic[filters[i]]
            modlist.append((smodel/sigmastack).ravel()[maskstack_r])
    
        modarr = np.array(modlist).T
    
        amps, chi = nnls(modarr, (datastack/sigmastack).ravel()[maskstack_r])
    
        # makes model rgb
        sci_list = []
        light_list = []
        source_list = []
        for band in filters:
            sci_list.append(images[band])
            lmodel = 0.*images[band]
            smodel = 0.*images[band]
            for light in light_ind_model:
                lmodel += light[band]
            light_list.append(lmodel)
            for source in source_ind_model:
                smodel += source[band]
            source_list.append(smodel)
        
        plotting_tools.make_full_rgb(sci_list, light_list, source_list, outname=rgbname)
    
    # writes a new configuration file
    def write_config_file(config, chain, ind=None):
    
        if ind is None:
            ind = chain['logp'].value.argmax()
            outconfname = config['output_dir']+configfile+'_ML'
        else:
            outconfname = config['output_dir']+configfile+'_%06d'%ind
    
        for par in par2index:
            pars[par2index[par]].value = chain[par].value.flatten()[ind]
    
        conflines = []
        confpars = ['data_dir', 'mask_dir', 'maskname', 'output_dir', 'filename', 'science_tag', 'err_tag', 'err_type', 'psf_tag', 'rmax', 'Nwalkers', 'Nsteps', 'main_band', 'filter_prefix', 'filter_suffix']
        for parname in confpars:
            if config[parname] is not None:
                conflines.append('%s: %s\n'%(parname, config[parname]))
        filtline = 'filters: '
        zpline = 'zeropoints: '
        nfilt = len(filters)
        for i in range(nfilt-1):
            filtline += '%s, '%filters[i]
            zpline += '%f, '%zp[filters[i]]
        filtline += filters[-1]+'\n'
        zpline += '%f\n'%zp[filters[-1]]
        conflines.append(filtline)
        conflines.append(zpline)
        
        filtline = 'fitbands: '
        nfilt = len(fitbands)
        for i in range(nfilt-1):
            filtline += '%s, '%fitbands[i]
        filtline += fitbands[-1]+'\n'
        conflines.append(filtline)
        
        conflines.append('\n')
        conflines.append('# MODELS\n')
        
        ncomp = 0
        for comp in config['light_components']:
    
            if comp['sed'] == 'freecolors':
                lightpars = sersicpars + config['colors']
            elif comp['sed'] == 'template':
                lightpars = sersicpars + ['zs', 'tn']
    
            conflines.append('\n')
            conflines.append('light_model Sersic %s\n'%comp['sed'])
            for par in lightpars:
                parname = 'light%d.%s'%(ncomp+1, par)
                if parname in par2index:
                    npar = par2index[parname]
                    conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
                else:
                    if config['light_components'][ncomp]['pars'][par]['link'] is None:
                        conflines.append('%s %f -1 -1 -1 0\n'%(par, config['light_components'][ncomp]['pars'][par]['value']))
                    else:
                        lname = config['light_components'][ncomp]['pars'][par]['link']
                        npar = par2index[lname]
                        conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
            for band in filters:
                conflines.append('mag_%s %3.2f\n'%(band, chain['light%d.mag_%s'%(ncomp+1, band)].value.flatten()[ind]))
            ncomp += 1
        
        ncomp = 0
    
        for comp in config['source_components']:
    
            if comp['sed'] == 'freecolors':
                sourcepars = sersicpars + config['colors']
            elif comp['sed'] == 'template':
                sourcepars = sersicpars + ['zs', 'tn']
    
            conflines.append('\n')
            conflines.append('source_model Sersic %s\n'%comp['sed'])
            for par in sourcepars:
                parname = 'source%d.%s'%(ncomp+1, par)
                if parname in par2index:
                    npar = par2index[parname]
                    conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
                else:
                    if config['source_components'][ncomp]['pars'][par]['link'] is None:
                        conflines.append('%s %f -1 -1 -1 0\n'%(par, config['source_components'][ncomp]['pars'][par]['value']))
                    else:
                        lname = config['source_components'][ncomp]['pars'][par]['link']
                        npar = par2index[lname]
                        conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
            for band in filters:
                conflines.append('mag_%s %3.2f\n'%(band, chain['source%d.mag_%s'%(ncomp+1, band)].value.flatten()[ind]))
    
            ncomp += 1
        
        ncomp = 0
        powpars = ['x', 'y', 'pa', 'q', 'b', 'eta']
        
        for lens in lens_models:
            conflines.append('\n')
            conflines.append('lens_model Powerlaw\n')
            for par in powpars:
                parname = 'lens%d.%s'%(ncomp+1, par)
                if parname in par2index:
                    npar = par2index[parname]
                    conflines.append('%s %f %f %f %f 1\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar]))
                else:
                    if config['lens_components'][ncomp]['pars'][par]['link'] is None:
                        conflines.append('%s %f -1 -1 -1 0\n'%(par, config['lens_components'][ncomp]['pars'][par]['value']))
                    else:
                        lname = config['lens_components'][ncomp]['pars'][par]['link']
                        npar = par2index[lname]
                        conflines.append('%s %f %f %f %f 1 %s\n'%(par, pars[npar].value, bounds[npar][0], bounds[npar][1], steps[npar], lname))
            ncomp += 1
        
        conflines.append('\n')
        conflines.append('logp %f\n'%(chain['logp'].value.flatten()[ind]))
        
        f = open(outconfname, 'w')
        f.writelines(conflines)
        f.close()
    
    if key != 's':
    
        start = []
        for j in range(npars):
            a, b = (bounds[j][0] - pars[j].value)/steps[j], (bounds[j][1] - pars[j].value)/steps[j]
            tmp = truncnorm.rvs(a, b, size=config['Nwalkers'])*steps[j] + pars[j].value
    
            start.append(tmp)
    
        start = np.array(start).T
    
        npars = len(pars)
       
        def logprior(allpars):
            for i in range(npars):
                if allpars[i] < bounds[i][0] or allpars[i] > bounds[i][1]:
                    return -np.inf
            return 0.
    
        nwalkers = len(start)
    
        fakemags = []
        for comp in light_models + source_models:
            magdic = {}
            for band in fitbands:
                magdic[band] = 99.
            fakemags.append(magdic)
    
        ncomp = len(light_models) + len(source_models)
    
        def logpfunc(allpars):
            lp = logprior(allpars)
            if not np.isfinite(lp):
                return -np.inf, fakemags
    
            for j in range(0, npars):
                pars[j].value = allpars[j]
    
            for lens in lens_models:
                lens.setPars()
    
            xl, yl = pylens.getDeflections(lens_models, (X, Y))
    
            modlist = []
            light_maglist = []
    
            for light, sed in zip(light_models, light_seds):
                lmodel = 0. * datastack
    
                light.setPars()
                sed.setPars()
                light.amp = 1.
                lpix = light.pixeval(X, Y)
                mags = {}
                for i in range(nfitbands):
                    scale = sed.scale(fitbands[i], config['main_band'])
                    lmodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(lpix, convol_matrix[fitbands[i]], False)[0]
                    mags[fitbands[i]] = -2.5*np.log10(scale) + zp[fitbands[i]] - zp[config['main_band']]
    
                modlist.append((lmodel/sigmastack).ravel()[maskstack_r])
                light_maglist.append(mags)
    
            source_maglist = []
            for source, sed in zip(source_models, source_seds):
                smodel = 0. * datastack
                source.setPars()
                source.amp = 1.
                sed.setPars()
                spix = source.pixeval(xl, yl)
                mags = {}
                for i in range(nfitbands):
                    scale = sed.scale(fitbands[i], config['main_band'])
                    smodel[i*ny: (i+1)*ny, :] = scale * convolve.convolve(spix, convol_matrix[fitbands[i]], False)[0]
                    mags[fitbands[i]] = -2.5*np.log10(scale) + zp[fitbands[i]] - zp[config['main_band']]
                modlist.append((smodel/sigmastack).ravel()[maskstack_r])
                source_maglist.append(mags)
    
            modarr = np.array(modlist).T
    
            if np.isnan(modarr).any():
    
                return -1e300, fakemags
    
            amps, chi = nnls(modarr, (datastack/sigmastack).ravel()[maskstack_r])
    
            maglist = []
            i = 0
            for light, mags in zip(light_models, light_maglist):
                if amps[i] > 0.:
                    light.amp *= amps[i]
                    mainmag = light.Mag(zp[config['main_band']])
                    for band in fitbands:
                        mags[band] += mainmag
                else:
                    for band in fitbands:
                        mags[band] = 99.
                maglist.append(mags)
                i += 1
    
            for source, mags in zip(source_models, source_maglist):
                if amps[i] > 0.:
                    source.amp *= amps[i]
                    mainmag = source.Mag(zp[config['main_band']])
                    for band in fitbands:
                        mags[band] += mainmag
                else:
                    for band in fitbands:
                        mags[band] = 99.
                maglist.append(mags)
                i += 1
    
            logp = -0.5*chi**2
            if logp != logp:
                return -np.inf, fakemags
    
            return logp, maglist
    
        sampler = emcee.EnsembleSampler(nwalkers, npars, logpfunc)
    
        print "fitting model..."
    
        sampler.run_mcmc(start, config['Nsteps'])
    
        chain = sampler.chain
        magschain = sampler.blobs
    
        outchain = {}
        outchain['logp'] = sampler.lnprobability
    
        for i in range(npars):
            outchain[index2par[i]] = chain[:, :, i]
    
        for i in range(nlight):
            for band in fitbands:
                outchain['light%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, config['Nsteps']))
    
        for i in range(nsource):
            for band in fitbands:
                outchain['source%d.mag_%s'%(i+1, band)] = np.zeros((nwalkers, config['Nsteps']))
    
        for i in range(config['Nsteps']):
            for j in range(nwalkers):
                for band in fitbands:
                    for l in range(nlight):
                        outchain['light%d.mag_%s'%(l+1, band)][j, i] = magschain[i][j][l][band]
                    for s in range(nsource):
                        outchain['source%d.mag_%s'%(s+1, band)][j, i] = magschain[i][j][nlight+s][band]
        
        h5py_file = h5py.File(config['output_dir']+configfile+'_chain.hdf5', 'w')
        for par in outchain:
            h5py_file.create_dataset(par, data=outchain[par])
    
        save_model(config, h5py_file)
        write_config_file(config, h5py_file)
    
        h5py_file.close()
    
    elif key == 's':
    
        chain = h5py.File(config['output_dir']+configfile+'_chain.hdf5', 'r')
    
        save_model(config, chain, saveind)
        write_config_file(config, chain, saveind)


