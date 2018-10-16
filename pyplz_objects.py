import numpy as np
from imageSim import SBModels, SEDModels, convolve
import pyplz_rgbtools
from pylens import pylens, MassModels
from scipy.optimize import nnls
from astropy.io import fits as pyfits
import os


rootdir = os.environ.get('PYPLZDIR')

def read_config(filename):

    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    config = {'data_dir':'./', 'mask_dir': None, 'output_dir':'./', 'sps_model_dir': None, 'filters': None, 'main_band': None, \
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
                    parval = lines[i].split('#')[0].split(':')[1].split('\n')[0].strip()
                    config[parname] = parval

        i += 1

    filtlist = []
    filternames = config['filters'].split(',')
    for name in filternames:
        filtlist.append(name.strip())
    config['filters'] = filtlist

    config['colors'] = []
    for band in config['filters']:
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

    config['zeropoints'] = np.array(config['zeropoints'].split(','), dtype='float')
    config['Nsteps'] = int(config['Nsteps'])
    config['Nwalkers'] = int(config['Nwalkers'])
    if config['burnin'] is not None:
        config['burnin'] = int(config['burnin'])
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

                if model_class in SBModels.parlists:
                    parnames = []
                    for parname in SBModels.parlists[model_class]:
                        parnames.append(parname)
                else:
                    df

                sps_model = None
                if sed_class == 'freecolors':
                    parnames += config['colors']
                elif sed_class == 'template':
                    parnames += ['zs', 'tn']
                elif sed_class == 'sps':
                    sps_model = line[3]
                    parnames += ['redshift', 'age', 'tau', 'logZ', 'logtau_V']
                else:
                    df

                npars = len(parnames)

                comp = {'class':model_class, 'pars':{}, 'sed': sed_class}
                if sps_model is not None:
                    comp['sps_model'] = sps_model

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
                model_class = line[1].strip()
                sed_class = line[2]

                if model_class in ['Sersic']:
                    parnames = []
                    for parname in SBModels.parlists[model_class]:
                        parnames.append(parname)
                else:
                    df

                if sed_class == 'freecolors':
                    parnames += config['colors']
                elif sed_class == 'template':
                    parnames += ['zs', 'tn']
                elif sed_class == 'sps':
                    parnames += ['redshift', 'age', 'logZ', 'tau', 'logtau_V']
                else:
                    df

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
                model_class = line[1].strip()

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


class Par:

    def __init__(self, name, lower=0., upper=1., value=0., step=1.):

        self.name = name
        self.lower = lower
        self.upper = upper
        self.value = value
        self.step = step

class PyPLZModel:

    def __init__(self, config):

        self.pars = []
        self.par2index = {}
        self.index2par = {}    

        self.sci = {}
        self.err = {}
        self.convol_matrix = {}
        self.zp = {}

        self.light_sb_models = []
        self.source_sb_models = []
        self.light_sed_models = []
        self.source_sed_models = []

        self.light_mags = []
        self.source_mags = []

        self.lens_models = []

        self.logp = None

        # defines bands 

        self.bands = []
        self.main_band = config['main_band']

        for band in config['filters']:
            self.bands.append(band)

        self.nbands = len(self.bands)

        # reads in data
    
        filtnames = {}
        psfdic = {}
        for i in range(self.nbands):
            band = self.bands[i]
        
            self.zp[band] = config['zeropoints'][i]
        
            filtname = config['filter_prefix'] + band + config['filter_suffix']
        
            filtnames[band] = filtname

            hdu = pyfits.open(os.path.expandvars(config['data_dir'])+'/'+config['filename']+'_%s'%band+config['science_tag'])[0]
        
            img = hdu.data.copy()
            subimg = img.copy()
            self.sci[band] = subimg
            suberr = pyfits.open(os.path.expandvars(config['data_dir'])+'/'+config['filename']+'_%s'%band+config['err_tag'])[0].data.copy()
        
            if config['err_type'] == 'VAR':
                self.err[band] = suberr**0.5
            elif config['err_type'] == 'SIGMA':
                self.err[band] = suberr.copy()
            else:
                df
        
            psf_file = pyfits.open(os.path.expandvars(config['data_dir'])+'/'+config['filename']+'_%s'%band+config['psf_tag'])
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
        
            psfdic[band] = psf
            self.convol_matrix[band] = convolve.convolve(self.sci[band], psf)[1]

        # prepares mask and data array

        ny, nx = self.sci[self.bands[0]].shape
        X, Y = np.meshgrid(np.arange(1.*nx), np.arange(1.*ny))

        self.ny = ny
        self.nx = nx
        self.X = X
        self.Y = Y
        self.R = ((X - nx/2)**2 + (Y - ny/2)**2)**0.5
        
        if config['mask_dir'] is None:
            mask_dir = './'
        else:
            mask_dir = os.path.expandvars(config['mask_dir'])
        
        if config['maskname'] is not None:
            MASK = pyfits.open(mask_dir+config['maskname'])[0].data.copy()
        else:
            MASK = np.ones(X.shape, dtype=int)
        
        if config['rmax'] is not None:
            MASK[self.R > config['rmax']] = 0
        
        mask = MASK > 0
        mask_r = mask.ravel()
        
        self.modelstack = np.zeros((self.nbands* ny, nx))
        self.scistack = 0. * self.modelstack
        self.errstack = 0. * self.modelstack
        self.maskstack = np.zeros((self.nbands * ny, nx), dtype=bool)

        for i in range(self.nbands):
            self.scistack[i*ny: (i+1)*ny, :] = self.sci[self.bands[i]]
            self.errstack[i*ny: (i+1)*ny, :] = self.err[self.bands[i]]
            self.maskstack[i*ny: (i+1)*ny, :] = mask

        self.maskstack_r = self.maskstack.ravel()
  
        #defines model parameters
        
        ncomp = 0
        npar = 0
        for comp in config['light_components']:
            ncomp += 1
            for par in comp['pars']:
                parpar = comp['pars'][par]
                if parpar['link'] is None and parpar['var'] == 1:
                    self.pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value'], step=parpar['step']))
                    self.par2index['light'+str(ncomp)+'.'+par] = npar
                    self.index2par[npar] = 'light'+str(ncomp)+'.'+par
                    npar += 1
        
        ncomp = 0
        for comp in config['source_components']:
            ncomp += 1
            for par in comp['pars']:
                parpar = comp['pars'][par]
                if parpar['link'] is None and parpar['var'] == 1:
                    self.pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value'], step=parpar['step']))
                    self.par2index['source'+str(ncomp)+'.'+par] = npar
                    self.index2par[npar] = 'source'+str(ncomp)+'.'+par
                    npar += 1
        
        ncomp = 0
        for comp in config['lens_components']:
            ncomp += 1
            for par in comp['pars']:
                parpar = comp['pars'][par]
                if parpar['link'] is None and parpar['var'] == 1:
                    self.pars.append(Par(par+str(ncomp), lower=parpar['low'], upper=parpar['up'], value=parpar['value'], step=parpar['step']))
                    self.par2index['lens'+str(ncomp)+'.'+par] = npar
                    self.index2par[npar] = 'lens'+str(ncomp)+'.'+par
                    npar += 1
        
        i = 0
        
        filtnames = {}
        for band in config['filters']:
        
            filtname = config['filter_prefix'] + band + config['filter_suffix']
        
            filtnames[band] = filtname
        
        ncomp = 1
        for comp in config['lens_components']:
        
            name = 'lens%d'%ncomp
        
            pars_here = {}
            for par in comp['pars']:
                if comp['pars'][par]['link'] is None:
                    if comp['pars'][par]['var'] == 1:
                        pars_here[par] = self.pars[self.par2index[name+'.'+par]]
                    elif comp['pars'][par]['var'] == 0:
                        pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                    else:
                        df
                else:
                    pars_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]
        
            ncomp += 1
        
            lens = MassModels.PowerLaw(name, pars_here)
            self.lens_models.append(lens)
        
        self.nlens = len(self.lens_models)

        ncomp = 1

        for comp in config['light_components']:
        
            name = 'light%d'%ncomp
        
            pars_here = {}
            sed_here = {}
            for par in comp['pars']:
                if par in SBModels.parlists[comp['class']]:
                    if comp['pars'][par]['link'] is None:
                        if comp['pars'][par]['var'] == 1:
                            pars_here[par] = self.pars[self.par2index[name+'.'+par]]
                        else:
                            pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                    else:
                        pars_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]
                else:
                    if comp['pars'][par]['link'] is None:
                        if comp['pars'][par]['var'] == 1:
                            sed_here[par] = self.pars[self.par2index[name+'.'+par]]
                        else:
                            sed_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                    else:
                        sed_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]
        
            if comp['class'] == 'Sersic':
                light = SBModels.Sersic('light%d'%ncomp, pars_here)
            elif comp['class'] == 'PointSource':
                light = SBModels.PointSource('light%d'%ncomp, psfdic, pars_here)
            else:
                df

            self.light_sb_models.append(light)

            if comp['sed'] == 'freecolors':
                sed = SEDModels.Colors('light_sed%d'%ncomp, sed_here, self.bands, self.main_band, self.zp)
            elif comp['sed'] == 'template':
                sed = SEDModels.Template('light_sed%d'%ncomp, sed_here, filtnames, self.zp)
            elif comp['sed'] == 'sps':
                modelname = config['sps_model_dir'] + comp['sps_model']
                sed = SEDModels.SPS('light_sed%d'%ncomp, sed_here, self.bands, self.zp, modelname)
            else:
                df
            self.light_sed_models.append(sed)
            self.light_mags.append({})

            ncomp += 1
         
        self.nlight = len(self.light_sb_models)

        ncomp = 1
        for comp in config['source_components']:
        
            name = 'source%d'%ncomp
        
            pars_here = {}
            sed_here = {}
            for par in comp['pars']:
                if par in SBModels.parlists['Sersic']:
                    if comp['pars'][par]['link'] is None:
                        if comp['pars'][par]['var'] == 1:
                            pars_here[par] = self.pars[self.par2index[name+'.'+par]]
                        else:
                            pars_here[par] = comp['pars'][par]['value']
                    else:
                        pars_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]
                else:
                    if comp['pars'][par]['link'] is None:
                        if comp['pars'][par]['var'] == 1:
                            sed_here[par] = self.pars[self.par2index[name+'.'+par]]
                        else:
                            sed_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                    else:
                        sed_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]

            if comp['class'] == 'Sersic':
                source = SBModels.Sersic('source%d'%ncomp, pars_here)
            else:
                df

            self.source_sb_models.append(source)
            if comp['sed'] == 'freecolors':
                sed = SEDModels.Colors('light_sed%d'%ncomp, sed_here, self.bands, self.main_band, zp=self.zp)
            elif comp['sed'] == 'template':
                sed = SEDModels.Template('light_sed%d'%ncomp, sed_here, filtnames, zp=self.zp)
            self.source_sed_models.append(sed)
            self.source_mags.append({})
         
            ncomp += 1

        self.nsource = len(self.source_sb_models)

    def update(self):

        for lens in self.lens_models:
            lens.setPars()

        for light_sb in self.light_sb_models:
            light_sb.setPars()

        for light_sed in self.light_sed_models:
            light_sed.setPars()

        for source_sb in self.source_sb_models:
            source_sb.setPars()

        for source_sed in self.source_sed_models:
            source_sed.setPars()

    def optimize_amp(self):

        light_ind_model = []
        source_ind_model = []
        
        xl, yl = pylens.getDeflections(self.lens_models, (self.X, self.Y))
        
        modlist = []

        for light, sed, mags in zip(self.light_sb_models, self.light_sed_models, self.light_mags):
            lmodel = 0. * self.scistack
            light.amp = 1.
            if light.__class__.__name__ == 'PointSource':
                for i in range(self.nbands):
                    scale = sed.scale(self.bands[i], self.main_band)
                    lmodel[i*self.ny: (i+1)*self.ny, :] = light.pixeval(self.X, self.Y, self.bands[i])
                    mags[self.bands[i]] = -2.5*np.log10(scale) + self.zp[self.bands[i]] - self.zp[self.main_band]
            else:
                lpix = light.pixeval(self.X, self.Y)
                for i in range(self.nbands):
                    scale = sed.scale(self.bands[i], self.main_band)
                    lmodel[i*self.ny: (i+1)*self.ny, :] = scale * convolve.convolve(lpix, self.convol_matrix[self.bands[i]], False)[0]
                    mags[self.bands[i]] = -2.5*np.log10(scale) + self.zp[self.bands[i]] - self.zp[self.main_band]
            modlist.append((lmodel/self.errstack).ravel()[self.maskstack_r])

        for source, sed, mags in zip(self.source_sb_models, self.source_sed_models, self.source_mags):
            smodel = 0. * self.scistack
            source.amp = 1.
            spix = source.pixeval(xl, yl)
            for i in range(self.nbands):
                scale = sed.scale(self.bands[i], self.main_band)
                smodel[i*self.ny: (i+1)*self.ny, :] = scale * convolve.convolve(spix, self.convol_matrix[self.bands[i]], False)[0]
                mags[self.bands[i]] = -2.5*np.log10(scale) + self.zp[self.bands[i]] - self.zp[self.main_band]
            modlist.append((smodel/self.errstack).ravel()[self.maskstack_r])
        
        modarr = np.array(modlist).T
        if np.isnan(modarr).any() or not np.isfinite(modarr).any():
            amps = np.ones(self.nlight + self.nsource)
            chi = 1e300
        else:
            amps, chi = nnls(modarr, (self.scistack/self.errstack).ravel()[self.maskstack_r])

        i = 0
        for light, sed, mags in zip(self.light_sb_models, self.light_sed_models, self.light_mags):
                
            if amps[i] > 0.:
                light.amp *= amps[i]
                mainmag = light.Mag(self.zp[self.main_band])
                for band in self.bands:
                    mags[band] += mainmag
                if sed.__class__.__name__ == 'SPS':
                    mags['mstar'] = 10.**(-2./5.*(mainmag - sed.mags[self.main_band]))
            else:
                for band in self.bands:
                    mags[band] = 99.
                if sed.__class__.__name__ == 'SPS':
                    mags['mstar'] = 0.
            i += 1

        for source, mags in zip(self.source_sb_models, self.source_mags):
            if amps[i] > 0.:
                source.amp *= amps[i]
                mainmag = source.Mag(self.zp[self.main_band])
                for band in self.bands:
                    mags[band] = mags[band] + mainmag
            else:
                for band in self.bands:
                    mags[band] = 99.
            i += 1

        self.logp = -0.5*chi**2
        return chi**2

    def get_restUV_mags(self):

        light_mags = []
        for light, sed in zip(self.light_sb_models, self.light_sed_models):
            if sed.__class__.__name__ == 'Template':
                mainmag = light.Mag(self.zp[self.main_band])
                scale = sed.restUV_scale(self.main_band)
                light_mags.append(mainmag - 2.5*np.log10(scale))
            else:
                light_mags.append(None)

        source_mags = []
        for source, sed in zip(self.source_sb_models, self.source_sed_models):
            if sed.__class__.__name__ == 'Template':
                mainmag = source.Mag(self.zp[self.main_band])
                scale = sed.restUV_scale(self.main_band)
                source_mags.append(mainmag - 2.5*np.log10(scale))
            else:
                source_mags.append(None)

        return light_mags, source_mags

    def save(self, outname, make_rgb=True):
    
        fitsname = outname+'.fits'
        rgbname = outname+'_rgb.png'
   
        hdr = pyfits.Header()
        hdr['logp'] = self.logp
    
        light_ind_model = []
        source_ind_model = []

        xl, yl = pylens.getDeflections(self.lens_models, (self.X, self.Y))
 
        n = 0
        for light, sed, mags in zip(self.light_sb_models, self.light_sed_models, self.light_mags):
        
            light_ind_dic = {}
    
            if light.__class__.__name__ == 'PointSource':
                for band in self.bands:
                    hdr['%s.mag_%s'%(light.name, band)] = mags[band]
                    scale = sed.scale(band, self.main_band)
                    light_ind_dic[band] = scale * light.pixeval(self.X, self.Y, band) 
            else:
                lpix = light.pixeval(self.X, self.Y)
                for band in self.bands:
                    hdr['%s.mag_%s'%(light.name, band)] = mags[band]
                    scale = sed.scale(band, self.main_band)
                    light_ind_dic[band] = scale * convolve.convolve(lpix, self.convol_matrix[band], False)[0]
    
            light_ind_model.append(light_ind_dic)
    
            n += 1

        for source, sed, mags in zip(self.source_sb_models, self.source_sed_models, self.source_mags):
        
            source_ind_dic = {}
    
            spix = source.pixeval(xl, yl)
        
            for band in self.bands:
                hdr['%s.mag_%s'%(source.name, band)] = mags[band]
                scale = sed.scale(band, self.main_band)
                source_ind_here = scale * convolve.convolve(spix, self.convol_matrix[band], False)[0]
    
                source_ind_dic[band] = source_ind_here
    
            source_ind_model.append(source_ind_dic)
        
            n += 1

        phdu = pyfits.PrimaryHDU(header=hdr)
    
        hdulist = pyfits.HDUList([phdu])
    
        for light, light_ind_dic in zip(self.light_sb_models, light_ind_model):
            
            for band in self.bands:
                hdu_here = pyfits.ImageHDU(data=light_ind_dic[band])
                hdu_here.header['EXTNAME'] = '%s_%s'%(light.name, band)
    
                hdulist.append(hdu_here)
    
        for source, source_ind_dic in zip(self.source_sb_models, source_ind_model):
            for band in self.bands:
                hdu_here = pyfits.ImageHDU(data=source_ind_dic[band])
                hdu_here.header['EXTNAME'] = '%s_%s'%(source.name, band)
    
                hdulist.append(hdu_here)
    
        hdulist.writeto(fitsname, overwrite=True)
        
        if make_rgb:

            # makes model rgb
            sci_list = []
            light_list = []
            source_list = []
            for band in self.bands:
                sci_list.append(self.sci[band])
                lmodel = 0.*self.sci[band]
                smodel = 0.*self.sci[band]
                for light in light_ind_model:
                    lmodel += light[band]
                light_list.append(lmodel)
                for source in source_ind_model:
                    smodel += source[band]
                source_list.append(smodel)
            
            pyplz_rgbtools.make_full_rgb(sci_list, light_list, source_list, outname=rgbname)

    def write_config_file(self, config, outname):
    
        conflines = []
        confpars = ['data_dir', 'mask_dir', 'maskname', 'output_dir', 'sps_model_dir', 'filename', 'science_tag', 'err_tag', 'err_type', 'psf_tag', 'rmax', 'Nwalkers', 'Nsteps', 'burnin', 'main_band', 'filter_prefix', 'filter_suffix']
        for parname in confpars:
            if config[parname] is not None:
                conflines.append('%s: %s\n'%(parname, config[parname]))
        filtline = 'filters: '
        zpline = 'zeropoints: '
        for i in range(self.nbands -1):
            filtline += '%s, '%self.bands[i]
            zpline += '%f, '%self.zp[self.bands[i]]
        filtline += self.bands[-1]+'\n'
        zpline += '%f\n'%self.zp[self.bands[-1]]
        conflines.append(filtline)
        conflines.append(zpline)
        
        conflines.append('\n')
        conflines.append('# MODELS\n')

        light_uvmags, source_uvmags = self.get_restUV_mags()

        ncomp = 0
        for comp, mags, uvmag in zip(config['light_components'], self.light_mags, light_uvmags):
    
            if comp['sed'] == 'freecolors':
                lightpars = SBModels.parlists[comp['class']] + config['colors']
            elif comp['sed'] == 'template':
                lightpars = SBModels.parlists[comp['class']] + ['zs', 'tn']
            elif comp['sed'] == 'sps':
                lightpars = SBModels.parlists[comp['class']] + ['redshift', 'age', 'logZ', 'tau', 'logtau_V']

            conflines.append('\n')
            modeltypeline = 'light_model %s %s'%(comp['class'], comp['sed'])
            if 'sps_model' in comp:
                modeltypeline += ' %s\n'%comp['sps_model']
            else:
                modeltypeline += '\n'
            conflines.append(modeltypeline)
            for par in lightpars:
                parname = 'light%d.%s'%(ncomp+1, par)
                if parname in self.par2index:
                    npar = self.par2index[parname]
                    conflines.append('%s %f %f %f %f 1\n'%(par, self.pars[npar].value, self.pars[npar].lower, self.pars[npar].upper, self.pars[npar].step))
                else:
                    if config['light_components'][ncomp]['pars'][par]['link'] is None:
                        conflines.append('%s %f -1 -1 -1 0\n'%(par, config['light_components'][ncomp]['pars'][par]['value']))
                    else:
                        lname = config['light_components'][ncomp]['pars'][par]['link']
                        npar = self.par2index[lname]
                        conflines.append('%s %f %f %f %f 1 %s\n'%(par, self.pars[npar].value, self.pars[npar].lower, self.pars[npar].upper, self.pars[npar].step, lname))
            for band in self.bands:
                conflines.append('mag_%s %3.2f\n'%(band, mags[band]))
            if uvmag is not None:
                conflines.append('uvmag %3.2f\n'%uvmag)
            if 'mstar' in mags:
                conflines.append('mstar %4.3f\n'%mags['mstar'])
            ncomp += 1
        
        ncomp = 0
    
        for comp, mags, uvmag in zip(config['source_components'], self.source_mags, source_uvmags):
    
            if comp['sed'] == 'freecolors':
                sourcepars = SBModels.parlists['Sersic'] + config['colors']
            elif comp['sed'] == 'template':
                sourcepars = SBModels.parlists['Sersic'] + ['zs', 'tn']
    
            conflines.append('\n')
            conflines.append('source_model Sersic %s\n'%comp['sed'])
            for par in sourcepars:
                parname = 'source%d.%s'%(ncomp+1, par)
                if parname in self.par2index:
                    npar = self.par2index[parname]
                    conflines.append('%s %f %f %f %f 1\n'%(par, self.pars[npar].value, self.pars[npar].lower, self.pars[npar].upper, self.pars[npar].step))
                else:
                    if config['source_components'][ncomp]['pars'][par]['link'] is None:
                        conflines.append('%s %f -1 -1 -1 0\n'%(par, config['source_components'][ncomp]['pars'][par]['value']))
                    else:
                        lname = config['source_components'][ncomp]['pars'][par]['link']
                        npar = self.par2index[lname]
                        conflines.append('%s %f %f %f %f 1 %s\n'%(par, self.pars[npar].value, self.pars[npar].lower, self.pars[npar].upper, self.pars[npar].step, lname))
            for band in self.bands:
                conflines.append('mag_%s %3.2f\n'%(band, mags[band]))
            if uvmag is not None:
                conflines.append('uvmag %3.2f\n'%uvmag)
    
            ncomp += 1
        
        ncomp = 0
        powpars = ['x', 'y', 'pa', 'q', 'b', 'eta']
        
        for lens in self.lens_models:
            conflines.append('\n')
            conflines.append('lens_model Powerlaw\n')
            for par in powpars:
                parname = 'lens%d.%s'%(ncomp+1, par)
                if parname in self.par2index:
                    npar = self.par2index[parname]
                    conflines.append('%s %f %f %f %f 1\n'%(par, self.pars[npar].value, self.pars[npar].lower, self.pars[npar].upper, self.pars[npar].step))
                else:
                    if config['lens_components'][ncomp]['pars'][par]['link'] is None:
                        conflines.append('%s %f -1 -1 -1 0\n'%(par, config['lens_components'][ncomp]['pars'][par]['value']))
                    else:
                        lname = config['lens_components'][ncomp]['pars'][par]['link']
                        npar = self.par2index[lname]
                        conflines.append('%s %f %f %f %f 1 %s\n'%(par, self.pars[npar].value, self.pars[npar].lower, self.pars[npar].upper, self.pars[npar].step, lname))
            ncomp += 1
        
        conflines.append('\n')
        conflines.append('logp %f\n'%self.logp)

        f = open(outname, 'w')
        f.writelines(conflines)
        f.close()
 
