import numpy as np
from pyplz_configtools import read_config
from imageSim import SBModels, SEDModels, convolve
import pyplz_rgbtools
from pylens import deflections as pylens, MassModels
from scipy.optimize import nnls
from astropy.io import fits as pyfits
import os


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

        self.light_mags = []
        self.source_mags = []

        self.lens_models = []

        self.logp = None

        # defines bands 

        self.bands = []

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
        
        self.mask = MASK > 0
        self.mask_r = self.mask.ravel()

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
            for par in comp['pars']:
                if par in SBModels.parlists[comp['class']]:
                    if comp['pars'][par]['link'] is None:
                        if comp['pars'][par]['var'] == 1:
                            pars_here[par] = self.pars[self.par2index[name+'.'+par]]
                        else:
                            pars_here[par] = Par(par, lower=comp['pars'][par]['value'], upper=comp['pars'][par]['value'], value=comp['pars'][par]['value'])
                    else:
                        pars_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]
       
            if comp['class'] == 'Sersic':
                light = SBModels.Sersic('light%d'%ncomp, pars_here)
            elif comp['class'] == 'PointSource':
                light = SBModels.PointSource('light%d'%ncomp, psfdic, pars_here)
            elif comp['class'] == 'Ring':
                light = SBModels.Ring('light%d'%ncomp, pars_here)
            else:
                df

            self.light_sb_models.append(light)
            self.light_mags.append({})

            ncomp += 1
         
        self.nlight = len(self.light_sb_models)

        ncomp = 1
        for comp in config['source_components']:
        
            name = 'source%d'%ncomp
        
            pars_here = {}
            for par in comp['pars']:
                if par in SBModels.parlists['Sersic']:
                    if comp['pars'][par]['link'] is None:
                        if comp['pars'][par]['var'] == 1:
                            pars_here[par] = self.pars[self.par2index[name+'.'+par]]
                        else:
                            pars_here[par] = comp['pars'][par]['value']
                    else:
                        pars_here[par] = self.pars[self.par2index[comp['pars'][par]['link']]]

            if comp['class'] == 'Sersic':
                source = SBModels.Sersic('source%d'%ncomp, pars_here)
            else:
                df

            self.source_sb_models.append(source)
            self.source_mags.append({})
         
            ncomp += 1

        self.nsource = len(self.source_sb_models)

    def update(self):

        for lens in self.lens_models:
            lens.setPars()

        for light_sb in self.light_sb_models:
            light_sb.setPars()

        for source_sb in self.source_sb_models:
            source_sb.setPars()

    def optimize_amp(self):

        light_ind_model = []
        source_ind_model = []
        
        xl, yl = pylens.getDeflections(self.lens_models, (self.X, self.Y))

        lpix = []
        spix = []

        for light in self.light_sb_models:
            light.amp = 1.
            lpix.append(light.pixeval(self.X, self.Y))

        for source in self.source_sb_models:
            source.amp = 1.
            spix.append(source.pixeval(xl, yl))

        chi2sum = 0.

        for band in self.bands:
            modlist = []
            for l in lpix:
                modlist.append((convolve.convolve(l, self.convol_matrix[band], False)[0]/self.err[band]).ravel()[self.mask_r])
            for s in spix:
                modlist.append((convolve.convolve(s, self.convol_matrix[band], False)[0]/self.err[band]).ravel()[self.mask_r])

            modarr = np.array(modlist).T
            if np.isnan(modarr).any() or not np.isfinite(modarr).any():
                self.logp = -0.5e300
                return 1e300
            else:
                amps, chi = nnls(modarr, (self.sci[band]/self.err[band]).ravel()[self.mask_r])
                chi2sum += chi**2

            i = 0
            for light, mags in zip(self.light_sb_models, self.light_mags):
                if amps[i] > 0.:
                    light.amp = amps[i]
                    mags[band] = light.Mag(self.zp[band])
                else:
                    mags[band] = 99.
                i += 1

            for source, mags in zip(self.source_sb_models, self.source_mags):
                if amps[i] > 0.:
                    source.amp = amps[i]
                    mags[band] = source.Mag(self.zp[band])
                else:
                    mags[band] = 99.
                i += 1

        self.logp = -0.5*chi2sum
        return chi2sum

    def save(self, outname, config, make_rgb=True):
    
        fitsname = outname+'.fits'
        rgbname = outname+'_rgb.png'
   
        hdr = pyfits.Header()
        hdr['logp'] = self.logp
    
        light_ind_model = []
        source_ind_model = []

        xl, yl = pylens.getDeflections(self.lens_models, (self.X, self.Y))
 
        n = 0
        for light, mags in zip(self.light_sb_models, self.light_mags):
        
            light_ind_dic = {}
    
            light.amp = 1.
            lpix = light.pixeval(self.X, self.Y)
            for band in self.bands:
                hdr['%s.mag_%s'%(light.name, band)] = mags[band]
                scale = 10.**(-2./5.*(mags[band] - light.Mag(self.zp[band])))

                light_ind_dic[band] = scale * convolve.convolve(lpix, self.convol_matrix[band], False)[0]
    
            light_ind_model.append(light_ind_dic)
    
            n += 1

        for source, mags in zip(self.source_sb_models, self.source_mags):

            source_ind_dic = {}
    
            source.amp = 1.
            spix = source.pixeval(xl, yl)
        
            for band in self.bands:
                hdr['%s.mag_%s'%(source.name, band)] = mags[band]
                scale = 10.**(-2./5.*(mags[band] - source.Mag(self.zp[band])))
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
            model_list = []
            for band in self.bands:
                sci_list.append(self.sci[band])
                comp_list = []
                for light in light_ind_model:
                    comp_list.append(light[band])
                for source in source_ind_model:
                    comp_list.append(source[band])
                model_list.append(comp_list)
            
            pyplz_rgbtools.make_full_rgb(sci_list, model_list, outname=rgbname, scheme=config['rgbscheme'], cuts=config['rgbcuts'], scales=config['rgbscales'])

    def write_config_file(self, config, outname):
    
        conflines = []
        confpars = ['data_dir', 'mask_dir', 'maskname', 'output_dir', 'filename', 'science_tag', 'err_tag', 'err_type', 'psf_tag', 'rmax', 'Nwalkers', 'Nsteps', 'burnin', 'filter_prefix', 'filter_suffix', 'modeltype', 'rgbscheme']
        for parname in confpars:
            if config[parname] is not None:
                conflines.append('%s: %s\n'%(parname, config[parname]))
        if config['rgbcuts'] is not None:
            rgbcutsline = 'rgbcuts:'
            for cut in config['rgbcuts']:
                rgbcutsline += ' %2.1f'%cut
            rgbcutsline += '\n'
            conflines.append(rgbcutsline)

        if config['rgbscales'] is not None:
            rgbscalesline = 'rgbscales:'
            for scale in config['rgbscales']:
                rgbscalesline += ' %f'%scale
            rgbscalesline += '\n'
            conflines.append(rgbscalesline)

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

        ncomp = 0
        for comp, mags in zip(config['light_components'], self.light_mags):

            lightpars = SBModels.parlists[comp['class']]    

            conflines.append('\n')
            modeltypeline = 'light_model %s MLflux\n'%(comp['class'])
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
            ncomp += 1
        
        ncomp = 0
    
        for comp, mags in zip(config['source_components'], self.source_mags):
            sourcepars = SBModels.parlists[comp['class']]    
   
            conflines.append('\n')
            conflines.append('source_model Sersic MLflux\n')
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
 
