import numpy as np
from imageSimg import SBModels, SEDModels


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
        self.psf = {}

        self.light_sb_models = []
        self.source_sb_models = []
        self.light_sed_models = []
        self.source_sed_models = []

        self.lens_models = []

        #defines model parameters
        
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
        
        filtdic = {}
        for band in config['filters']:
        
            zp[band] = config['zeropoints'][i]
        
            filtname = rootdir+'pylens/filters/%s%s%s'%(config['filter_prefix'], band, config['filter_suffix'])
        
            f = open(filtname, 'r')
            ftable = np.loadtxt(f)
            f.close()
        
            filtdic[band] = (ftable[:, 0], ftable[:, 1])
        
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
        
        ncomp = 1

        for comp in config['light_components']:
        
            name = 'light%d'%ncomp
        
            pars_here = {}
            sed_here = {}
            for par in comp['pars']:
                if par in SBModels.sersicpars:
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
        
            light = SBModels.Sersic('light%d'%ncomp, pars_here)
            self.light_sb_models.append(light)

            if comp['sed'] == 'freecolors':
                sed = SEDModels.Colors('light_sed%d'%ncomp, sed_here, zp)
            elif comp['sed'] == 'template':
                sed = SEDModels.Template('light_sed%d'%ncomp, sed_here, filtdic, zp)
            self.light_sed_models.append(sed)

            ncomp += 1
         
        ncomp = 1
        for comp in config['source_components']:
        
            name = 'source%d'%ncomp
        
            pars_here = {}
            for par in comp['pars']:
                if par in SBModels.sersicpars:
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

            source = SBModels.Sersic('source%d'%ncomp, pars_here)
            self.source_sb_models.append(source)
            if comp['sed'] == 'freecolors':
                sed = SEDModels.Colors('light_sed%d'%ncomp, sed_here, filtdic, zp=zp)
            elif comp['sed'] == 'template':
                sed = SEDModels.Template('light_sed%d'%ncomp, sed_here, filtdic, zp=zp)
            self.source_sed_models.append(sed)
         
            ncomp += 1

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

    def optimize_amp():

