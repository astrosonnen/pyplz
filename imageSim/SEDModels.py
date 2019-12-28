import numpy as np
from scipy.interpolate import splrep, splev, splint
import os
import pickle
import pyplz_cosmology
import pyplz_dustlaw


tempdir = os.environ.get('PYPLZDIR') + '/templates/'
filtdir = os.environ.get('PYPLZDIR') + '/filters/'

parlists = {'template': ['tn', 'zs'], 'sps': ['redshift', 'age', 'tau', 'logZ', 'logtau_V']}

def extinct(wav, ebmv, Rv=3.1):
    return pyplz_dustlaw.Alambda(wav, Rv=Rv)*Rv*ebmv
 
def etau_madau(wl, z):
    """
    Madau 1995 extinction for a galaxy spectrum at redshift z 
    defined on a wavelength grid wl
    """
    n=len(wl)
    l=np.array([1216.,1026.,973.,950.])
    xe=1.+z
    
    #If all the spectrum is redder than (1+z)*wl_lyman_alfa 
    if wl[0]> l[0]*xe: return np.zeros(n)+1.
    
    #Madau coefficients
    c=np.array([3.6e-3,1.7e-3,1.2e-3,9.3e-4])
    ll=912.
    tau=wl*0.
    i1=np.searchsorted(wl,ll)
    i2=n-1
    #Lyman series absorption
    for i in range(len(l)):
        i2=np.searchsorted(wl[i1:i2],l[i]*xe)
        tau[i1:i2]=tau[i1:i2]+c[i]*(wl[i1:i2]/l[i])**3.46

    if ll*xe < wl[0]:
        return np.exp(-tau)

    #Photoelectric absorption
    xe=1.+z
    i2=np.searchsorted(wl,ll*xe)
    xc=wl[i1:i2]/ll
    xc3=xc**3
    tau[i1:i2]=tau[i1:i2]+\
                (0.25*xc3*(xe**.46-xc**0.46)\
                 +9.4*xc**1.5*(xe**0.18-xc**0.18)\
                 -0.7*xc3*(xc**(-1.32)-xe**(-1.32))\
                 -0.023*(xe**1.68-xc**1.68))

    tau = np.clip(tau, 0, 700)
    return np.exp(-tau)
    # if tau>700. : return 0.
    # else: return exp(-tau)


class Colors:

    def __init__(self, name, pars, bands, main_band, zp):

        self.keys = []
        if not 'main_mag' in pars:
            pars['main_mag'] = 0.

        for par in pars:
            self.keys.append(par)

        self.vmap = {}
        self.pars = pars
        for key in self.keys:
            if self.pars[key].__class__.__name__ == 'Par':
                self.vmap[key] = self.pars[key]
            else:
                self.__setattr__(key, self.pars[key])

        self.bands = bands
        self.main_band = main_band

        self.mags = {}
        for band in self.bands:
            self.mags[band] = None

        self.name = name
        self.zp = zp

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key, self.vmap[key].value)
        for band in self.bands:
            if band == self.main_band:
                self.mags[band] = self.main_mag
            else:
                self.mags[band] = self.__dict__['%s-%s'%(band, self.main_band)] + self.main_mag

    def scale(self, band2, band1):
        return 10.**(-(self.mags[band2] - self.mags[band1] - self.zp[band2] + self.zp[band1])/2.5)

    def restUV_scale(self, refband):
        return None


class Template:

    def __init__(self, name, pars, filtnames, zp, normrange=(4000., 5000.)):

        splinedic = {}
        rangedic = {}

        wavmax = 0.
        for band in filtnames:
        
            f = open(filtdir+filtnames[band], 'r')
            ftable = np.loadtxt(f)
            f.close()
        
            splinedic[band] = splrep(ftable[:, 0], ftable[:, 1])
            rangedic[band] = (ftable[0, 0], ftable[-1, 0])
            if rangedic[band][1] > wavmax:
                wavmax = rangedic[band][1]
    
        self.filtsplines = splinedic
        self.filtranges = rangedic
    
        f = open(tempdir+'CWWSB4.list', 'r')
        lines = f.readlines()
        f.close()
    
        templates = []
        
        for line in lines:
            line = line.rstrip()
            f = open(tempdir+line, 'r')
            tab = np.loadtxt(f)
            f.close()
    
            wav = tab[:, 0]
            flambda = tab[:, 1]
    
            wrange = wav < wavmax
    
            norm = 1./np.median(flambda[(wav > normrange[0]) & (wav < normrange[1])])
    
            templates.append((wav[wrange], norm * flambda[wrange]))
    
        self.templates = templates

        self.keys = ['zs', 'tn']

        self.vmap = {}
        self.pars = pars
        for key in self.keys:
            if self.pars[key].__class__.__name__ == 'Par':
                self.vmap[key] = self.pars[key]
            else:
                self.__setattr__(key, self.pars[key])

        self.name = name
        self.zp = zp

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key, self.vmap[key].value)

    def scale(self, band2, band1, madau=True):
        zs = self.zs

        a2 = self.tn%1
        a1 = 1. - a2
        n1 = int(np.floor(self.tn))
        n2 = int(np.ceil(self.tn))

        bands_here = [band2, band1]
        temp_coeffs = [a1, a2]
        temp_ind = [n1, n2]

        temp_fband = []
        for n in temp_ind:

            flambda = self.templates[n][1].copy()
            wav = self.templates[n][0].copy()

            wobs = wav * (1. + zs)

            fband_here = []

            for band in bands_here:
                wrange = (wobs > self.filtranges[band][0]) & (wobs < self.filtranges[band][1])
                wband = wobs[wrange]
                flambda_here = flambda[wrange]

                if madau:
                    madau_corr = etau_madau(wband, zs)
                    flambda_here *= madau_corr

                nu_here = 1./wband
                fnu_here = flambda_here/nu_here**2
                weights = splev(wband, self.filtsplines[band])

                nu_here = nu_here[::-1]
                fnu_here = fnu_here[::-1]
                weights = weights[::-1]

                num_integrand = splrep(nu_here, weights * fnu_here / nu_here, k=1)
                den_integrand = splrep(nu_here, weights / nu_here, k=1) 

                num = splint(nu_here[0], nu_here[-1], num_integrand)
                den = splint(nu_here[0], nu_here[-1], den_integrand)

                fband_here.append(num/den)

            temp_fband.append(fband_here)

        ratio = (a1 * temp_fband[0][0] + a2 * temp_fband[1][0]) / (a1 * temp_fband[0][1] + a2 * temp_fband[1][1])

        return ratio * 10.**((self.zp[band2] - self.zp[band1])/2.5)

    def restUV_scale(self, refband, uvrange=(1500., 1700.), madau=True):
        zs = self.zs

        a2 = self.tn%1
        a1 = 1. - a2
        n1 = int(np.floor(self.tn))
        n2 = int(np.ceil(self.tn))

        temp_coeffs = [a1, a2]
        temp_ind = [n1, n2]

        temp_fnus = []
        for n in temp_ind:

            flambda = self.templates[n][1].copy()
            wav = self.templates[n][0].copy()

            wobs = wav * (1. + zs)

            fnus_here = []

            # computes fnu at rest-frame UV

            wuvrange = (wav > uvrange[0]) & (wav < uvrange[1])
            wuvband = wobs[wuvrange]
            fuvband = flambda[wuvrange]

            if madau:
                madau_corr = etau_madau(wuvband, zs)
                fuvband *= madau_corr

            fuvnu = fuvband * wuvband**2
            fnus_here.append(np.median(fuvnu))

            # computes fnu in reference filter

            wrefrange = (wobs > self.filtranges[refband][0]) & (wobs < self.filtranges[refband][1])
            wrefband = wobs[wrefrange]
            frefband = flambda[wrefrange]

            if madau:
                madau_corr = etau_madau(wrefband, zs)
                frefband *= madau_corr

            frefnu = frefband * wrefband**2
            weights = splev(wrefband, self.filtsplines[refband])
            fnus_here.append((frefnu * weights).sum()/weights.sum())

            temp_fnus.append(fnus_here)

        ratio = (a1 * temp_fnus[0][0] + a2 * temp_fnus[1][0]) / (a1 * temp_fnus[0][1] + a2 * temp_fnus[1][1])

        return ratio 


class SPS:

    def __init__(self, name, pars, bands, zp, modelname, filtnames, ebmv=0., restmodelname=None):

        self.keys = ['logmass', 'redshift', 'age', 'tau', 'logZ', 'logtau_V']

        self.vmap = {}
        if not 'logmass' in pars:
            pars['logmass'] = 1.
        self.pars = pars

        for key in self.keys:
            if self.pars[key].__class__.__name__ == 'Par':
                self.vmap[key] = self.pars[key]
            else:
                self.__setattr__(key, self.pars[key])

        self.name = name
        self.zp = zp
        self.bands = bands


        self.mags = {}
        for band in self.bands:
            self.mags[band] = None

        self.extinction_corr = {}

        for band in self.bands:
        
            f = open(filtdir+filtnames[band], 'r')
            ftable = np.loadtxt(f)
            f.close()
        
            filtwav = (ftable[:, 0] * ftable[:, 1]).sum() / ftable[:, 1].sum()
            self.extinction_corr[band] = extinct(filtwav, ebmv)

        f = open(modelname, 'rb')
        self.model = pickle.load(f)
        f.close()

        if restmodelname is not None:
            f = open(restmodelname, 'r')
            self.restmodel = pickle.load(f)
            f.close()
        else:
            self.restmodel = None

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key, self.vmap[key].value)

        # calculates fluxes

        pnt = {}
        pnt['redshift'] = [self.redshift]
        pnt['age'] = [self.age]
        pnt['tau'] = [self.tau]
        pnt['Z'] = [10.**self.logZ]
        pnt['tau_V'] = [10.**self.logtau_V]

        for band in self.bands:
            self.mags[band] = self.model.eval(pnt, band)[0] - 2.5*self.logmass + self.extinction_corr[band]

    def scale(self, band2, band1):
        return 10.**(-2./5*(self.mags[band2] - self.zp[band2] - self.mags[band1] + self.zp[band1]))


