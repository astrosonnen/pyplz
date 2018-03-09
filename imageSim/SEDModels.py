import numpy as np
from scipy.interpolate import splrep, splev
import os


tempdir = os.environ.get('PYPLZDIR') + '/templates/'

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

    def __init__(self, name, pars, zp):

        self.keys = []
        for par in pars:
            self.keys.append(par)

        self.vmap = {}
        self.pars = pars
        for key in self.keys:
            self.vmap[key] = self.pars[key]

        self.name = name
        self.zp = zp

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key, self.vmap[key].value)

    def scale(self, band2, band1):
        if band2 == band1:
            return 1.
        else:
            col = self.__dict__['%s-%s'%(band2, band1)]
            return 10.**(-(col - self.zp[band2] + self.zp[band1])/2.5)


class Template:

    def __init__(self, name, pars, filters, zp, normrange=(4000., 5000.)):

        splinedic = {}
        rangedic = {}
    
        wavmax = 0.
        for band in filters:
            splinedic[band] = splrep(filters[band][0], filters[band][1])
            rangedic[band] = (filters[band][0][0], filters[band][0][-1])
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
            self.vmap[key] = self.pars[key]

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

        temp_fnus = []
        for n in temp_ind:

            flambda = self.templates[n][1].copy()
            wav = self.templates[n][0].copy()

            wobs = wav * (1. + zs)

            fnus_here = []
            for band in bands_here:
                wrange = (wobs > self.filtranges[band][0]) & (wobs < self.filtranges[band][1])
                wband = wobs[wrange]
                fband = flambda[wrange]

                if madau:
                    madau_corr = etau_madau(wband, zs)
                    fband *= madau_corr

                fnu = fband * wband**2
                weights = splev(wband, self.filtsplines[band])
                fnus_here.append((fnu * weights).sum()/weights.sum())

            temp_fnus.append(fnus_here)

        ratio = (a1 * temp_fnus[0][0] + a2 * temp_fnus[1][0]) / (a1 * temp_fnus[0][1] + a2 * temp_fnus[1][1])

        return ratio * 10.**((self.zp[band2] - self.zp[band1])/2.5)


