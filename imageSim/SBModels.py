from imageSim import SBProfiles as SBProfiles
from math import pi
import numpy as np
from scipy.interpolate import splrep, splev
from scipy import ndimage
import os


tempdir = os.environ.get('PYPLZDIR') + '/templates/'
parlists = {'Sersic': ['x', 'y', 'q', 're', 'pa', 'n'], 'Sersic_e1e2': ['x', 'y', 'e1', 'e2', 're', 'n'], 'PointSource': ['x', 'y'],\
            'Ring': ['x', 'y', 'q', 'rr', 'hi', 'ho', 'pa']}

def cnts2mag(cnts, zp):
    from math import log10
    return -2.5*log10(cnts) + zp


class SBModel:
    def __init__(self, name, pars, convolve=0):
        if 'amp' not in pars.keys() and 'logamp' not in pars.keys():
            pars['amp'] = 1.
        self.keys = pars.keys()
        keylist = sorted(self.keys)
        if keylist not in self._SBkeys:
            import sys
            print(keylist)
            print('Not all (or too many) parameters were defined!')
            sys.exit()
        self._baseProfile.__init__(self)
        self.vmap = {}
        self.pars = pars
        for key in keylist:
            try:
                v = self.pars[key].value
                self.vmap[key] = self.pars[key]
            except:
                self.__setattr__(key,self.pars[key])
        self.setPars()
        self.name = name
        self.convolve = convolve


    def __setattr__(self,key,value):
        if key=='pa':
            self.__dict__['pa'] = value
            if value is not None:
                self.__dict__['theta'] = value*pi/180.
        elif key=='theta':
            if value is not None:
                self.__dict__['pa'] = value*180./pi
            self.__dict__['theta'] = value
        elif key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value


    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

class PixelizedModel:
    def __init__(self, psfs):
        self.psfs = {}
        for band in psfs:
            self.psfs[band] = psfs[band].copy()
            self.psfs[band] /= psfs[band].sum()
        self.setCentroid()
        self.createModel()
        self.x = None
        self.y = None
        self.amp = 1.

    def setCentroid(self):
        self.x0 = {}
        self.y0 = {}
        for band in self.psfs:
            y, x = np.indices(self.psfs[band].shape).astype(np.float32)
            self.x0[band] = (x*self.psfs[band]).sum()
            self.y0[band] = (y*self.psfs[band]).sum()

    def createModel(self, order=1):
        self.model = {}
        for band in self.psfs:
            if order==1:
                self.model[band] = self.psfs[band].copy()
            else:
                self.model[band] = ndimage.spline_filter(self.psfs[band], output=np.float64, order=order)
        self.order = order

    def pixeval(self, x, y, band):
        X = x-self.x+self.x0[band]
        Y = y-self.y+self.y0[band]
        psf = ndimage.map_coordinates(self.model[band], [Y,X], prefilter=False)
        psf /= psf.sum()
        return self.amp*psf


class Sersic(SBModel,SBProfiles.Sersic):
    _baseProfile = SBProfiles.Sersic
    _SBkeys = [['amp','n','pa','q','re','x','y'],
                ['logamp','n','pa','q','re','x','y'],
                ['amp','n','q','re','theta','x','y'],
                ['logamp','n','q','re','theta','x','y'], \
                ['amp', 'e1', 'e2', 'n', 're', 'x', 'y']]

    def __init__(self,name,pars,convolve=0):
        SBModel.__init__(self,name,pars,convolve)

    def getMag(self,amp,zp):
        from scipy.special import gamma
        from math import exp,pi
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        cnts = (re**2)*amp*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)

    def setAmpFromMag(self,mag,zp):
        from math import exp,log10,pi
        from scipy.special import gamma
        cnts = 10**(-0.4*(mag-zp))
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        self.amp = cnts/((re**2)*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi)


class Sersic_e1e2(SBModel,SBProfiles.Sersic_e1e2):
    _baseProfile = SBProfiles.Sersic_e1e2
    _SBkeys = [['amp', 'e1', 'e2', 'n', 're', 'x', 'y']]

    def __init__(self,name,pars,convolve=0):
        SBModel.__init__(self,name,pars,convolve)

    def getMag(self,amp,zp):
        from scipy.special import gamma
        from math import exp,pi
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        cnts = (re**2)*amp*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi
        return cnts2mag(cnts,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)

    def setAmpFromMag(self,mag,zp):
        from math import exp,log10,pi
        from scipy.special import gamma
        cnts = 10**(-0.4*(mag-zp))
        n = self.n
        re = self.re
        k = 2.*n-1./3+4./(405.*n)+46/(25515.*n**2)
        self.amp = cnts/((re**2)*exp(k)*n*(k**(-2*n))*gamma(2*n)*2*pi)


class PointSource(PixelizedModel):
    def __init__(self, name, model, pars):
        if 'amp' not in pars.keys():
            pars['amp'] = 1.
        self.keys = pars.keys()
        keylist = sorted(self.keys)
        if keylist!=['amp', 'x', 'y']:
            import sys
            print('Not all (or too many) parameters were defined!')
            print(self.keys)
            sys.exit()
        PixelizedModel.__init__(self, model)
        self.vmap = {}
        self.pars = pars
        for key in keylist:
            try:
                v = self.pars[key].value
                self.vmap[key] = self.pars[key]
            except:
                self.__setattr__(key,self.pars[key])
        self.setPars()
        self.name = name

    def __setattr__(self,key,value):
        if key=='logamp':
            if value is not None:
                self.__dict__['amp'] = 10**value
        else:
            self.__dict__[key] = value

    def setPars(self):
        for key in self.vmap:
            self.__setattr__(key,self.vmap[key].value)

    def getMag(self,amp,zp):
        return cnts2mag(amp,zp)

    def Mag(self,zp):
        return self.getMag(self.amp,zp)


class Ring(SBModel, SBProfiles.Ring):
    _baseProfile = SBProfiles.Ring
    _SBkeys = [['amp', 'hi', 'ho', 'pa', 'q', 'rr', 'x', 'y']]

    def __init__(self, name, pars, convolve=0):
        SBModel.__init__(self, name, pars, convolve)

    def getMag(self, amp, zp):
        cnts = amp
        return cnts2mag(cnts, zp)

    def Mag(self, zp):
        return self.getMag(self.amp, zp)


class Spiral(SBModel, SBProfiles.Spiral):
    _baseProfile = SBProfiles.Spiral
    _SBkeys = [['A', 'B', 'N', 'amp', 'h', 'omega', 'pa', 'q', 'rmax', 'x', 'y']]

    def __init__(self, name, pars, convolve=0):
        SBModel.__init__(self, name, pars, convolve)

    def getMag(self, amp, zp):
        cnts = amp
        return cnts2mag(cnts, zp)

    def Mag(self, zp):
        return self.getMag(self.amp, zp)


