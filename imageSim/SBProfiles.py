import numpy as np
from scipy import interpolate
"""
MINIMAL ERROR CHECKING!
"""
def cnts2mag(cnts,zp):
    from math import log10
    return -2.5*log10(cnts) + zp


class Sersic:
    def __init__(self,x=None,y=None,q=None,pa=None,re=None,amp=None,n=None):
        self.x = x
        self.y = y
        self.q = q
        self.pa = pa
        self.re = re
        self.amp = amp
        self.n = n
        self.convolve = True

    def eval(self,r):
        k = 2.*self.n-1./3+4./(405.*self.n)+46/(25515.*self.n**2)
        R = r/self.re
        return self.amp*np.exp(-k*(R**(1./self.n) - 1.))

    def coords2rcirc(self,x,y):
        from math import pi,cos as COS,sin as SIN
        shape = x.shape
        x = x.ravel()
        y = y.ravel()

        cos = COS(self.pa*pi/180.)
        sin = SIN(self.pa*pi/180.)
        xp = (x-self.x)*cos+(y-self.y)*sin
        yp = (y-self.y)*cos-(x-self.x)*sin
        r = (self.q*xp**2+yp**2/self.q)**0.5
        return r

    def pixeval(self,x,y,scale=1,csub=23):
        from scipy import interpolate
        from math import pi,cos as COS,sin as SIN

        shape = x.shape
        x = x.ravel()
        y = y.ravel()

        cos = COS(self.pa*pi/180.)
        sin = SIN(self.pa*pi/180.)
        xp = (x-self.x)*cos+(y-self.y)*sin
        yp = (y-self.y)*cos-(x-self.x)*sin
        r = (self.q*xp**2+yp**2/self.q)**0.5

        k = 2.*self.n-1./3+4./(405.*self.n)+46/(25515.*self.n**2)
        R = np.logspace(-5.,4.,451) # 50 pnts / decade
        s0 = np.exp(-k*(R**(1./self.n) - 1.))

        # Determine corrections for curvature
        rpow = R**(1./self.n - 1.)
        term1 = (k*rpow/self.n)**2
        term2 = k*(self.n-1.)*rpow/(R*self.n**2)
        wid = scale/self.re
        corr = (term1+term2)*wid**3/6.
        try:
            minR = R[abs(corr)<0.005].min()
        except:
            minR = 0

        # Evaluate model!
        model = interpolate.splrep(R,s0,k=3,s=0)
        R0 = r/self.re
        s = interpolate.splev(R0,model)*scale**2
        if self.n<=1. or minR==0:
            return self.amp*s.reshape(shape)
        model2 = interpolate.splrep(R,s0*R*self.re**2,k=3,s=0)
        coords = np.where(R0<minR)[0]
        c = (np.indices((csub,csub)).astype(np.float32)-csub/2)*scale/csub
        for i in coords:
            # The central pixels are tricky because we can't assume that we
            #   are integrating in delta-theta segments of an annulus; these
            #   pixels are treated separately by sub-sampling with ~500 pixels
            if R0[i]<3*scale/self.re: # the pixels within 3*scale are evaluated by sub-sampling
                s[i] = 0.
                y0 = c[1]+y[i]
                x0 = c[0]+x[i]
                xp = (x0-self.x)*cos+(y0-self.y)*sin
                yp = (y0-self.y)*cos-(x0-self.x)*sin
                r0 = (self.q*xp**2+yp**2/self.q)**0.5/self.re
                s[i] = interpolate.splev(r0.ravel(),model).mean()*scale**2
                continue
            lo = R0[i]-0.5*scale/self.re
            hi = R0[i]+0.5*scale/self.re
            angle = (scale/self.re)/R0[i]
            s[i] = angle*interpolate.splint(lo,hi,model2)
        return self.amp*s.reshape(shape)


class Ring:
    def __init__(self, x=None, y=None, q=None, pa=None, rr=None, amp=None, hi=None, ho=None):
        self.x = x
        self.y = y
        self.q = q
        self.pa = pa
        self.rr = rr
        self.amp = amp
        self.hi = hi
        self.ho = ho
        self.convolve = True

    def pixeval(self, x, y):
        from math import pi,cos as COS,sin as SIN
        shape = x.shape
        x = x.ravel()
        y = y.ravel()

        cos = COS(self.pa*pi/180.)
        sin = SIN(self.pa*pi/180.)
        xp = (x-self.x)*cos+(y-self.y)*sin
        yp = (y-self.y)*cos-(x-self.x)*sin
        r = (self.q*xp**2+yp**2/self.q)**0.5

        inner = 0.5*(np.sign(self.rr - r) + 1.)*np.exp(-(self.rr - r)/self.hi)
        outer = 0.5*(np.sign(r - self.rr) + 1.)*np.exp(-(r - self.rr)/self.ho)

        return self.amp*(inner + outer).reshape(shape)


