import numpy as np
from scipy.integrate import quad


# assuming flat cosmology

Mpc = 3.08568025e24
c = 2.99792458e10
H0 = 70.
omegaL = 0.7
omegaM = 0.3
yr = 365*3600*24.

def comovd(z): # comoving distance in Mpc
    I = quad(lambda z: 1./(omegaL + omegaM*(1+z)**3)**0.5, 0., z)
    return c/(H0*10.**5)*I[0]

def uniage(z): # age of the Universe in Gyr
    fint = lambda z: 1./(1+z)/(omegaL + omegaM*(1+z)**3)**0.5
    return quad(fint, z, np.inf)[0]/H0/yr*Mpc/10.**5/10.**9

